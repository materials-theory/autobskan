"""Pre-SCF protection against untranslated bSKAN symmetry operations.

bSKAN 3.8 finds a space-group operation ``{R | t}``, but its wavefunction
expansion retains only ``R``.  A nonzero ``t`` can therefore create repeated
or ghost contrast in an STM image.  This module mirrors the relevant basis
matching in ``PGBSYM`` and finds a common origin ``o`` satisfying

    (I - R) o = t  (modulo lattice translations).

When such an origin exists, the structure must be translated before the SCF
calculation so every downstream wavefunction file uses the same origin.
"""

from __future__ import annotations

import itertools
import math
import os
import shlex
import tempfile
from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum
from io import StringIO
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io.vasp import read_vasp, write_vasp
from scipy.spatial import cKDTree

try:
    import spglib
except ImportError as exc:  # pragma: no cover - dependency error is environment-specific
    raise ImportError(
        "AutoBSKAN symmetry pre-check requires spglib. Install autobskan with "
        "its declared dependencies."
    ) from exc


# PGBSYM compares every fractional component against this value and regards a
# translation as zero only when the sum of its absolute components is <= TOL.
BSKAN_SYMMETRY_TOLERANCE = 1.0e-4
BSKAN_SURFACE_NORMAL_TOLERANCE_DEGREES = 1.0e-4


class PrecheckStatus(str, Enum):
    READY = "READY"
    SHIFT_REQUIRED = "SHIFT_REQUIRED"
    UNSAFE = "UNSAFE"
    UNSUPPORTED_GEOMETRY = "UNSUPPORTED_GEOMETRY"


@dataclass(frozen=True)
class SymmetryOperation:
    """One basis-preserving operation in fractional coordinates."""

    rotation: np.ndarray
    translation: np.ndarray


@dataclass(frozen=True)
class PrecheckResult:
    """Result of a bSKAN symmetry-origin check."""

    status: PrecheckStatus
    operations: tuple[SymmetryOperation, ...]
    nonzero_operation_indices: tuple[int, ...]
    origin_fractional: np.ndarray | None
    ase_displacement_angstrom: np.ndarray | None
    max_residual_translation: float
    tolerance: float
    unsafe_reason: str | None = None

    @property
    def operation_count(self) -> int:
        return len(self.operations)

    @property
    def nonzero_operation_count(self) -> int:
        return len(self.nonzero_operation_indices)


class ExistingScfRequiresRestartError(RuntimeError):
    """Raised when existing SCF data cannot safely feed bSKAN."""


@dataclass(frozen=True)
class _LoadedPoscar:
    atoms: Atoms
    type_ids: np.ndarray
    symbol_count: tuple[tuple[str, int], ...]
    selective_flags: tuple[tuple[str, str, str], ...] | None
    title: str


def _absolute(path: os.PathLike[str] | str) -> Path:
    return Path(os.path.abspath(os.path.expanduser(os.fspath(path))))


def _centered_fractional(values: np.ndarray | Sequence[float]) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    return np.mod(array + 0.5, 1.0) - 0.5


def _unit_fractional(values: np.ndarray | Sequence[float]) -> np.ndarray:
    wrapped = np.mod(np.asarray(values, dtype=np.float64), 1.0)
    # np.mod(-epsilon, 1) can round to exactly 1.0. It is the same periodic
    # point as zero, but scipy's periodic KD-tree requires [0, boxsize).
    wrapped[np.isclose(wrapped, 1.0, atol=1.0e-12, rtol=0.0)] = 0.0
    return wrapped


def _integer_tokens(line: str) -> list[int] | None:
    tokens = line.split()
    if not tokens:
        return None
    try:
        values = [int(token) for token in tokens]
    except ValueError:
        return None
    return values if all(value >= 0 for value in values) else None


def _load_poscar(path: os.PathLike[str] | str) -> _LoadedPoscar:
    source = _absolute(path)
    try:
        lines = source.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError as exc:
        raise ValueError(f"Cannot read structure file: {source}") from exc
    if len(lines) < 7:
        raise ValueError(f"Incomplete VASP structure file: {source}")

    atoms = read_vasp(source)
    count_line_index = 5
    counts = _integer_tokens(lines[count_line_index])
    symbols: list[str]
    if counts is None:
        symbols = lines[5].split()
        count_line_index = 6
        counts = _integer_tokens(lines[count_line_index])
        if counts is None or len(symbols) != len(counts):
            raise ValueError(f"Invalid VASP species/count records: {source}")
    else:
        chemical_symbols = atoms.get_chemical_symbols()
        offsets = np.cumsum([0, *counts[:-1]])
        symbols = [chemical_symbols[int(index)] for index in offsets]

    if sum(counts) != len(atoms):
        raise ValueError(
            f"VASP atom counts total {sum(counts)}, but ASE read {len(atoms)} atoms: {source}"
        )

    cursor = count_line_index + 1
    if cursor >= len(lines):
        raise ValueError(f"VASP coordinate mode is missing: {source}")
    has_selective_dynamics = lines[cursor].lstrip().lower().startswith("s")
    if has_selective_dynamics:
        cursor += 1
    if cursor >= len(lines) or not lines[cursor].lstrip().lower().startswith(
        ("d", "c", "k")
    ):
        raise ValueError(f"Invalid VASP coordinate mode: {source}")

    selective_flags: tuple[tuple[str, str, str], ...] | None = None
    if has_selective_dynamics:
        position_start = cursor + 1
        if position_start + len(atoms) > len(lines):
            raise ValueError(f"VASP atom positions are incomplete: {source}")
        parsed_flags: list[tuple[str, str, str]] = []
        for line_number, line in enumerate(
            lines[position_start : position_start + len(atoms)],
            start=position_start + 1,
        ):
            tokens = line.split()
            if len(tokens) < 6:
                raise ValueError(
                    "Selective-dynamics flags are missing at "
                    f"{source}:{line_number}"
                )
            flags = tuple(token[:1].upper() for token in tokens[3:6])
            if any(flag not in {"T", "F"} for flag in flags):
                raise ValueError(
                    "Invalid selective-dynamics flag at "
                    f"{source}:{line_number}; expected T or F"
                )
            parsed_flags.append(flags)
        selective_flags = tuple(parsed_flags)

    type_ids = np.repeat(np.arange(len(counts), dtype=np.int64), counts)
    return _LoadedPoscar(
        atoms=atoms,
        type_ids=type_ids,
        symbol_count=tuple(zip(symbols, counts)),
        selective_flags=selective_flags,
        title=lines[0],
    )


def _unique_lattice_rotations(cell: np.ndarray, lattice_symprec: float) -> list[np.ndarray]:
    symmetry = spglib.get_symmetry(
        (np.asarray(cell, dtype=np.float64), np.zeros((1, 3)), np.ones(1, dtype=int)),
        symprec=lattice_symprec,
    )
    if symmetry is None:
        raise ValueError("spglib could not determine the lattice point group")

    rotations: list[np.ndarray] = []
    for candidate in symmetry["rotations"]:
        rotation = np.asarray(candidate, dtype=np.int64)
        if not any(np.array_equal(rotation, known) for known in rotations):
            rotations.append(rotation)
    if not rotations:
        raise ValueError("No lattice symmetry operation was found")
    return rotations


def _basis_operations(
    atoms: Atoms,
    type_ids: np.ndarray,
    *,
    tolerance: float,
    lattice_symprec: float,
) -> tuple[SymmetryOperation, ...]:
    positions = _unit_fractional(atoms.get_scaled_positions(wrap=False))
    types = np.asarray(type_ids, dtype=np.int64)
    if types.shape != (len(atoms),):
        raise ValueError("type_ids must contain one bSKAN type index per atom")

    type_trees = {
        int(atom_type): cKDTree(positions[types == atom_type], boxsize=1.0)
        for atom_type in np.unique(types)
    }
    first_type_candidates = np.flatnonzero(types == types[0])
    operations: list[SymmetryOperation] = []

    for rotation in _unique_lattice_rotations(atoms.cell.array, lattice_symprec):
        transformed = _unit_fractional(positions @ rotation.T)
        for target_index in first_type_candidates:
            translation = _centered_fractional(
                positions[target_index] - transformed[0]
            )
            accepted = True
            for atom_type, tree in type_trees.items():
                targets = _unit_fractional(
                    transformed[types == atom_type] + translation
                )
                distances, _indices = tree.query(
                    targets,
                    k=1,
                    p=np.inf,
                    distance_upper_bound=tolerance * (1.0 + 1.0e-12),
                )
                if not np.all(np.isfinite(distances)):
                    accepted = False
                    break
            if accepted:
                operations.append(
                    SymmetryOperation(
                        rotation=rotation.copy(),
                        translation=translation.copy(),
                    )
                )
                break

    if not operations:
        raise ValueError("The identity operation did not preserve the atomic basis")
    return tuple(operations)


def _independent_row_indices(matrix: np.ndarray) -> list[int]:
    selected: list[int] = []
    rank = 0
    for index, row in enumerate(matrix):
        trial = matrix[selected + [index]]
        trial_rank = int(np.linalg.matrix_rank(trial, tol=1.0e-10))
        if trial_rank > rank:
            selected.append(index)
            rank = trial_rank
        if rank == 3:
            break
    return selected


def _nearest_integer(values: np.ndarray) -> np.ndarray:
    """Fortran NINT-compatible nearest integers, including half-way values."""

    array = np.asarray(values, dtype=np.float64)
    return np.copysign(np.floor(np.abs(array) + 0.5), array)


def _origin_candidates(
    operations: Sequence[SymmetryOperation],
) -> list[np.ndarray]:
    identity = np.eye(3, dtype=np.float64)
    matrix = np.vstack([identity - operation.rotation for operation in operations])
    translation = np.concatenate([operation.translation for operation in operations])
    selected = _independent_row_indices(matrix)
    if not selected:
        return [np.zeros(3, dtype=np.float64)]

    basis = matrix[selected]
    integer_ranges: list[range] = []
    for row_index in selected:
        row = matrix[row_index]
        lower = float(np.minimum(row, 0.0).sum())
        upper = float(np.maximum(row, 0.0).sum())
        offset = float(translation[row_index])
        first = math.ceil(lower - offset - 1.0e-10)
        last = math.floor(upper - offset + 1.0e-10)
        integer_ranges.append(range(first, last + 1))

    candidates: list[np.ndarray] = []
    selected_translation = translation[selected]
    for integers in itertools.product(*integer_ranges):
        seed = np.linalg.lstsq(
            basis,
            selected_translation + np.asarray(integers, dtype=np.float64),
            rcond=None,
        )[0]
        branch = _nearest_integer(matrix @ seed - translation)
        refined = np.linalg.lstsq(matrix, translation + branch, rcond=None)[0]
        candidates.append(np.mod(refined, 1.0))
    return candidates


def _find_common_origin(
    operations: Sequence[SymmetryOperation],
    cell: np.ndarray,
    tolerance: float,
) -> tuple[np.ndarray | None, float]:
    identity = np.eye(3, dtype=np.float64)
    matrix = np.vstack([identity - operation.rotation for operation in operations])
    translation = np.concatenate([operation.translation for operation in operations])

    valid: list[tuple[tuple[float, ...], np.ndarray, float]] = []
    smallest_residual = math.inf
    for candidate in _origin_candidates(operations):
        residual = _centered_fractional(matrix @ candidate - translation).reshape(-1, 3)
        per_operation = np.sum(np.abs(residual), axis=1)
        maximum = float(np.max(per_operation, initial=0.0))
        smallest_residual = min(smallest_residual, maximum)
        if maximum > tolerance:
            continue

        centered = _centered_fractional(candidate)
        cartesian_norm = float(np.linalg.norm(centered @ np.asarray(cell)))
        # Prefer the shortest physical displacement, then a stable [0, 1)
        # representative. This selects o=(0.25, 0, 0) for the p2s mirror.
        score = (cartesian_norm, maximum, *np.mod(candidate, 1.0).tolist())
        valid.append((score, centered, maximum))

    if not valid:
        return None, smallest_residual
    _score, origin, maximum = min(valid, key=lambda item: item[0])
    return origin, maximum


def _surface_normal_geometry_error(
    cell: np.ndarray,
    tolerance_degrees: float,
) -> str | None:
    lattice = np.asarray(cell, dtype=np.float64)
    a_vector, b_vector, c_vector = lattice
    ab_normal = np.cross(a_vector, b_vector)
    normal_length = float(np.linalg.norm(ab_normal))
    c_length = float(np.linalg.norm(c_vector))
    if normal_length <= np.finfo(float).eps or c_length <= np.finfo(float).eps:
        return (
            "The lattice is degenerate, so the ab-plane normal and c-axis "
            "alignment cannot be determined. The bSKAN pre-check cannot "
            "guarantee a safe symmetry expansion for this geometry."
        )

    normal_unit = ab_normal / normal_length
    normal_projection = float(np.dot(c_vector, normal_unit))
    in_plane_vector = c_vector - normal_projection * normal_unit
    in_plane_length = float(np.linalg.norm(in_plane_vector))
    alignment = float(
        np.clip(abs(normal_projection) / c_length, 0.0, 1.0)
    )
    deviation_degrees = float(np.degrees(np.arccos(alignment)))
    if deviation_degrees <= tolerance_degrees:
        return None

    return (
        "The c lattice vector is not perpendicular to the ab plane "
        f"(deviation from the ab-plane normal: {deviation_degrees:.8f} deg; "
        f"in-plane c component: {in_plane_length:.8f} A). The current bSKAN "
        "pre-check assumes c is parallel to a x b, so it cannot guarantee "
        "that the symmetry-expanded STM result is free of repetition or "
        "ghost-contrast errors."
    )


def precheck_bskan_structure(
    structure: os.PathLike[str] | str | Atoms,
    *,
    type_ids: Sequence[int] | None = None,
    tolerance: float = BSKAN_SYMMETRY_TOLERANCE,
    lattice_symprec: float = 1.0e-5,
    surface_normal_tolerance_degrees: float = (
        BSKAN_SURFACE_NORMAL_TOLERANCE_DEGREES
    ),
) -> PrecheckResult:
    """Classify a structure before it is used for a VASP SCF calculation.

    Paths are parsed with their POSCAR count groups because bSKAN treats each
    count group as a separate atom type. ``Atoms`` inputs use atomic numbers
    unless explicit ``type_ids`` are supplied.
    """

    if (
        tolerance <= 0.0
        or lattice_symprec <= 0.0
        or surface_normal_tolerance_degrees <= 0.0
    ):
        raise ValueError("Symmetry tolerances must be positive")

    if isinstance(structure, Atoms):
        atoms = structure.copy()
        resolved_types = (
            np.asarray(type_ids, dtype=np.int64)
            if type_ids is not None
            else np.asarray(atoms.numbers, dtype=np.int64)
        )
    else:
        loaded = _load_poscar(structure)
        atoms = loaded.atoms
        resolved_types = (
            np.asarray(type_ids, dtype=np.int64)
            if type_ids is not None
            else loaded.type_ids
        )

    geometry_error = _surface_normal_geometry_error(
        atoms.cell.array,
        surface_normal_tolerance_degrees,
    )
    if geometry_error is not None:
        return PrecheckResult(
            status=PrecheckStatus.UNSUPPORTED_GEOMETRY,
            operations=(),
            nonzero_operation_indices=(),
            origin_fractional=None,
            ase_displacement_angstrom=None,
            max_residual_translation=math.nan,
            tolerance=tolerance,
            unsafe_reason=geometry_error,
        )

    operations = _basis_operations(
        atoms,
        resolved_types,
        tolerance=tolerance,
        lattice_symprec=lattice_symprec,
    )
    nonzero = tuple(
        index
        for index, operation in enumerate(operations)
        if float(np.sum(np.abs(operation.translation))) > tolerance
    )
    if not nonzero:
        return PrecheckResult(
            status=PrecheckStatus.READY,
            operations=operations,
            nonzero_operation_indices=(),
            origin_fractional=np.zeros(3, dtype=np.float64),
            ase_displacement_angstrom=np.zeros(3, dtype=np.float64),
            max_residual_translation=0.0,
            tolerance=tolerance,
        )

    origin, residual = _find_common_origin(operations, atoms.cell.array, tolerance)
    if origin is None:
        intrinsic_operations = []
        for index in nonzero:
            individual_origin, _individual_residual = _find_common_origin(
                (operations[index],), atoms.cell.array, tolerance
            )
            if individual_origin is None:
                intrinsic_operations.append(index + 1)
        if intrinsic_operations:
            operation_text = ", ".join(str(index) for index in intrinsic_operations)
            unsafe_reason = (
                "Operation(s) "
                f"{operation_text} contain an intrinsic translation component "
                "that an origin shift cannot remove, as occurs for a glide or "
                "screw operation."
            )
        else:
            unsafe_reason = (
                "The translated operations require mutually incompatible origins; "
                "each can be shifted individually, but no single global origin "
                "makes all translations zero."
            )
        return PrecheckResult(
            status=PrecheckStatus.UNSAFE,
            operations=operations,
            nonzero_operation_indices=nonzero,
            origin_fractional=None,
            ase_displacement_angstrom=None,
            max_residual_translation=residual,
            tolerance=tolerance,
            unsafe_reason=unsafe_reason,
        )

    displacement = -(origin @ np.asarray(atoms.cell.array, dtype=np.float64))
    return PrecheckResult(
        status=PrecheckStatus.SHIFT_REQUIRED,
        operations=operations,
        nonzero_operation_indices=nonzero,
        origin_fractional=origin,
        ase_displacement_angstrom=displacement,
        max_residual_translation=residual,
        tolerance=tolerance,
    )


def apply_bskan_origin_shift(atoms: Atoms, result: PrecheckResult) -> Atoms:
    """Return an ASE-translated copy for a ``SHIFT_REQUIRED`` result."""

    if result.status is not PrecheckStatus.SHIFT_REQUIRED:
        raise ValueError("An origin shift can only be applied to SHIFT_REQUIRED structures")
    if result.ase_displacement_angstrom is None:
        raise ValueError("The pre-check result contains no ASE displacement")

    translated = atoms.copy()
    translated.translate(result.ase_displacement_angstrom)
    translated.wrap(pbc=True, eps=1.0e-12)
    return translated


def write_prechecked_poscar(
    source_path: os.PathLike[str] | str,
    output_path: os.PathLike[str] | str,
    *,
    force: bool = False,
    tolerance: float = BSKAN_SYMMETRY_TOLERANCE,
    lattice_symprec: float = 1.0e-5,
    surface_normal_tolerance_degrees: float = (
        BSKAN_SURFACE_NORMAL_TOLERANCE_DEGREES
    ),
) -> tuple[Path, PrecheckResult, PrecheckResult]:
    """Apply an ASE origin shift, write a POSCAR, and verify it is ``READY``."""

    source = _absolute(source_path)
    output = _absolute(output_path)
    if source == output:
        raise ValueError("Refusing to overwrite the SCF input structure in place")
    if output.exists() and not force:
        raise FileExistsError(f"Output already exists (use --force to replace it): {output}")

    loaded = _load_poscar(source)
    before = precheck_bskan_structure(
        loaded.atoms,
        type_ids=loaded.type_ids,
        tolerance=tolerance,
        lattice_symprec=lattice_symprec,
        surface_normal_tolerance_degrees=surface_normal_tolerance_degrees,
    )
    if before.status is not PrecheckStatus.SHIFT_REQUIRED:
        raise ValueError(f"Structure status is {before.status.value}; no removable shift to apply")

    translated = apply_bskan_origin_shift(loaded.atoms, before)
    after = precheck_bskan_structure(
        translated,
        type_ids=loaded.type_ids,
        tolerance=tolerance,
        lattice_symprec=lattice_symprec,
        surface_normal_tolerance_degrees=surface_normal_tolerance_degrees,
    )
    if after.status is not PrecheckStatus.READY:
        raise RuntimeError(
            "Post-translation verification failed; no structure file was written"
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="\n",
            dir=output.parent,
            prefix=f".{output.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_name = handle.name
            write_vasp(
                handle,
                translated,
                direct=True,
                sort=False,
                symbol_count=list(loaded.symbol_count),
            )
        temporary = Path(temporary_name)
        written_lines = temporary.read_text(encoding="utf-8").splitlines(keepends=True)
        if written_lines:
            written_lines[0] = loaded.title + "\n"
            temporary.write_text("".join(written_lines), encoding="utf-8", newline="\n")
        os.replace(temporary, output)
    except Exception:
        if temporary_name is not None:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass
        raise
    return output, before, after


def _render_bskan_asample(loaded: _LoadedPoscar) -> str:
    """Render one native bSKAN ASAMPLE in normalized VASP 4 form."""

    atoms = loaded.atoms.copy()
    atoms.set_constraint()
    buffer = StringIO()
    write_vasp(
        buffer,
        atoms,
        direct=True,
        sort=False,
        symbol_count=list(loaded.symbol_count),
    )
    lines = buffer.getvalue().splitlines()
    if len(lines) < 8 + len(atoms):
        raise RuntimeError("ASE returned an incomplete VASP structure")

    lines[0] = loaded.title
    del lines[5]  # Native bSKAN ASAMPLE uses VASP 4 counts without symbols.
    if not lines[6].lstrip().lower().startswith("d"):
        raise RuntimeError("ASE did not render Direct coordinates")
    lines.insert(6, "Selective dynamics")

    flags = loaded.selective_flags
    if flags is None:
        flags = tuple(("T", "T", "T") for _ in atoms)
    if len(flags) != len(atoms):
        raise RuntimeError("Selective-dynamics flag count does not match atom count")

    position_start = 8
    for index, atom_flags in enumerate(flags):
        lines[position_start + index] = (
            f"{lines[position_start + index].rstrip()}   "
            + "   ".join(atom_flags)
        )
    return "\n".join(lines) + "\n"


def _validate_bskan_asample_text(
    text: str,
    *,
    expected_counts: Sequence[int],
    expected_flags: Sequence[Sequence[str]],
) -> None:
    lines = text.splitlines()
    atom_count = sum(expected_counts)
    if len(lines) < 8 + atom_count:
        raise RuntimeError("Generated ASAMPLE is incomplete")
    if _integer_tokens(lines[5]) != list(expected_counts):
        raise RuntimeError("Generated ASAMPLE is not VASP 4 count-only format")
    if not lines[6].lstrip().lower().startswith("s"):
        raise RuntimeError("Generated ASAMPLE is missing Selective dynamics")
    if not lines[7].lstrip().lower().startswith("d"):
        raise RuntimeError("Generated ASAMPLE does not use Direct coordinates")

    for index, expected in enumerate(expected_flags):
        tokens = lines[8 + index].split()
        if len(tokens) < 6:
            raise RuntimeError("Generated ASAMPLE position is missing T/F flags")
        try:
            coordinates = np.asarray(tokens[:3], dtype=np.float64)
        except ValueError as exc:
            raise RuntimeError("Generated ASAMPLE contains an invalid coordinate") from exc
        if not np.all(np.isfinite(coordinates)):
            raise RuntimeError("Generated ASAMPLE contains a non-finite coordinate")
        if tuple(tokens[3:6]) != tuple(expected):
            raise RuntimeError("Generated ASAMPLE changed selective-dynamics flags")


def write_bskan_asample(
    source_path: os.PathLike[str] | str,
    output_path: os.PathLike[str] | str | None = None,
    *,
    force: bool = False,
) -> Path:
    """Convert a POSCAR to native bSKAN ``ASAMPLE`` format.

    The output omits the VASP 5 element-symbol record, always contains
    ``Selective dynamics``, and always uses ``Direct`` coordinates. Existing
    selective flags are retained; an unconstrained input receives ``T T T``
    for every atom.
    """

    source = _absolute(source_path)
    output = _absolute(output_path or source.with_name("ASAMPLE"))
    if source == output:
        raise ValueError("Refusing to overwrite the source structure in place")
    if output.exists() and not force:
        raise FileExistsError(f"Output already exists (use force=True to replace it): {output}")

    loaded = _load_poscar(source)
    flags = loaded.selective_flags
    if flags is None:
        flags = tuple(("T", "T", "T") for _ in loaded.atoms)
    counts = tuple(count for _symbol, count in loaded.symbol_count)
    rendered = _render_bskan_asample(loaded)
    _validate_bskan_asample_text(
        rendered,
        expected_counts=counts,
        expected_flags=flags,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="\n",
            dir=output.parent,
            prefix=f".{output.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_name = handle.name
            handle.write(rendered)
        os.replace(temporary_name, output)
    except Exception:
        if temporary_name is not None:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass
        raise
    return output


def _structures_match(left: Atoms, right: Atoms, tolerance: float = 1.0e-5) -> bool:
    if len(left) != len(right):
        return False
    if left.get_chemical_symbols() != right.get_chemical_symbols():
        return False
    if not np.allclose(left.cell.array, right.cell.array, atol=tolerance, rtol=0.0):
        return False
    difference = _centered_fractional(
        left.get_scaled_positions(wrap=False) - right.get_scaled_positions(wrap=False)
    )
    return bool(np.max(np.abs(difference), initial=0.0) <= tolerance)


def validate_existing_scf_for_bskan(
    scf_directory: os.PathLike[str] | str,
) -> PrecheckResult:
    """Reject existing SCF data that need an origin shift or are inconsistent."""

    directory = _absolute(scf_directory)
    poscar = directory / "POSCAR"
    if not poscar.is_file():
        raise ExistingScfRequiresRestartError(f"SCF POSCAR is missing: {poscar}")

    result = precheck_bskan_structure(poscar)
    check_command = shlex.join(["autobskan", "pre-check", str(poscar)])
    if result.status is PrecheckStatus.UNSUPPORTED_GEOMETRY:
        raise ExistingScfRequiresRestartError(
            "Existing SCF geometry is outside the validated bSKAN pre-check "
            "domain. The c lattice vector must be perpendicular to the ab "
            "plane, but this condition is not satisfied, so symmetry safety "
            "cannot be guaranteed. Run:\n"
            f"  {check_command}\n"
            f"Reported reason: {result.unsafe_reason}\n"
            "Construct an equivalent surface cell with c parallel to a x b, "
            "run pre-check again, and start a new SCF calculation. Do not reuse "
            "the current WAVECAR, CHGCAR, or WAVSAMPLE."
        )
    if result.status is PrecheckStatus.SHIFT_REQUIRED:
        shifted = poscar.with_name(f"{poscar.name}-autobskan_shifted.vasp")
        raise ExistingScfRequiresRestartError(
            "Existing SCF data use an origin that is unsafe for bSKAN: "
            f"{result.nonzero_operation_count} symmetry operation(s) have a "
            "nonzero translation. Do not translate POSCAR/ASAMPLE alone or reuse "
            "the current WAVECAR, CHGCAR, or WAVSAMPLE. Run:\n"
            f"  {check_command}\n"
            f"This automatically creates:\n  {shifted}\n"
            "Use the generated structure as the POSCAR of a new SCF "
            "calculation and restart AutoBSKAN from that SCF directory."
        )
    if result.status is PrecheckStatus.UNSAFE:
        raise ExistingScfRequiresRestartError(
            "Existing SCF data contain bSKAN symmetry translations that no single "
            "origin shift can remove (for example a glide or screw component). "
            "Do not reuse this SCF for bSKAN. Run:\n"
            f"  {check_command}\n"
            "Then rebuild the calculation with a symmetry-free/full-BZ workflow "
            "or a corrected bSKAN symmetry implementation."
        )

    chgcar = directory / "CHGCAR"
    if chgcar.is_file():
        try:
            chgcar_atoms = read_vasp(chgcar)
        except Exception as exc:
            raise ExistingScfRequiresRestartError(
                f"Cannot validate the structure header in SCF CHGCAR: {chgcar}"
            ) from exc
        poscar_atoms = read_vasp(poscar)
        if not _structures_match(poscar_atoms, chgcar_atoms):
            raise ExistingScfRequiresRestartError(
                "SCF POSCAR and CHGCAR structures do not match. This commonly "
                "means POSCAR was shifted after SCF while old wavefunction data "
                "were retained. Run SCF again from the pre-checked POSCAR before "
                "starting AutoBSKAN."
            )
    return result


__all__ = [
    "BSKAN_SURFACE_NORMAL_TOLERANCE_DEGREES",
    "BSKAN_SYMMETRY_TOLERANCE",
    "ExistingScfRequiresRestartError",
    "PrecheckResult",
    "PrecheckStatus",
    "SymmetryOperation",
    "apply_bskan_origin_shift",
    "precheck_bskan_structure",
    "validate_existing_scf_for_bskan",
    "write_bskan_asample",
    "write_prechecked_poscar",
]

"""Create bSKAN ``CURRENT`` files from Chen ``CURSAVE`` output.

The implementation mirrors the formatting and coordinate conversion in
``evaluate.F`` from bSKAN 3.8.  ``ASAMPLE`` is the only required explicit
argument; calculation metadata are discovered beside it by default.
"""

from __future__ import annotations

import math
import os
import re
import tempfile
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# VASP writes WAVSAMPLE coordinates in bohr with this conversion, while
# bSKAN 3.8 converts its single-precision values back with AUTOA below.
VASP_AUTOA = 0.529177249
BSKAN_AUTOA = np.float32(0.529177)
BSKAN_DELZ_BOHR = np.float32(0.1)
_NUMBER = re.compile(
    r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?$"
)


@dataclass(frozen=True)
class _Structure:
    title: str
    lattice_angstrom: np.ndarray
    counts: tuple[int, ...]
    direct_positions: np.ndarray


@dataclass(frozen=True)
class _CurSave:
    points_angstrom: np.ndarray
    values: np.ndarray

    @property
    def nz(self) -> int:
        return int(self.values.shape[1])


@dataclass(frozen=True)
class _ScanGeometry:
    xl_bohr: np.ndarray
    yl_bohr: np.ndarray
    ra1_bohr: np.float32
    ra2_bohr: np.float32
    wave_lattice_bohr: np.ndarray


def _absolute_without_resolving(path: os.PathLike[str] | str) -> Path:
    """Return an absolute path while preserving a calculation-local symlink."""

    return Path(os.path.abspath(os.path.expanduser(os.fspath(path))))


def _numeric_tokens(line: str) -> list[float] | None:
    tokens = line.split()
    if not tokens or not all(_NUMBER.fullmatch(token) for token in tokens):
        return None
    return [float(token.replace("D", "E").replace("d", "e")) for token in tokens]


def _parse_title(line: str) -> str:
    stripped = line.strip()
    if not stripped:
        return ""
    if stripped[0] in "'\"" and stripped[-1:] == stripped[0]:
        return stripped[1:-1][:40]
    # bSKAN uses list-directed CHARACTER input: blanks and commas separate
    # unquoted values, and a slash terminates the input record.
    return re.split(r"[\s,/]", stripped, maxsplit=1)[0][:40]


def _parse_structure(path: Path) -> _Structure:
    try:
        lines = path.read_text(encoding="ascii", errors="strict").splitlines()
    except OSError as exc:
        raise ValueError(f"Cannot read ASAMPLE: {path}") from exc

    if len(lines) < 8:
        raise ValueError(f"ASAMPLE is incomplete: {path}")

    try:
        scale_tokens = lines[1].split()
        if len(scale_tokens) != 1:
            raise ValueError("only a scalar ASAMPLE scale is supported")
        scale = float(scale_tokens[0].replace("D", "E").replace("d", "e"))
        if scale <= 0.0:
            raise ValueError("ASAMPLE scale must be positive")
        lattice = np.asarray(
            [[float(value.replace("D", "E").replace("d", "e")) for value in lines[i].split()[:3]]
             for i in range(2, 5)],
            dtype=np.float64,
        )
        if lattice.shape != (3, 3):
            raise ValueError("ASAMPLE must contain three 3-vector lattice rows")
    except (IndexError, ValueError) as exc:
        raise ValueError(f"Invalid ASAMPLE lattice in {path}: {exc}") from exc

    cursor = 5
    count_tokens = lines[cursor].split()
    if not count_tokens or not all(re.fullmatch(r"[+-]?\d+", token) for token in count_tokens):
        # Accept VASP 5 POSCAR input as a convenience; native ASAMPLE omits symbols.
        cursor += 1
        if cursor >= len(lines):
            raise ValueError(f"ASAMPLE atom counts are missing: {path}")
        count_tokens = lines[cursor].split()
    if not count_tokens or not all(re.fullmatch(r"[+-]?\d+", token) for token in count_tokens):
        raise ValueError(f"Invalid ASAMPLE atom counts: {path}")

    counts = tuple(int(token) for token in count_tokens)
    if any(count < 0 for count in counts) or sum(counts) <= 0:
        raise ValueError(f"ASAMPLE atom counts must be positive: {path}")

    cursor += 1
    if cursor < len(lines) and lines[cursor].lstrip().lower().startswith("s"):
        cursor += 1
    if cursor >= len(lines) or not lines[cursor].lstrip().lower().startswith("d"):
        raise ValueError("bSKAN CURRENT conversion requires Direct coordinates")
    cursor += 1

    natoms = sum(counts)
    if cursor + natoms > len(lines):
        raise ValueError(f"ASAMPLE contains fewer than {natoms} atom positions: {path}")
    try:
        positions = np.asarray(
            [[float(value.replace("D", "E").replace("d", "e")) for value in lines[i].split()[:3]]
             for i in range(cursor, cursor + natoms)],
            dtype=np.float64,
        )
    except ValueError as exc:
        raise ValueError(f"Invalid ASAMPLE atom position in {path}") from exc
    if positions.shape != (natoms, 3):
        raise ValueError(f"Each ASAMPLE atom position must have three values: {path}")

    return _Structure(
        title=_parse_title(lines[0]),
        lattice_angstrom=lattice * scale,
        counts=counts,
        direct_positions=positions,
    )


def _parse_cursave(path: Path) -> _CurSave:
    try:
        lines = path.read_text(encoding="ascii", errors="strict").splitlines()
    except OSError as exc:
        raise ValueError(f"Cannot read CURSAVE: {path}") from exc

    points: list[list[float]] = []
    blocks: list[list[float]] = []
    active: list[float] | None = None

    for line_number, line in enumerate(lines, start=1):
        values = _numeric_tokens(line)
        if values is None:
            # The bSKAN file starts with the text line "Backup file".
            if line.strip() and (points or line_number != 1):
                raise ValueError(f"Unexpected CURSAVE text at line {line_number}: {line.strip()!r}")
            continue
        if len(values) == 2:
            if active is not None and not active:
                raise ValueError(f"CURSAVE block before line {line_number} has no current values")
            points.append(values)
            active = []
            blocks.append(active)
        elif len(values) == 1 and active is not None:
            active.append(values[0])
        else:
            raise ValueError(
                f"Invalid CURSAVE numeric record at line {line_number}; expected x y or one current value"
            )

    if not points or active is None:
        raise ValueError(f"CURSAVE contains no scan points: {path}")
    if not active:
        raise ValueError(f"The final CURSAVE block has no current values: {path}")
    lengths = {len(block) for block in blocks}
    if len(lengths) != 1:
        raise ValueError(f"CURSAVE scan blocks have inconsistent z-grid lengths: {sorted(lengths)}")

    return _CurSave(
        points_angstrom=np.asarray(points, dtype=np.float64),
        # bSKAN stores and writes these values as REAL.
        values=np.asarray(blocks, dtype=np.float32),
    )


def _parse_inscan(path: Path | None) -> dict[str, tuple[float, ...]]:
    if path is None or not path.is_file():
        return {}
    options: dict[str, tuple[float, ...]] = {}
    for raw_line in path.read_text(encoding="ascii", errors="replace").splitlines():
        line = raw_line.split("#", 1)[0]
        if "=" not in line:
            continue
        key, raw_values = line.split("=", 1)
        values = _numeric_tokens(raw_values)
        if values is not None:
            options[key.strip().upper()] = tuple(values)
    return options


def _f32_dot(left: Sequence[float], right: Sequence[float]) -> np.float32:
    result = np.float32(0.0)
    for lhs, rhs in zip(left, right):
        result = np.float32(result + np.float32(np.float32(lhs) * np.float32(rhs)))
    return result


def _scan_geometry(lattice_angstrom: np.ndarray, cell_repeat: Sequence[float]) -> _ScanGeometry:
    if len(cell_repeat) != 2 or any(value <= 0.0 for value in cell_repeat):
        raise ValueError("CELL must contain two positive repeat factors")

    wave = np.asarray(lattice_angstrom / VASP_AUTOA, dtype=np.float32)
    a1 = wave[0, :2].copy()
    a2 = wave[1, :2].copy()
    cell = np.asarray(cell_repeat, dtype=np.float32)
    tolerance = np.float32(1.0e-3)

    if abs(float(_f32_dot(a1, a2))) <= float(tolerance):
        xl = np.asarray(a1 * cell[0], dtype=np.float32)
        yl = np.asarray(a2 * cell[1], dtype=np.float32)
    else:
        xl = np.asarray((a1 + a2) * cell[0], dtype=np.float32)
        yl = np.asarray((a1 - a2) * cell[1], dtype=np.float32)
        cross_z = np.float32(
            np.float32(xl[0] * yl[1]) - np.float32(xl[1] * yl[0])
        )
        if cross_z < 0.0:
            xl, yl = yl.copy(), xl.copy()

        if abs(float(_f32_dot(xl, yl))) > float(tolerance):
            xl = np.asarray(a1 * cell[0], dtype=np.float32)
            projection = np.float32(
                np.sqrt(_f32_dot(a2, a2))
                * np.sqrt(
                    np.float32(1.0)
                    - np.float32(
                        _f32_dot(a1, a2)
                        / np.sqrt(_f32_dot(a1, a1))
                        / np.sqrt(_f32_dot(a2, a2))
                    ) ** np.float32(2.0)
                )
                * cell[1]
            )
            yl = np.asarray([np.float32(0.0), projection], dtype=np.float32)

    if abs(float(_f32_dot(xl, yl))) > float(tolerance):
        raise ValueError("bSKAN cannot construct a perpendicular scan area from this lattice")

    ra1 = np.float32(np.sqrt(_f32_dot(xl, xl)))
    ra2 = np.float32(np.sqrt(_f32_dot(yl, yl)))
    return _ScanGeometry(xl, yl, ra1, ra2, wave)


def _factor_pairs(number: int) -> Iterable[tuple[int, int]]:
    for divisor in range(1, math.isqrt(number) + 1):
        if number % divisor:
            continue
        quotient = number // divisor
        yield divisor, quotient
        if divisor != quotient:
            yield quotient, divisor


def _infer_xy_shape(points: np.ndarray) -> tuple[int, int]:
    """Infer bSKAN's ix-outer/iy-inner raster dimensions from CURSAVE."""

    npoints = int(points.shape[0])
    candidates: list[tuple[float, int, int]] = []
    for nx, ny in _factor_pairs(npoints):
        if nx > 201 or ny > 201 or (nx == 1 and ny == 1):
            continue
        ix = np.repeat(np.arange(nx, dtype=np.float64), ny)
        iy = np.tile(np.arange(ny, dtype=np.float64), nx)
        design = np.column_stack((np.ones(npoints), ix, iy))
        fitted = design @ np.linalg.lstsq(design, points, rcond=None)[0]
        rms = float(np.sqrt(np.mean(np.square(points - fitted))))
        candidates.append((rms, nx, ny))

    if not candidates:
        raise ValueError(f"Cannot factor {npoints} CURSAVE points into a bSKAN grid")
    rms, nx, ny = min(candidates)
    if rms > 2.0e-3:
        raise ValueError(
            f"CURSAVE coordinates do not form a regular bSKAN raster (best RMS error {rms:.3g} A)"
        )
    return nx, ny


def _pivot_origin_bohr(
    options: dict[str, tuple[float, ...]],
    geometry: _ScanGeometry,
    first_point_angstrom: np.ndarray,
) -> np.ndarray:
    pivot = options.get("PIVOT")
    if pivot is not None and len(pivot) >= 2:
        p = np.asarray(pivot[:2], dtype=np.float32)
        origin = np.zeros(2, dtype=np.float32)
        for component in range(2):
            origin[component] = np.float32(
                np.float32(geometry.wave_lattice_bohr[0, component] * p[0])
                + np.float32(geometry.wave_lattice_bohr[1, component] * p[1])
            )
        rendered = np.asarray(origin * BSKAN_AUTOA, dtype=np.float64)
        if np.linalg.norm(rendered - first_point_angstrom) <= 1.0e-3:
            return np.asarray(origin, dtype=np.float64)

    # CURSAVE is rounded, but this fallback is exact for AutoBSKAN's PIVOT=0.
    return np.asarray(first_point_angstrom / float(BSKAN_AUTOA), dtype=np.float64)


def _fallback_z_vacuum(structure: _Structure) -> float:
    cartesian = structure.direct_positions @ structure.lattice_angstrom
    return float(np.max(cartesian[:, 2]) + 0.5)


def _format_header(
    structure: _Structure,
    geometry: _ScanGeometry,
    grid_shape: tuple[int, int, int],
    z_vacuum_angstrom: float,
    origin_bohr: Sequence[float],
    delz_bohr: float,
) -> str:
    nx, ny, nz = grid_shape
    if min(nx, ny, nz) <= 0:
        raise ValueError("CURRENT grid dimensions must be positive")

    autoa_double = float(BSKAN_AUTOA)
    lattice_bohr = structure.lattice_angstrom.T / autoa_double
    z_vacuum_bohr = float(
        np.float32(np.float32(z_vacuum_angstrom) / BSKAN_AUTOA)
    )
    delz = np.float32(delz_bohr)
    map_matrix = np.zeros((3, 3), dtype=np.float64)
    map_matrix[:2, 0] = np.asarray(geometry.xl_bohr, dtype=np.float64)
    map_matrix[:2, 1] = np.asarray(geometry.yl_bohr, dtype=np.float64)
    map_matrix[2, 2] = float(np.float32(np.float32(nz) * delz))
    inverse_map = np.linalg.inv(map_matrix)

    atom_types = np.repeat(np.arange(len(structure.counts)), structure.counts)
    converted: list[np.ndarray] = []
    kept_counts = [0] * len(structure.counts)
    origin = np.asarray(origin_bohr, dtype=np.float64)
    for atom_type, direct in zip(atom_types, structure.direct_positions):
        cartesian = lattice_bohr @ direct
        relative = cartesian - np.asarray(
            [origin[0], origin[1], z_vacuum_bohr], dtype=np.float64
        )
        mapped = inverse_map @ relative
        if mapped[2] < 0.0:
            converted.append(mapped)
            kept_counts[int(atom_type)] += 1

    if not converted:
        raise ValueError("No ASAMPLE atoms lie below ZVACUUM")

    x_length = np.float32(geometry.ra1_bohr * BSKAN_AUTOA)
    y_length = np.float32(geometry.ra2_bohr * BSKAN_AUTOA)
    z_length = np.float32(np.float32(np.float32(nz) * delz) * BSKAN_AUTOA)

    lines = [
        f" {structure.title:<40}\n",
        "   1.0   \n",
        f"{float(x_length):12.6f}{0.0:12.6f}{0.0:12.6f}\n",
        f"{0.0:12.6f}{float(y_length):12.6f}{0.0:12.6f}\n",
        f"{0.0:12.6f}{0.0:12.6f}{float(z_length):12.6f}\n",
        "".join(f"{count:5d}" for count in kept_counts) + "\n",
        " Direct\n",
    ]
    lines.extend(
        "".join(f"{float(value):12.6f}" for value in coordinate) + "\n"
        for coordinate in converted
    )
    lines.extend(("  \n", f"{nx:5d}{ny:5d}{nz:5d}\n"))
    return "".join(lines)


def build_current_head_from_asample(
    asample_path: os.PathLike[str] | str,
    grid_shape: Sequence[int],
    *,
    inscan_path: os.PathLike[str] | str | None = None,
    cell_repeat: Sequence[float] | None = None,
    z_vacuum_angstrom: float | None = None,
    scan_origin_bohr: Sequence[float] = (0.0, 0.0),
    delz_bohr: float = float(BSKAN_DELZ_BOHR),
) -> str:
    """Return a source-compatible bSKAN ``CURRENT`` header.

    ``grid_shape`` is ``(nx, ny, nz)``.  When available, ``CELL`` and
    ``ZVACUUM`` are read from a sibling ``INSCAN`` file.
    """

    asample = _absolute_without_resolving(asample_path)
    structure = _parse_structure(asample)
    if len(grid_shape) != 3:
        raise ValueError("grid_shape must contain nx, ny, and nz")
    shape = tuple(int(value) for value in grid_shape)

    if inscan_path is None:
        candidate = asample.parent / "INSCAN"
        inscan = candidate if candidate.is_file() else None
    else:
        inscan = _absolute_without_resolving(inscan_path)
    options = _parse_inscan(inscan)

    if cell_repeat is None:
        cell_repeat = options.get("CELL", (1.0, 1.0))
    geometry = _scan_geometry(structure.lattice_angstrom, cell_repeat)
    if z_vacuum_angstrom is None:
        raw_z = options.get("ZVACUUM")
        z_vacuum_angstrom = raw_z[0] if raw_z else _fallback_z_vacuum(structure)

    return _format_header(
        structure,
        geometry,
        shape,
        z_vacuum_angstrom,
        scan_origin_bohr,
        delz_bohr,
    )


def _format_fortran_e14_4(value: float) -> str:
    real_value = float(np.float32(value))
    if not math.isfinite(real_value):
        raise ValueError(f"Non-finite CURSAVE value: {value}")
    if real_value == 0.0:
        core = "0.0000E+00"
    else:
        scientific = f"{real_value:.3E}"
        mantissa_text, exponent_text = scientific.split("E")
        mantissa = float(mantissa_text) / 10.0
        exponent = int(exponent_text) + 1
        core = f"{mantissa:.4f}E{exponent:+03d}"
    if len(core) > 14:
        raise ValueError(f"Value {value} does not fit bSKAN E14.4 format")
    return core.rjust(14)


def _format_current_body(values: np.ndarray, nx: int, ny: int, nz: int) -> str:
    if values.shape != (nx * ny, nz):
        raise ValueError(
            f"CURSAVE has shape {values.shape}, expected ({nx * ny}, {nz})"
        )
    ordered = values.reshape(nx, ny, nz).transpose(2, 1, 0).ravel()
    chunks: list[str] = []
    for index, value in enumerate(ordered, start=1):
        chunks.append(_format_fortran_e14_4(float(value)))
        if index % 5 == 0:
            chunks.append("\n")
    if len(ordered) % 5:
        chunks.append("\n")
    return "".join(chunks)


def asample_to_current(
    asample_path: os.PathLike[str] | str = "ASAMPLE",
    *,
    cursave_path: os.PathLike[str] | str | None = None,
    inscan_path: os.PathLike[str] | str | None = None,
    current_path: os.PathLike[str] | str | None = None,
    cell_repeat: Sequence[float] | None = None,
    z_vacuum_angstrom: float | None = None,
    delz_bohr: float = float(BSKAN_DELZ_BOHR),
) -> Path:
    """Convert a Chen ``CURSAVE`` into a VASP-style bSKAN ``CURRENT``.

    The normal call is simply ``asample_to_current("ASAMPLE")`` from a
    completed Chen bias directory.  ``CURSAVE`` and ``INSCAN`` are then read
    from that directory and ``CURRENT`` is written there atomically.
    """

    asample = _absolute_without_resolving(asample_path)
    cursave = (
        asample.parent / "CURSAVE"
        if cursave_path is None
        else _absolute_without_resolving(cursave_path)
    )
    inscan = (
        asample.parent / "INSCAN"
        if inscan_path is None
        else _absolute_without_resolving(inscan_path)
    )
    output = (
        asample.parent / "CURRENT"
        if current_path is None
        else _absolute_without_resolving(current_path)
    )

    structure = _parse_structure(asample)
    saved = _parse_cursave(cursave)
    nx, ny = _infer_xy_shape(saved.points_angstrom)
    options = _parse_inscan(inscan)
    repeats = cell_repeat if cell_repeat is not None else options.get("CELL", (1.0, 1.0))
    geometry = _scan_geometry(structure.lattice_angstrom, repeats)
    origin = _pivot_origin_bohr(options, geometry, saved.points_angstrom[0])
    raw_z = options.get("ZVACUUM")
    z_vacuum = (
        z_vacuum_angstrom
        if z_vacuum_angstrom is not None
        else (raw_z[0] if raw_z else _fallback_z_vacuum(structure))
    )

    header = _format_header(
        structure,
        geometry,
        (nx, ny, saved.nz),
        float(z_vacuum),
        origin,
        delz_bohr,
    )
    body = _format_current_body(saved.values, nx, ny, saved.nz)

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="ascii",
            newline="\n",
            dir=output.parent,
            prefix=f".{output.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary_name = handle.name
            handle.write(header)
            handle.write(body)
        os.replace(temporary_name, output)
    except Exception:
        if temporary_name is not None:
            try:
                os.unlink(temporary_name)
            except FileNotFoundError:
                pass
        raise
    return output


__all__ = ["asample_to_current", "build_current_head_from_asample"]

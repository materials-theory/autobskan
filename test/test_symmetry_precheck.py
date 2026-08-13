from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.cell import Cell
from ase.constraints import FixAtoms, FixScaled
from ase.io.vasp import write_vasp

from autobskan.calculation import connect
from autobskan.calculation.symmetry_precheck import (
    ExistingScfRequiresRestartError,
    PrecheckStatus,
    apply_bskan_origin_shift,
    precheck_bskan_structure,
    validate_existing_scf_for_bskan,
    write_bskan_asample,
    write_prechecked_poscar,
)
from autobskan.main import main as cli_main


def _shifted_mirror_structure() -> Atoms:
    return Atoms(
        "CuCuOO",
        scaled_positions=[
            [0.10, 0.17, 0.31],
            [0.40, 0.17, 0.31],
            [0.13, 0.37, 0.43],
            [0.37, 0.37, 0.43],
        ],
        cell=np.diag([4.0, 5.0, 20.0]),
        pbc=True,
    )


def _write_poscar(path: Path, atoms: Atoms, symbol_count=None) -> None:
    write_vasp(
        path,
        atoms,
        direct=True,
        sort=False,
        symbol_count=symbol_count,
    )


def test_p2s_style_mirror_shift_is_found_and_removed_with_ase(monkeypatch):
    atoms = _shifted_mirror_structure()
    before = precheck_bskan_structure(atoms)

    assert before.status is PrecheckStatus.SHIFT_REQUIRED
    np.testing.assert_allclose(before.origin_fractional, [0.25, 0.0, 0.0])
    np.testing.assert_allclose(before.ase_displacement_angstrom, [-1.0, 0.0, 0.0])

    calls = []
    ase_translate = Atoms.translate

    def record_translate(self, displacement):
        calls.append(np.asarray(displacement, dtype=float).copy())
        return ase_translate(self, displacement)

    monkeypatch.setattr(Atoms, "translate", record_translate)
    shifted = apply_bskan_origin_shift(atoms, before)
    after = precheck_bskan_structure(shifted)

    assert len(calls) == 1
    np.testing.assert_allclose(calls[0], before.ase_displacement_angstrom)
    assert after.status is PrecheckStatus.READY
    assert after.nonzero_operation_count == 0
    np.testing.assert_allclose(shifted.cell.array, atoms.cell.array)
    assert shifted.get_chemical_symbols() == atoms.get_chemical_symbols()


@pytest.mark.parametrize(
    "cell",
    [
        Cell.fromcellpar([4, 5, 7, 90, 90, 73]).array,
        Cell.fromcellpar([4, 5, 7, 90, 90, 104]).array,
        np.diag([4.0, 5.0, 7.0]),
        np.diag([4.0, 4.0, 7.0]),
        Cell.fromcellpar([4, 4, 7, 90, 90, 120]).array,
        np.diag([4.0, 4.0, 4.0]),
    ],
    ids=["oblique-73", "oblique-104", "orthorhombic", "tetragonal", "hexagonal", "cubic"],
)
def test_common_origin_solver_covers_supported_ab_plane_geometries(cell):
    atoms = Atoms(
        "Cu",
        scaled_positions=[[0.137, 0.219, 0.311]],
        cell=cell,
        pbc=True,
    )

    before = precheck_bskan_structure(atoms)
    assert before.status is PrecheckStatus.SHIFT_REQUIRED

    after = precheck_bskan_structure(apply_bskan_origin_shift(atoms, before))
    assert after.status is PrecheckStatus.READY
    assert after.nonzero_operation_count == 0


@pytest.mark.parametrize(
    "cell",
    [
        Cell.fromcellpar([4, 5, 7, 78, 82, 73]).array,
        Cell.fromcellpar([4, 5, 7, 90, 104, 90]).array,
    ],
    ids=["triclinic", "tilted-c-monoclinic"],
)
def test_tilted_c_axis_is_reported_as_unsupported_geometry(cell):
    atoms = Atoms(
        "Cu",
        scaled_positions=[[0.137, 0.219, 0.311]],
        cell=cell,
        pbc=True,
    )

    result = precheck_bskan_structure(atoms)

    assert result.status is PrecheckStatus.UNSUPPORTED_GEOMETRY
    assert result.operation_count == 0
    assert result.origin_fractional is None
    assert result.ase_displacement_angstrom is None
    assert "not perpendicular to the ab plane" in result.unsafe_reason
    assert "cannot guarantee" in result.unsafe_reason


@pytest.mark.parametrize(
    "scaled_positions",
    [
        # Mirror normal to x plus a y/2 glide.
        [[0.10, 0.10, 0.30], [0.90, 0.60, 0.30], [0.23, 0.31, 0.41], [0.77, 0.81, 0.41]],
        # Two-fold rotation around z plus a z/2 screw component.
        [[0.10, 0.20, 0.10], [0.90, 0.80, 0.60], [0.23, 0.31, 0.21], [0.77, 0.69, 0.71]],
    ],
    ids=["glide", "screw"],
)
def test_intrinsic_translation_is_not_silently_shifted(scaled_positions):
    atoms = Atoms(
        "CuCuOO",
        scaled_positions=scaled_positions,
        cell=np.diag([4.0, 5.0, 20.0]),
        pbc=True,
    )

    result = precheck_bskan_structure(atoms)

    assert result.status is PrecheckStatus.UNSAFE
    assert result.origin_fractional is None
    assert result.ase_displacement_angstrom is None
    with pytest.raises(ValueError, match="SHIFT_REQUIRED"):
        apply_bskan_origin_shift(atoms, result)


def test_poscar_count_groups_match_bskan_atom_types(tmp_path):
    atoms = Atoms(
        "CuCu",
        scaled_positions=[[0.10, 0.20, 0.30], [0.37, 0.43, 0.47]],
        cell=np.diag([4.0, 5.0, 20.0]),
        pbc=True,
    )
    poscar = tmp_path / "POSCAR"
    _write_poscar(poscar, atoms, symbol_count=[("Cu", 1), ("Cu", 1)])

    # bSKAN reads the two count groups as distinct types, even though ASE's
    # atomic numbers are equal. The mirror therefore must not be admitted.
    assert precheck_bskan_structure(poscar).status is PrecheckStatus.READY
    assert precheck_bskan_structure(atoms).status is PrecheckStatus.SHIFT_REQUIRED


def test_asample_writer_uses_vasp4_selective_direct_and_default_ttt(tmp_path):
    atoms = Atoms(
        "CuCuO",
        scaled_positions=[
            [0.13, 0.21, 0.31],
            [0.47, 0.56, 0.43],
            [0.78, 0.34, 0.67],
        ],
        cell=Cell.fromcellpar([4.0, 5.0, 8.0, 78.0, 83.0, 72.0]),
        pbc=True,
    )
    source = tmp_path / "POSCAR"
    write_vasp(
        source,
        atoms,
        direct=False,
        sort=False,
        symbol_count=[("Cu", 2), ("O", 1)],
    )

    output = write_bskan_asample(source)
    lines = output.read_text(encoding="utf-8").splitlines()

    assert output == tmp_path / "ASAMPLE"
    assert lines[5].split() == ["2", "1"]
    assert lines[6] == "Selective dynamics"
    assert lines[7] == "Direct"
    coordinates = np.asarray(
        [[float(value) for value in line.split()[:3]] for line in lines[8:11]]
    )
    np.testing.assert_allclose(coordinates, atoms.get_scaled_positions(wrap=False))
    assert [line.split()[3:6] for line in lines[8:11]] == [
        ["T", "T", "T"],
        ["T", "T", "T"],
        ["T", "T", "T"],
    ]


def test_asample_writer_preserves_existing_selective_flags(tmp_path):
    atoms = Atoms(
        "CuCuO",
        scaled_positions=[
            [0.11, 0.22, 0.33],
            [0.42, 0.51, 0.62],
            [0.73, 0.36, 0.81],
        ],
        cell=np.diag([4.0, 5.0, 20.0]),
        pbc=True,
    )
    atoms.set_constraint(
        [
            FixScaled(0, mask=(True, False, True)),
            FixAtoms(indices=[1]),
        ]
    )
    source = tmp_path / "POSCAR"
    _write_poscar(source, atoms, symbol_count=[("Cu", 2), ("O", 1)])

    output = write_bskan_asample(source)
    lines = output.read_text(encoding="utf-8").splitlines()

    assert [line.split()[3:6] for line in lines[8:11]] == [
        ["F", "T", "F"],
        ["F", "F", "F"],
        ["T", "T", "T"],
    ]


def test_calculation_asample_converter_uses_checked_writer(tmp_path, monkeypatch):
    source = tmp_path / "POSCAR"
    atoms = Atoms(
        "Cu", scaled_positions=[[0.0, 0.0, 0.0]], cell=[4, 5, 20], pbc=True
    )
    write_vasp(source, atoms, direct=False, sort=False)
    monkeypatch.chdir(tmp_path)

    output = connect.contcar_to_asample(source)

    assert output == tmp_path / "ASAMPLE"
    lines = output.read_text(encoding="utf-8").splitlines()
    assert lines[5].split() == ["1"]
    assert lines[6:8] == ["Selective dynamics", "Direct"]
    assert lines[8].split()[3:6] == ["T", "T", "T"]


@pytest.mark.skip(reason="CLI integration is delivered in the follow-up application PR")
def test_cli_precheck_automatically_writes_verified_structure(tmp_path, capsys):
    source = tmp_path / "POSCAR"
    output = tmp_path / "POSCAR-autobskan_shifted.vasp"
    asample = tmp_path / "ASAMPLE"
    _write_poscar(source, _shifted_mirror_structure())

    exit_code = cli_main(["pre-check", str(source)])
    output_text = capsys.readouterr().out

    assert exit_code == 0
    assert "[1/4] Reading VASP structure" in output_text
    assert "Status: SHIFT_REQUIRED" in output_text
    assert "ASE Atoms.translate()" in output_text
    assert "Shift applied and verified: created" in output_text
    assert f"Generated file: {output}" in output_text
    assert f"No ASAMPLE was generated or modified: {asample}" in output_text
    assert output.is_file()
    assert not asample.exists()
    assert precheck_bskan_structure(output).status is PrecheckStatus.READY

    second_exit = cli_main(["pre-check", str(source)])
    second_output = capsys.readouterr().out
    assert second_exit == 0
    assert "Shift applied and verified: updated" in second_output


@pytest.mark.skip(reason="CLI integration is delivered in the follow-up application PR")
def test_cli_ready_does_not_generate_a_shifted_file(tmp_path, capsys):
    source = tmp_path / "POSCAR"
    shifted = tmp_path / "POSCAR-autobskan_shifted.vasp"
    asample = tmp_path / "ASAMPLE"
    atoms = Atoms(
        "Cu", scaled_positions=[[0.0, 0.0, 0.0]], cell=[4, 5, 20], pbc=True
    )
    _write_poscar(source, atoms)

    exit_code = cli_main(["pre-check", str(source), "--verbose"])
    output_text = capsys.readouterr().out

    assert exit_code == 0
    assert "Status: READY" in output_text
    assert "Shift not required" in output_text
    assert f"No shifted file was generated: {shifted}" in output_text
    assert f"Generated ASAMPLE: {asample}" in output_text
    assert not shifted.exists()
    assert asample.is_file()
    asample_lines = asample.read_text(encoding="utf-8").splitlines()
    assert asample_lines[5].split() == ["1"]
    assert asample_lines[6:8] == ["Selective dynamics", "Direct"]
    assert asample_lines[8].split()[3:6] == ["T", "T", "T"]

    quiet_exit = cli_main(["pre-check", str(source), "--quiet"])
    quiet_output = capsys.readouterr().out
    assert quiet_exit == 0
    assert "[1/4]" not in quiet_output
    assert "Status: READY" in quiet_output
    assert f"Updated ASAMPLE: {asample}" in quiet_output


@pytest.mark.skip(reason="CLI integration is delivered in the follow-up application PR")
def test_cli_unsafe_explains_failure_without_generating_file(tmp_path, capsys):
    source = tmp_path / "POSCAR"
    shifted = tmp_path / "POSCAR-autobskan_shifted.vasp"
    asample = tmp_path / "ASAMPLE"
    atoms = Atoms(
        "CuCuOO",
        scaled_positions=[
            [0.10, 0.10, 0.30],
            [0.90, 0.60, 0.30],
            [0.23, 0.31, 0.41],
            [0.77, 0.81, 0.41],
        ],
        cell=np.diag([4.0, 5.0, 20.0]),
        pbc=True,
    )
    _write_poscar(source, atoms)

    exit_code = cli_main(["pre-check", str(source)])
    output_text = capsys.readouterr().out

    assert exit_code == 3
    assert "Status: UNSAFE" in output_text
    assert "Simple shift cannot resolve" in output_text
    assert "glide or screw" in output_text
    assert "Nonzero operations" in output_text
    assert f"No shifted file was generated: {shifted}" in output_text
    assert f"No ASAMPLE was generated or modified: {asample}" in output_text
    assert not shifted.exists()
    assert not asample.exists()


@pytest.mark.skip(reason="CLI integration is delivered in the follow-up application PR")
def test_cli_rejects_tilted_c_before_generating_files(tmp_path, capsys):
    source = tmp_path / "POSCAR"
    shifted = tmp_path / "POSCAR-autobskan_shifted.vasp"
    asample = tmp_path / "ASAMPLE"
    atoms = Atoms(
        "Cu",
        scaled_positions=[[0.0, 0.0, 0.0]],
        cell=Cell.fromcellpar([4, 5, 20, 90, 104, 90]),
        pbc=True,
    )
    _write_poscar(source, atoms)

    exit_code = cli_main(["pre-check", str(source)])
    output_text = capsys.readouterr().out

    assert exit_code == 2
    assert "Status: UNSUPPORTED_GEOMETRY" in output_text
    assert "basis operations: not evaluated" in output_text
    assert "Pre-check error: c must be perpendicular to the ab plane" in output_text
    assert "cannot guarantee" in output_text
    assert f"No shifted file was generated: {shifted}" in output_text
    assert f"No ASAMPLE was generated or modified: {asample}" in output_text
    assert not shifted.exists()
    assert not asample.exists()


def test_existing_scf_guard_requires_restart_from_precheck(tmp_path):
    scf = tmp_path / "finished_scf"
    scf.mkdir()
    _write_poscar(scf / "POSCAR", _shifted_mirror_structure())

    with pytest.raises(ExistingScfRequiresRestartError) as error:
        validate_existing_scf_for_bskan(scf)

    message = str(error.value)
    assert "autobskan pre-check" in message
    assert "new SCF" in message
    assert "Do not translate POSCAR/ASAMPLE alone" in message
    assert "WAVECAR, CHGCAR, or WAVSAMPLE" in message


def test_existing_scf_guard_rejects_tilted_c_geometry(tmp_path):
    scf = tmp_path / "finished_scf"
    scf.mkdir()
    atoms = Atoms(
        "Cu",
        scaled_positions=[[0.0, 0.0, 0.0]],
        cell=Cell.fromcellpar([4, 5, 20, 90, 104, 90]),
        pbc=True,
    )
    _write_poscar(scf / "POSCAR", atoms)

    with pytest.raises(ExistingScfRequiresRestartError) as error:
        validate_existing_scf_for_bskan(scf)

    message = str(error.value)
    assert "outside the validated bSKAN pre-check domain" in message
    assert "c lattice vector must be perpendicular" in message
    assert "autobskan pre-check" in message
    assert "start a new SCF calculation" in message


def test_calculation_entry_stops_before_creating_work_directories(tmp_path, monkeypatch):
    scf = tmp_path / "finished_scf"
    scf.mkdir()
    _write_poscar(scf / "POSCAR", _shifted_mirror_structure())
    input_file = tmp_path / "bskan.in"
    input_file.write_text(
        "TASK = CALCULATION\n"
        "VASP_COMMAND = true\n"
        "BSKAN_COMMAND = true\n"
        f"SCF_PATH = {scf}\n"
        f"TIP_PATH = {tmp_path / 'tip'}\n"
        "METHOD = TH\n"
        "BIAS = 0.0\n",
        encoding="ascii",
    )
    monkeypatch.chdir(tmp_path)

    with pytest.raises(ExistingScfRequiresRestartError, match="new SCF"):
        connect.main(input_file)

    assert not (tmp_path / "1_nonscf").exists()
    assert not (tmp_path / "2_bskan").exists()


@pytest.mark.skip(reason="CLI integration is delivered in the follow-up application PR")
def test_calculation_cli_prints_clean_restart_error(tmp_path, monkeypatch, capsys):
    scf = tmp_path / "finished_scf"
    scf.mkdir()
    _write_poscar(scf / "POSCAR", _shifted_mirror_structure())
    input_file = tmp_path / "bskan.in"
    input_file.write_text(
        "TASK = CALCULATION\n"
        "VASP_COMMAND = true\n"
        "BSKAN_COMMAND = true\n"
        f"SCF_PATH = {scf}\n"
        f"TIP_PATH = {tmp_path / 'tip'}\n"
        "METHOD = TH\n"
        "BIAS = 0.0\n",
        encoding="ascii",
    )
    monkeypatch.chdir(tmp_path)

    exit_code = cli_main(["--input", str(input_file)])
    captured = capsys.readouterr()

    assert exit_code == 4
    assert captured.out == ""
    assert captured.err.startswith("[ERROR] Existing SCF data")
    assert "autobskan pre-check" in captured.err
    assert "new SCF" in captured.err


def test_existing_scf_guard_detects_poscar_chgcar_mismatch(tmp_path):
    scf = tmp_path / "finished_scf"
    scf.mkdir()
    poscar_atoms = Atoms(
        "Cu", scaled_positions=[[0.0, 0.0, 0.0]], cell=[4, 5, 20], pbc=True
    )
    chgcar_atoms = poscar_atoms.copy()
    chgcar_atoms.translate([0.2, 0.0, 0.0])
    _write_poscar(scf / "POSCAR", poscar_atoms)
    # read_vasp only needs the VASP structure header for this consistency check.
    _write_poscar(scf / "CHGCAR", chgcar_atoms)

    with pytest.raises(ExistingScfRequiresRestartError, match="do not match"):
        validate_existing_scf_for_bskan(scf)


def test_writer_refuses_in_place_modification(tmp_path):
    poscar = tmp_path / "POSCAR"
    _write_poscar(poscar, _shifted_mirror_structure())

    with pytest.raises(ValueError, match="in place"):
        write_prechecked_poscar(poscar, poscar)

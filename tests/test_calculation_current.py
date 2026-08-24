from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from autobskan.calculation import (
    asample_to_current,
    build_current_head_from_asample,
)

ROOT = Path(__file__).resolve().parents[1]


def _write_minimal_calculation(directory: Path, cursave: str) -> Path:
    asample = directory / "ASAMPLE"
    asample.write_text(
        "synthetic\n"
        "1.0\n"
        "2.0 0.0 0.0\n"
        "0.0 2.0 0.0\n"
        "0.0 0.0 10.0\n"
        "1\n"
        "Selective dynamics\n"
        "Direct\n"
        "0.0 0.0 0.2 T T T\n",
        encoding="ascii",
    )
    (directory / "INSCAN").write_text(
        "CHEN\n"
        "PIVOT = 0.0 0.0\n"
        "CELL = 1.0 1.0\n"
        "ZVACUUM = 3.0\n",
        encoding="ascii",
    )
    (directory / "CURSAVE").write_text(cursave, encoding="ascii")
    return asample


def test_real_bskan_header_is_byte_identical():
    asample = ROOT / "examples/2_hexagonal/1_p2.vasp"
    current = ROOT / "examples/2_hexagonal/1_p2_nonscf_for_STM_bskan_th_0.5_V_current"

    generated = build_current_head_from_asample(asample, (34, 60, 101))
    reference = "".join(current.read_text(encoding="ascii").splitlines(keepends=True)[:30])

    assert generated.encode("ascii") == reference.encode("ascii")
    assert len(generated.encode("ascii")) == 983
    assert hashlib.sha256(generated.encode("ascii")).hexdigest() == (
        "f9b9dc2469cd1204d12f8a41f5e4a1b77a9596a999ec32be1ec11b62496d0e40"
    )


def test_asample_only_argument_discovers_metadata_and_reorders_cursave(tmp_path):
    asample = _write_minimal_calculation(
        tmp_path,
        "Backup file\n"
        "0.0000 0.0000\n1.0\n2.0\n"
        "0.0000 2.0000\n3.0\n4.0\n"
        "2.0000 0.0000\n5.0\n6.0\n"
        "2.0000 2.0000\n7.0\n8.0\n",
    )

    output = asample_to_current(asample)
    lines = output.read_text(encoding="ascii").splitlines(keepends=True)
    grid_line = next(index for index, line in enumerate(lines) if line.split() == ["2", "2", "2"])

    assert output == tmp_path / "CURRENT"
    assert "".join(lines[grid_line + 1 :]) == (
        "    0.1000E+01    0.5000E+01    0.3000E+01    0.7000E+01    0.2000E+01\n"
        "    0.6000E+01    0.4000E+01    0.8000E+01\n"
    )


def test_incomplete_cursave_fails_with_a_clear_error(tmp_path):
    asample = _write_minimal_calculation(
        tmp_path,
        "Backup file\n"
        "0.0000 0.0000\n1.0\n2.0\n"
        "0.0000 2.0000\n3.0\n",
    )

    with pytest.raises(ValueError, match="inconsistent z-grid lengths"):
        asample_to_current(asample)

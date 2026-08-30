from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from autobskan.calculation import (
    asample_to_current,
    build_current_head_from_asample,
    connect,
)


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


def test_synthetic_bskan_header_is_byte_identical(tmp_path):
    asample = _write_minimal_calculation(tmp_path, "Backup file\n")
    generated = build_current_head_from_asample(asample, (2, 2, 3))
    reference = (
        " synthetic                               \n"
        "   1.0   \n"
        "    1.999999    0.000000    0.000000\n"
        "    0.000000    1.999999    0.000000\n"
        "    0.000000    0.000000    0.158753\n"
        "    1\n"
        " Direct\n"
        "    0.000000    0.000000   -6.299089\n"
        "  \n"
        "    2    2    3\n"
    )

    assert generated.encode("ascii") == reference.encode("ascii")
    assert len(generated.encode("ascii")) == 233
    assert hashlib.sha256(generated.encode("ascii")).hexdigest() == (
        "ba17e0cad5c68644c12b17eaeeda810d62a8049a60645f03c043b43a254f4e59"
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


def test_connect_ctoc_uses_asample_converter_directly(tmp_path, monkeypatch):
    asample = tmp_path / "ASAMPLE"
    expected = tmp_path / "CURRENT"
    calls = []

    def fake_converter(path):
        calls.append(Path(path))
        return expected

    monkeypatch.setattr(connect, "asample_to_current", fake_converter)

    assert connect.ctoc(asample) == expected
    assert calls == [asample]

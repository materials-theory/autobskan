from __future__ import annotations

from types import SimpleNamespace

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from autobskan.gui.gui import _bskanin_export_text
from autobskan.image.serial_run import (
    _constant_current_iso_values,
    single_current,
)
from autobskan.io.input import Options


def _synthetic_current():
    current_grid = np.asarray(
        [
            [[6.0, 7.0], [8.0, 9.0]],
            [[4.0, 4.5], [5.0, 5.5]],
            [[2.0, 2.0], [2.0, 2.0]],
        ],
        dtype=float,
    )
    return SimpleNamespace(
        cur_3d=current_grid,
        iso_min=2.0,
        iso_max=6.0,
        topmost_atom=-0.25,
        cell=np.diag([2.0, 2.0, 3 * 0.052918]),
        cellpar=np.asarray([2.0, 2.0, 3 * 0.052918, 90.0, 90.0, 90.0]),
        grids=np.asarray([2, 2, 3]),
        filename="synthetic_CURRENT",
    )


def test_iso_auto_true_uses_logscale_decades():
    values = _constant_current_iso_values(
        {"ISO_AUTO": True, "ISO": 5},
        iso_min=1.0e2,
        iso_max=1.0e6,
    )

    assert values == pytest.approx([1.0e3, 1.0e4, 1.0e5])


def test_iso_auto_logscale_returns_decades_as_actual_current_values():
    values = _constant_current_iso_values(
        {"ISO_AUTO": "LOGSCALE", "ISO": 5},
        iso_min=1.0e2,
        iso_max=1.0e6,
    )

    assert values == pytest.approx([1.0e3, 1.0e4, 1.0e5])


def test_iso_auto_linear_returns_evenly_spaced_current_values_inside_range():
    values = _constant_current_iso_values(
        {"ISO_AUTO": "LINEAR", "ISO": 5},
        iso_min=2.0,
        iso_max=6.0,
    )

    assert len(values) == 5
    assert np.all(np.asarray(values) > 2.0)
    assert np.all(np.asarray(values) < 6.0)
    assert np.diff(values) == pytest.approx(np.full(4, 2.0 / 3.0))


def test_iso_auto_logscale_reports_range_without_an_internal_decade():
    with pytest.raises(ValueError, match="no power of ten strictly inside"):
        _constant_current_iso_values(
            {"ISO_AUTO": "LOGSCALE", "ISO": 5},
            iso_min=2.0,
            iso_max=6.0,
        )


@pytest.mark.parametrize("tag", ["TRUE", "LOGSCALE", "L"])
def test_options_maps_true_and_logscale_tags_to_logscale(tmp_path, tag):
    input_path = tmp_path / f"{tag.lower()}-bskan.in"
    input_path.write_text(f"ISO_AUTO = {tag}\n", encoding="ascii")

    options = Options(input_path).option_dict

    assert options["ISO_AUTO"] == "LOGSCALE"


@pytest.mark.parametrize("tag", ["LINEAR", "EQUAL"])
def test_options_maps_equal_spacing_tags_to_linear(tmp_path, tag):
    input_path = tmp_path / f"{tag.lower()}-bskan.in"
    input_path.write_text(
        f"ISO_AUTO = {tag}\nISO = 5\n",
        encoding="ascii",
    )

    options = Options(input_path).option_dict

    assert options["ISO_AUTO"] == "LINEAR"
    assert options["ISO"] == 5


def test_single_current_iso_auto_generates_every_image_and_warns_once(tmp_path, capsys):
    current = _synthetic_current()
    options = Options(None)
    options.option_dict.update(
        {
            "IMAGE_MODE": "CONSTANT CURRENT",
            "ISO_AUTO": "LINEAR",
            "ISO": 5,
            "POSCAR": "",
            "ITERATION": [1, 1],
            "DISPLAY_ATOMS": False,
            "DISPLAY_CELL": False,
            "DISPLAY_CBAR": False,
            "BLUR_SIGMA": 0.0,
            "CONTOUR_RESOLUTION": 20,
        }
    )

    outputs = single_current(current, options, image_dir=tmp_path / "images")

    assert len(outputs) == 5
    assert all(path.endswith(".png") for path in outputs)
    assert all((tmp_path / "images" / path.split("/")[-1]).stat().st_size > 0 for path in outputs)
    captured = capsys.readouterr().out
    assert captured.count("CURRENT header check") == 1


@pytest.mark.parametrize(
    ("input_source", "expected_key", "unexpected_key"),
    [("BSKAN", "CURRENT", "VOLUME"), ("VASP", "VOLUME", "CURRENT")],
)
def test_gui_bskanin_export_uses_source_specific_absolute_paths(
    tmp_path,
    input_source,
    expected_key,
    unexpected_key,
):
    volume = tmp_path / ("sample_current" if input_source == "BSKAN" else "PARCHG")
    structure = tmp_path / "POSCAR"
    volume.write_text("data", encoding="ascii")
    structure.write_text("structure", encoding="ascii")

    content = _bskanin_export_text(
        simulation="STM",
        input_source=input_source,
        volume_path=str(volume),
        mode="CONSTANT CURRENT",
        iso_value=1.0e-5,
        fit_radius=0.5,
        cmap_name="afmhot",
        contrast=0.0,
        brightness=0.0,
        poscar_path=str(structure),
        iteration=(3, 3),
        gaussian_sigma=0.0,
        gamma=84.315324,
        display_atoms=False,
        layers=1,
        radius_type="ATOMIC",
        size_ratio=10,
        display_cell=True,
    )

    assert f"{expected_key} = {volume.resolve()}" in content
    assert f"{unexpected_key} =" not in content
    assert f"POSCAR = {structure.resolve()}" in content
    assert "ATOM_ADDINFO" not in content
    assert "CONTOUR_RESOLUTION" not in content

    exported_input = tmp_path / f"{input_source.lower()}-bskan.in"
    exported_input.write_text(content, encoding="ascii")
    parsed = Options(exported_input).option_dict
    parsed_path_key = "CURRENT" if input_source == "BSKAN" else "VOLUME"
    assert parsed[parsed_path_key] == str(volume.resolve())
    assert parsed["POSCAR"] == str(structure.resolve())
    assert parsed["ITERATION"] == [3, 3]

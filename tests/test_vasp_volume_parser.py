from __future__ import annotations

import numpy as np

from autobskan.image import AR


def _volume_text(grids, first, second=None, row_width=5):
    lines = [
        "test system",
        "1.0",
        "1.0 0.0 0.0",
        "0.0 1.0 0.0",
        "0.0 0.0 2.0",
        "H",
        "1",
        "Direct",
        "0.0 0.0 0.0",
        "",
        " ".join(str(value) for value in grids),
    ]

    def append_dataset(values):
        for start in range(0, len(values), row_width):
            lines.append(
                " ".join(f"{value:.10E}" for value in values[start : start + row_width])
            )

    append_dataset(first)
    if second is not None:
        lines.extend(
            [
                "augmentation occupancies 1 1",
                " ".join(str(value) for value in grids),
            ]
        )
        append_dataset(second)
    return "\n".join(lines) + "\n"


def test_chgcar_streams_non_divisible_vasp_grid(tmp_path):
    path = tmp_path / "LOCPOT"
    values = np.arange(8, dtype=float) + 0.25
    path.write_text(_volume_text((2, 2, 2), values), encoding="ascii")

    volumetric = AR.Chgcar(str(path))

    assert volumetric.grids == [2, 2, 2]
    assert volumetric.pot.shape == (2, 2, 2)
    assert not volumetric.magnetic
    np.testing.assert_allclose(volumetric.pot.reshape(-1), values)


def test_chgcar_streams_second_locpot_dataset(tmp_path):
    path = tmp_path / "LOCPOT"
    first = np.arange(8, dtype=float)
    second = first + 10.0
    path.write_text(
        _volume_text((2, 2, 2), first, second=second),
        encoding="ascii",
    )

    volumetric = AR.Chgcar(str(path))

    assert volumetric.magnetic
    np.testing.assert_allclose(volumetric.pot[0].reshape(-1), first)
    np.testing.assert_allclose(volumetric.pot[1].reshape(-1), second)

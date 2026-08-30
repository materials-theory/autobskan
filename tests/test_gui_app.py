from __future__ import annotations

import base64
import threading
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
import pytest
from ase import Atoms
from dash.development.base_component import Component

import autobskan
from autobskan import main
from autobskan.gui import gui as gui_module
from autobskan.gui.cache import SurfaceCache
from autobskan.gui.gui import (
    _add_line_selection_traces,
    _BrowserLifecycle,
    _cached_surface_figure,
    _clip_line_points,
    _display_surface,
    _effective_color_limits,
    _effective_input_source,
    _effective_path,
    _export_colorbar_figure,
    _export_file_details,
    _gaussian_blur_surface,
    _is_managed_upload,
    _line_profile,
    _line_profile_csv,
    _load_current,
    _named_volume_source,
    _overlay_traces,
    _save_upload,
    _scalar_field_figure,
    _structure_override_is_confirmed,
    _surface_apparent_barrier_from_parchg,
    _surface_from_figure,
    _surface_stm,
    _validate_export_dimensions,
    _validate_simulation_input,
    _validate_volume_source,
    build_app,
)
from autobskan.image import lwfplot
from autobskan.io.input import Options


def _walk_components(node):
    if not isinstance(node, Component):
        return
    yield node
    children = getattr(node, "children", None)
    if isinstance(children, (list, tuple)):
        for child in children:
            yield from _walk_components(child)
    elif children is not None:
        yield from _walk_components(children)


def _component_by_id(layout, component_id):
    for component in _walk_components(layout):
        if getattr(component, "id", None) == component_id:
            return component
    raise AssertionError(f"Missing component id: {component_id}")


def _write_synthetic_current(path):
    nx, ny, nz = 4, 3, 5
    lateral_scale = 1.0 + 0.01 * np.arange(nx * ny).reshape(ny, nx)
    decay = np.asarray([10.0, 7.0, 4.5, 2.5, 1.0])[:, None, None]
    values = (decay * lateral_scale).reshape(-1)
    rows = [
        " ".join(f"{value:.8e}" for value in values[index : index + 5])
        for index in range(0, values.size, 5)
    ]
    header = (
        "Synthetic CURRENT\n"
        "1.0\n"
        "2.0 0.0 0.0\n"
        "0.5 1.5 0.0\n"
        "0.0 0.0 5.2918\n"
        "1\n"
        "Selective dynamics\n"
        "Direct\n"
        "0.25 0.25 -0.094485 T T T\n"
        "\n"
        f"{nx} {ny} {nz}\n"
    )
    path.write_text(header + "\n".join(rows) + "\n", encoding="ascii")


def test_layout_ids_are_unique_and_debug_console_is_hidden():
    app = build_app(debug_mode=False)
    app._setup_server()
    ids = [
        str(component.id)
        for component in _walk_components(app.layout)
        if getattr(component, "id", None) is not None
    ]
    duplicates = [key for key, count in Counter(ids).items() if count > 1]

    assert not duplicates
    assert _component_by_id(app.layout, "render-button") is not None
    assert _component_by_id(app.layout, "main-image") is not None
    assert _component_by_id(app.layout, "input-source") is not None
    assert _component_by_id(app.layout, "display-quality").value == "BALANCED"
    assert _component_by_id(app.layout, "export-format").value == "png"
    assert _component_by_id(app.layout, "export-width").value == 4096
    assert _component_by_id(app.layout, "export-height").value == 2048
    assert _component_by_id(app.layout, "colorbar-width").value == 1600
    assert _component_by_id(app.layout, "colorbar-height").value == 80
    assert _component_by_id(app.layout, "render-overlay") is not None
    assert _component_by_id(app.layout, "user-error-banner").role == "alert"
    assert _component_by_id(app.layout, "upload-volume-button") is not None
    assert _component_by_id(app.layout, "upload-volume-progress") is not None
    assert _component_by_id(app.layout, "controls-toggle-button") is not None
    assert _component_by_id(app.layout, "new-window-button") is not None
    assert _component_by_id(app.layout, "controls-backdrop") is not None
    assert _component_by_id(app.layout, "structure-override-modal") is not None
    assert _component_by_id(app.layout, "embedded-structure-banner") is not None
    assert _component_by_id(app.layout, "save-line-profile-button") is not None
    assert _component_by_id(app.layout, "download-line-profile") is not None
    assert _component_by_id(app.layout, "fermi-level-input").value is None
    assert _component_by_id(app.layout, "fit-radius-slider").value == 0.5
    assert _component_by_id(app.layout, "fit-radius-input").value == 0.5
    assert _component_by_id(app.layout, "isosurface-range-default-store").data == -2.5
    assert _component_by_id(app.layout, "color-range-mode").value == "AUTO"
    assert _component_by_id(app.layout, "color-vmin-input").disabled is True
    assert _component_by_id(app.layout, "color-vmax-input").disabled is True
    assert _component_by_id(app.layout, "control-tabs").value == "data"
    simulation_control = _component_by_id(app.layout, "simulation-type")
    assert [
        entry["value"]
        for entry in simulation_control.data
    ] == ["STM", "PHI_APP", "LWF"]
    assert simulation_control.data[2]["label"].children == r"$\Phi_{\rm loc}$"

    debug_console = next(
        component
        for component in _walk_components(app.layout)
        if getattr(component, "className", None) == "debug-console"
    )
    assert debug_console.style["display"] == "none"

    error_callback = next(
        callback
        for callback in app._callback_list
        if callback["clientside_function"]
        == {"namespace": "autobskan", "function_name": "surfaceErrors"}
    )
    assert "user-error-banner.children" in error_callback["output"]


def test_surface_cache_uses_lru_and_idle_expiration():
    now = [0.0]
    cache = SurfaceCache(
        max_entries=2,
        max_bytes=1024 * 1024,
        ttl_seconds=10.0,
        clock=lambda: now[0],
    )
    first = {"z": np.zeros((2, 2))}
    second = {"z": np.ones((2, 2))}
    third = {"z": np.full((2, 2), 2.0)}

    first_key = cache.put(first)
    second_key = cache.put(second)
    assert cache.get(first_key) is first

    third_key = cache.put(third)
    assert cache.get(second_key) is None
    assert cache.get(first_key) is first
    assert cache.get(third_key) is third

    now[0] = 11.0
    assert cache.get(first_key) is None
    assert cache.get(third_key) is None
    assert cache.entry_count == 0


def test_large_vasp_volume_load_is_coalesced(tmp_path, monkeypatch):
    volume_path = tmp_path / "LOCPOT"
    volume_path.write_text("placeholder", encoding="ascii")
    sentinel = object()
    calls = []

    def fake_chgcar(path):
        calls.append(path)
        time.sleep(0.05)
        return sentinel

    gui_module._load_locpot_cached.cache_clear()
    monkeypatch.setattr(gui_module.AR, "Chgcar", fake_chgcar)
    try:
        with ThreadPoolExecutor(max_workers=4) as pool:
            loaded = list(
                pool.map(
                    gui_module._load_locpot,
                    [str(volume_path)] * 4,
                )
            )
    finally:
        gui_module._load_locpot_cached.cache_clear()

    assert loaded == [sentinel] * 4
    assert calls == [str(volume_path)]


def test_chunked_volume_upload_writes_without_base64_payload():
    app = build_app(debug_mode=False)
    client = app.server.test_client()
    payload = b"0123456789abcdef"

    started = client.post(
        "/_autobskan-upload/start",
        json={"filename": "LOCPOT", "size": len(payload)},
    )
    assert started.status_code == 200
    upload_id = started.get_json()["upload_id"]

    first = client.post(
        f"/_autobskan-upload/chunk/{upload_id}",
        data=payload[:7],
        headers={"X-Upload-Offset": "0"},
    )
    second = client.post(
        f"/_autobskan-upload/chunk/{upload_id}",
        data=payload[7:],
        headers={"X-Upload-Offset": "7"},
    )
    assert first.get_json()["received"] == 7
    assert second.get_json()["received"] == len(payload)

    finished = client.post(f"/_autobskan-upload/finish/{upload_id}")
    assert finished.status_code == 200
    uploaded_path = Path(finished.get_json()["path"])
    try:
        assert _is_managed_upload(str(uploaded_path))
        assert uploaded_path.read_bytes() == payload
    finally:
        uploaded_path.unlink(missing_ok=True)


def test_render_callback_has_live_visualization_triggers():
    app = build_app(debug_mode=False)
    render_callback = next(
        callback
        for output, callback in app.callback_map.items()
        if "main-image.figure" in output
        and "colorbar-graph.figure" in output
        and "plot-meta-store.data" in output
    )

    assert [entry["id"] for entry in render_callback["inputs"]] == [
        "render-button",
        "volume-path",
        "volume-upload-ready-store",
        "structure-path",
        "structure-upload-ready-store",
        "structure-override-decision-store",
        "simulation-type",
        "input-source",
        "stm-mode",
        "fermi-level-input",
        "isosurface-render-store",
        "view-options",
        "scalebar-length",
        "xy-upsample",
        "display-quality",
        "render-request-store",
    ]
    render_input_ids = {entry["id"] for entry in render_callback["inputs"]}
    render_state_ids = {entry["id"] for entry in render_callback["state"]}
    assert "line-points-store" not in render_input_ids
    assert {"colormap-dropdown", "brightness-slider", "contrast-slider"}.isdisjoint(
        render_input_ids
    )
    assert {
        "colormap-dropdown",
        "brightness-slider",
        "contrast-slider",
        "isosurface-slider",
        "atom-radius-scale",
        "gaussian-sigma",
        "fit-radius-slider",
        "gamma-input",
        "layers-input",
        "repeat-x",
        "repeat-y",
        "atom-radius-type",
        "grid-line-width",
        "grid-line-color",
        "color-range-mode",
        "color-vmin-input",
        "color-vmax-input",
    }.issubset(render_state_ids)

    slider_inputs = {
        "isosurface-slider": "isosurface-input",
        "fit-radius-slider": "fit-radius-input",
        "brightness-slider": "brightness-input",
        "contrast-slider": "contrast-input",
        "gaussian-sigma": "gaussian-sigma-input",
        "atom-radius-scale": "atom-radius-scale-input",
    }
    for slider_id, input_id in slider_inputs.items():
        assert _component_by_id(app.layout, slider_id).updatemode == "drag"
        assert _component_by_id(app.layout, input_id).value == _component_by_id(
            app.layout, slider_id
        ).value


def test_appearance_updates_are_clientside_and_exports_use_full_surface_cache():
    app = build_app(debug_mode=False)
    appearance_callback = next(
        callback
        for callback in app._callback_list
        if callback["clientside_function"]
        == {"namespace": "autobskan", "function_name": "updateAppearance"}
    )
    assert appearance_callback["clientside_function"] == {
        "namespace": "autobskan",
        "function_name": "updateAppearance",
    }
    assert [entry["id"] for entry in appearance_callback["inputs"]] == [
        "appearance-store",
        "brightness-slider",
        "color-range-mode",
        "color-vmin-input",
        "color-vmax-input",
    ]
    assert "line-profile-graph.figure" in appearance_callback["output"]
    with pytest.raises(AssertionError):
        _component_by_id(app.layout, "blur-render-request-store")

    export_callback = next(
        callback
        for output, callback in app.callback_map.items()
        if "download-image.data" in output
    )
    export_state_ids = [entry["id"] for entry in export_callback["state"]]
    assert "surface-key-store" in export_state_ids
    assert {"export-format", "export-width", "export-height"}.issubset(
        export_state_ids
    )
    assert "export-scale" not in export_state_ids
    assert "main-image" not in export_state_ids

    loading_callback = next(
        callback
        for callback in app._callback_list
        if callback["clientside_function"]
        == {"namespace": "autobskan", "function_name": "beginRender"}
    )
    assert "volume-path" in {entry["id"] for entry in loading_callback["inputs"]}
    assert "isosurface-render-store" in {
        entry["id"] for entry in loading_callback["inputs"]
    }
    assert "render-request-store" in {
        entry["id"] for entry in loading_callback["inputs"]
    }
    assert {
        "layers-input",
        "repeat-x",
        "repeat-y",
        "atom-radius-type",
        "grid-line-width",
        "grid-line-color",
    }.isdisjoint({entry["id"] for entry in loading_callback["inputs"]})
    render_control_callback = next(
        callback
        for callback in app._callback_list
        if callback["clientside_function"]
        == {"namespace": "autobskan", "function_name": "scheduleRenderControls"}
    )
    assert [entry["id"] for entry in render_control_callback["inputs"]] == [
        "isosurface-slider",
        "gaussian-sigma",
        "atom-radius-scale",
        "fit-radius-slider",
    ]
    assert [entry["id"] for entry in render_control_callback["state"]] == [
        "view-options",
        "simulation-type",
    ]
    context_render_callback = next(
        callback
        for callback in app._callback_list
        if callback["clientside_function"]
        == {"namespace": "autobskan", "function_name": "scheduleContextRender"}
    )
    assert [entry["id"] for entry in context_render_callback["inputs"]] == [
        "gamma-input",
        "layers-input",
        "repeat-x",
        "repeat-y",
        "atom-radius-type",
        "grid-line-width",
        "grid-line-color",
    ]
    range_apply_callback = next(
        callback
        for callback in app._callback_list
        if callback["clientside_function"]
        == {"namespace": "autobskan", "function_name": "applySliderRange"}
    )
    assert [entry["id"] for entry in range_apply_callback["inputs"]] == [
        "isosurface-slider",
        "isosurface-slider",
        "isosurface-range-default-store",
    ]
    range_server_callback = next(
        callback
        for callback in app._callback_list
        if "isosurface-help.children" in callback["output"]
    )
    assert "isosurface-range-default-store.data" in range_server_callback["output"]
    assert "isosurface-slider.value" not in range_server_callback["output"]

    volume_parse_callback = next(
        callback
        for callback in app._callback_list
        if "volume-status.children" in callback["output"]
    )
    assert volume_parse_callback["running"] == {
        "running": {"render-overlay.className": "render-overlay is-active"},
        "runningOff": {"render-overlay.className": "render-overlay"},
    }


def test_display_surface_caps_browser_pixels_without_changing_physical_source():
    x_axis = np.linspace(-3.0, 7.0, 800)
    y_axis = np.linspace(1.0, 5.0, 400)
    X, Y = np.meshgrid(x_axis, y_axis)
    source = np.sin(0.8 * X) + 0.25 * np.cos(1.3 * Y)
    original = source.copy()

    display_x, display_y, display_z = _display_surface(
        x_axis,
        y_axis,
        source,
        "FAST",
    )

    assert display_z.size <= 200_000
    assert display_z.shape == (display_y.size, display_x.size)
    assert display_x[[0, -1]].tolist() == pytest.approx(x_axis[[0, -1]])
    assert display_y[[0, -1]].tolist() == pytest.approx(y_axis[[0, -1]])
    np.testing.assert_array_equal(source, original)

    full_x, full_y, full_z = _display_surface(x_axis, y_axis, source, "FULL")
    np.testing.assert_array_equal(full_x, x_axis)
    np.testing.assert_array_equal(full_y, y_axis)
    np.testing.assert_array_equal(full_z, source)


def test_line_profile_uses_rectilinear_interpolation_without_griddata(monkeypatch):
    x_axis = np.linspace(0.0, 5.0, 200)
    y_axis = np.linspace(-1.0, 3.0, 160)
    X, Y = np.meshgrid(x_axis, y_axis)
    surface = X + 2.0 * Y

    def fail_griddata(*_args, **_kwargs):
        raise AssertionError("rectilinear line profile must not call griddata")

    monkeypatch.setattr(gui_module, "griddata", fail_griddata)
    figure = _line_profile(
        x_axis,
        y_axis,
        surface,
        [0.5, -0.5],
        [4.5, 2.5],
        "Line profile",
        "A",
        y_range=(-2.0, 12.0),
    )

    assert len(figure.data) == 1
    assert len(figure.data[0].x) == 1200
    assert float(figure.data[0].y[0]) == pytest.approx(-0.5)
    assert float(figure.data[0].y[-1]) == pytest.approx(9.5)
    assert list(figure.layout.yaxis.range) == pytest.approx([-2.0, 12.0])
    assert figure.layout.yaxis.autorange is False


def test_gaussian_blur_preserves_periodic_coordinates_and_peak_location():
    z = np.zeros((32, 48), dtype=float)
    z[3, 5] = 1.0
    blurred = _gaussian_blur_surface(z, 2.0)

    assert np.unravel_index(np.nanargmax(blurred), blurred.shape) == (3, 5)
    np.testing.assert_allclose(
        _gaussian_blur_surface(np.roll(z, (7, 11), axis=(0, 1)), 2.0),
        np.roll(blurred, (7, 11), axis=(0, 1)),
        atol=1e-14,
    )

    x_axis = np.linspace(-2.0, 8.0, z.shape[1])
    y_axis = np.linspace(1.5, 7.0, z.shape[0])
    X, Y = np.meshgrid(x_axis, y_axis)
    colorscale = [[0.0, "#000000"], [1.0, "#ffffff"]]
    plain, *plain_axes = _scalar_field_figure(
        X,
        Y,
        z,
        colorscale,
        0.0,
        1.0,
        False,
        0.0,
    )
    filtered, *filtered_axes = _scalar_field_figure(
        X,
        Y,
        z,
        colorscale,
        0.0,
        1.0,
        True,
        2.0,
    )

    np.testing.assert_array_equal(plain.data[0].x, filtered.data[0].x)
    np.testing.assert_array_equal(plain.data[0].y, filtered.data[0].y)
    np.testing.assert_allclose(plain_axes[2:6], filtered_axes[2:6])
    np.testing.assert_allclose(filtered.data[0].z, blurred)
    assert len(filtered.layout.images) == 0
    assert plain_axes[-1] is False
    assert filtered_axes[-1] is True


def test_line_profile_csv_contains_the_plotted_x_and_y_columns():
    figure = _line_profile(
        np.linspace(0.0, 2.0, 5),
        np.linspace(0.0, 1.0, 3),
        np.arange(15, dtype=float).reshape(3, 5),
        [0.0, 0.0],
        [2.0, 1.0],
        "Line profile",
        "a.u.",
    )
    payload = _line_profile_csv(figure.to_plotly_json())
    rows = payload.strip().splitlines()

    assert rows[0] == "x,y"
    assert len(rows) == len(figure.data[0].x) + 1
    first_x, first_y = [float(value) for value in rows[1].split(",")]
    last_x, last_y = [float(value) for value in rows[-1].split(",")]
    assert first_x == pytest.approx(float(figure.data[0].x[0]))
    assert first_y == pytest.approx(float(figure.data[0].y[0]))
    assert last_x == pytest.approx(float(figure.data[0].x[-1]))
    assert last_y == pytest.approx(float(figure.data[0].y[-1]))


def test_line_points_are_preserved_and_clipped_for_a_reloaded_surface():
    points = [[-1.0, 2.0], [12.0, 9.0]]
    clipped = _clip_line_points(
        points,
        np.linspace(0.0, 10.0, 20),
        np.linspace(3.0, 8.0, 10),
    )
    assert clipped == [[0.0, 3.0], [10.0, 8.0]]


def test_structure_override_confirmation_is_scoped_to_the_vasp_volume():
    decision = {
        "path": "/tmp/POSCAR",
        "volume_path": "/tmp/LOCPOT-A",
        "decision": "proceed",
    }
    assert _structure_override_is_confirmed(
        "/tmp/POSCAR",
        decision,
        "/tmp/LOCPOT-A",
    )
    assert not _structure_override_is_confirmed(
        "/tmp/POSCAR",
        decision,
        "/tmp/LOCPOT-B",
    )


def test_manual_color_limits_override_brightness_and_invalid_values_fall_back():
    values = np.asarray([-2.0, 4.0])

    assert _effective_color_limits(values, 0.5, "AUTO") == pytest.approx(
        (-5.0, 4.0)
    )
    assert _effective_color_limits(
        values,
        0.5,
        "MANUAL",
        0.25,
        1.75,
    ) == pytest.approx((0.25, 1.75))
    assert _effective_color_limits(
        values,
        0.5,
        "MANUAL",
        2.0,
        1.0,
    ) == pytest.approx((-2.0, 4.0))


def test_colorbar_export_is_an_opaque_tickless_strip_without_margins():
    figure = _export_colorbar_figure(
        [[0.0, "#000000"], [1.0, "#ffffff"]],
        -1.0,
        2.0,
    )

    assert figure.layout.paper_bgcolor == "#ffffff"
    assert figure.layout.plot_bgcolor == "#ffffff"
    assert figure.layout.margin.to_plotly_json() == {
        "b": 0,
        "l": 0,
        "pad": 0,
        "r": 0,
        "t": 0,
    }
    assert figure.layout.xaxis.visible is False
    assert figure.layout.yaxis.visible is False
    assert figure.data[0].showscale is False
    assert np.asarray(figure.data[0].z).shape == (2, 1024)
    assert float(figure.data[0].zmin) == pytest.approx(-1.0)
    assert float(figure.data[0].zmax) == pytest.approx(2.0)


def test_export_dimensions_and_output_formats_are_explicit_and_bounded():
    assert _validate_export_dimensions(4096, 2048) == (4096, 2048)
    assert _export_file_details("jpeg", "surface") == (
        "jpeg",
        "surface.jpg",
        "image/jpeg",
    )
    assert _export_file_details("svg", "surface") == (
        "svg",
        "surface.svg",
        "image/svg+xml",
    )
    with pytest.raises(ValueError, match="60,000,000"):
        _validate_export_dimensions(10000, 10000)
    with pytest.raises(ValueError, match="Output format"):
        _export_file_details("tiff", "surface")


def test_browser_lifecycle_stops_only_after_the_last_window_closes():
    stopped = threading.Event()
    lifecycle = _BrowserLifecycle(
        enabled=True,
        shutdown_callback=stopped.set,
        grace_seconds=0.03,
    )

    assert lifecycle.register("window-a-123") == 1
    assert lifecycle.register("window-b-123") == 2
    assert lifecycle.close("window-a-123") == 1
    assert not stopped.wait(0.06)
    assert lifecycle.close("window-b-123") == 0
    assert stopped.wait(0.3)


def test_browser_lifecycle_tracks_reconnected_streams_independently():
    stopped = threading.Event()
    lifecycle = _BrowserLifecycle(
        enabled=True,
        shutdown_callback=stopped.set,
        grace_seconds=0.02,
    )

    lifecycle.register("window-a-123", "stream-old")
    lifecycle.register("window-a-123", "stream-new")
    lifecycle.close("window-a-123", "stream-old")
    assert lifecycle.active_count == 1
    assert not stopped.wait(0.06)

    lifecycle.close("window-a-123", "stream-new")
    assert stopped.wait(0.3)


def test_browser_lifecycle_routes_and_close_beacon_asset():
    stopped = threading.Event()
    app = build_app(
        debug_mode=False,
        shutdown_on_last_client=True,
        shutdown_callback=stopped.set,
        shutdown_grace=0.02,
    )
    client = app.server.test_client()

    client_id = "browser-client-123"
    assert client.post(f"/_autobskan-client/open/{client_id}").status_code == 200
    stream_response = client.get(
        f"/_autobskan-client/events/{client_id}", buffered=False
    )
    assert stream_response.status_code == 200
    assert stream_response.mimetype == "text/event-stream"
    assert next(stream_response.response) == b": keepalive\n\n"
    assert client.post(f"/_autobskan-client/close/{client_id}").status_code == 200
    assert stopped.wait(0.3)
    stream_response.close()
    assert client.post("/_autobskan-client/open/not-valid").status_code == 400

    lifecycle_source = (
        Path(__file__).resolve().parents[1]
        / "autobskan"
        / "gui"
        / "assets"
        / "client_lifecycle.js"
    ).read_text(encoding="utf-8")
    assert 'navigator.sendBeacon(closeUrl, "")' in lifecycle_source
    assert 'new window.EventSource(endpoint("events"))' in lifecycle_source
    assert 'addEventListener("pagehide"' in lifecycle_source


def test_cached_export_figure_uses_full_resolution_surface():
    full_z = np.arange(60, dtype=float).reshape(6, 10)
    surface = {
        "x_axis": np.linspace(0.0, 9.0, 10),
        "y_axis": np.linspace(0.0, 5.0, 6),
        "z": full_z,
        "display_z": full_z[::2, ::2],
        "overlay": [],
        "shapes": [],
        "blur_applied": False,
        "gaussian_sigma": 0.0,
    }
    figure = _cached_surface_figure(
        surface=surface,
        colorscale=[[0.0, "#000000"], [1.0, "#ffffff"]],
        z_min=0.0,
        z_max=59.0,
        cmap_name="Greys",
        contrast=0.0,
    )

    assert np.asarray(figure.data[0].z).shape == full_z.shape


def test_cached_export_uses_the_same_filtered_surface_as_analysis():
    raw = np.zeros((16, 20), dtype=float)
    raw[3, 4] = 1.0
    filtered = _gaussian_blur_surface(raw, 1.5)
    surface = {
        "x_axis": np.linspace(0.0, 2.0, raw.shape[1]),
        "y_axis": np.linspace(0.0, 1.0, raw.shape[0]),
        "z": raw,
        "analysis_z": filtered,
        "overlay": [],
        "shapes": [],
        "blur_applied": True,
        "gaussian_sigma": 1.5,
    }

    figure = _cached_surface_figure(
        surface=surface,
        colorscale=[[0.0, "#000000"], [1.0, "#ffffff"]],
        z_min=float(np.min(filtered)),
        z_max=float(np.max(filtered)),
        cmap_name="Greys",
        contrast=0.0,
    )

    np.testing.assert_allclose(np.asarray(figure.data[0].z), filtered)


def test_line_selection_overlay_shows_p1_before_drawing_the_line():
    one_point = go.Figure()
    _add_line_selection_traces(one_point, [[1.25, 2.5]])

    assert len(one_point.data) == 1
    assert one_point.data[0].mode == "markers+text"
    assert list(one_point.data[0].text) == ["P1"]
    assert list(one_point.data[0].x) == [1.25]
    assert list(one_point.data[0].y) == [2.5]

    two_points = go.Figure()
    _add_line_selection_traces(two_points, [[1.25, 2.5], [3.0, 4.0]])

    assert len(two_points.data) == 2
    assert two_points.data[0].mode == "lines"
    assert two_points.data[1].mode == "markers+text"
    assert list(two_points.data[1].text) == ["P1", "P2"]


def test_line_selection_uses_an_immediate_clientside_callback():
    app = build_app(debug_mode=False)
    callback = next(
        callback
        for callback in app._callback_list
        if "line-points-store.data" in callback["output"]
        and "main-image.figure" in callback["output"]
    )

    assert [entry["id"] for entry in callback["inputs"]] == [
        "main-image",
        "clear-line-button",
        "line-profile-enabled",
    ]
    assert callback["clientside_function"] == {
        "namespace": "autobskan",
        "function_name": "updateLineSelection",
    }

    volume_reset_callbacks = [
        candidate
        for output, candidate in app.callback_map.items()
        if "line-points-store.data" in output
        and {
            entry["id"] for entry in candidate["inputs"]
        }.intersection({"volume-path", "volume-upload-ready-store"})
    ]
    assert volume_reset_callbacks == []


def test_new_window_and_keyboard_panel_controls():
    app = build_app(debug_mode=False)
    new_window = _component_by_id(app.layout, "new-window-button")
    assert new_window.href == "/"
    assert new_window.target == "_blank"
    assert "noopener" in new_window.rel

    panel = next(
        callback
        for callback in app._callback_list
        if callback["clientside_function"]
        == {"namespace": "autobskan", "function_name": "toggleControls"}
    )
    assert "app-interactive-wrap.className" in panel["output"]
    assert [entry["id"] for entry in panel["inputs"]][-1] == (
        "controls-keyboard-store"
    )

    performance_source = (
        Path(__file__).resolve().parents[1]
        / "autobskan"
        / "gui"
        / "assets"
        / "gui_performance.js"
    ).read_text(encoding="utf-8")
    assert 'document.getElementById("controls-close-button")' in performance_source
    assert 'document.getElementById("controls-toggle-button")' in performance_source


def test_file_selection_never_overrides_explicit_data_source():
    app = build_app(debug_mode=False)

    source_callbacks = [
        callback
        for output, callback in app.callback_map.items()
        if "input-source.value" in output
    ]

    assert len(source_callbacks) == 1
    assert [entry["id"] for entry in source_callbacks[0]["inputs"]] == [
        "set-default-button"
    ]
    assert _effective_input_source("BSKAN", "CURRENT") == "BSKAN"
    assert _effective_input_source("VASP", "CURRENT") == "VASP"


def test_known_volume_names_produce_clear_source_validation():
    assert _named_volume_source("CURRENT") == "BSKAN"
    assert _named_volume_source("volume_20260715_CURRENT") == "BSKAN"
    assert _named_volume_source("sample_PARCHG") == "VASP"
    assert _named_volume_source("LOCPOT") == "VASP"

    with pytest.raises(ValueError, match="identified as bSKAN"):
        _validate_volume_source("VASP", "CURRENT")
    with pytest.raises(ValueError, match="identified as VASP"):
        _validate_volume_source("BSKAN", "PARCHG")


def test_simulation_input_validation_separates_phi_app_and_lwf():
    assert _validate_simulation_input("PHI_APP", "BSKAN", "CURRENT") == "BSKAN"
    assert _validate_simulation_input("PHI_APP", "VASP", "PARCHG") == "VASP"
    assert _validate_simulation_input("LWF", "VASP", "LOCPOT") == "VASP"

    with pytest.raises(ValueError, match="requires PARCHG"):
        _validate_simulation_input("PHI_APP", "VASP", "LOCPOT")
    with pytest.raises(ValueError, match="requires a LOCPOT"):
        _validate_simulation_input("LWF", "VASP", "PARCHG")
    with pytest.raises(ValueError, match="requires VASP LOCPOT"):
        _validate_simulation_input("LWF", "BSKAN", "CURRENT")


def test_local_work_function_is_locpot_slice_minus_fermi_level():
    assert lwfplot.APPARENT_BARRIER_COEFFICIENT == pytest.approx(0.952495)

    class DummyVolume:
        cell = np.diag([2.0, 3.0, 4.0])
        pot = np.stack(
            [np.full((2, 3), value, dtype=float) for value in [4.0, 5.0, 6.0, 7.0]]
        )

    raw = lwfplot.surface_potential_slice(DummyVolume(), 1.5)
    local_wf = lwfplot.local_work_function_slice(
        DummyVolume(),
        1.5,
        fermi_level=1.25,
    )

    np.testing.assert_allclose(local_wf, raw - 1.25)


def test_fermi_level_is_read_from_last_outcar_entry(tmp_path):
    outcar = tmp_path / "OUTCAR"
    outcar.write_text(
        " E-fermi : 4.1000 XC(G=0): 0.0\n"
        " Fermi energy: 4.2750 eV\n",
        encoding="ascii",
    )
    locpot = tmp_path / "LOCPOT"
    locpot.write_text("placeholder", encoding="ascii")

    assert lwfplot.fermi_energy_from_outcar(outcar) == pytest.approx(4.275)
    assert lwfplot.resolve_fermi_level(str(locpot)) == pytest.approx(4.275)
    assert lwfplot.resolve_fermi_level(str(locpot), 3.9) == pytest.approx(3.9)


def test_input_options_preserve_phi_app_and_lwf_fermi_level(tmp_path):
    phi_input = tmp_path / "phi.in"
    phi_input.write_text("SIMULATION = PHI_APP\n", encoding="ascii")
    lwf_input = tmp_path / "lwf.in"
    lwf_input.write_text(
        "SIMULATION = LWF\n"
        "FERMI_LEVEL = 4.321\n"
        "ISO_AUTO = FALSE\n"
        "ISO = 2.5\n",
        encoding="ascii",
    )

    assert Options(phi_input).option_dict["SIMULATION"] == "PHI_APP"
    lwf_options = Options(lwf_input).option_dict
    assert lwf_options["SIMULATION"] == "LWF"
    assert lwf_options["FERMI_LEVEL"] == pytest.approx(4.321)
    assert lwf_options["ISO"] == pytest.approx([2.5])


def test_invalid_typed_path_does_not_fall_back_to_stale_upload(tmp_path):
    uploaded = tmp_path / "uploaded_CURRENT"
    uploaded.write_text("uploaded", encoding="ascii")

    assert _effective_path("", str(uploaded)) == str(uploaded)
    assert _effective_path(str(tmp_path / "missing"), str(uploaded)) is None


def test_upload_reuses_identical_file_in_launch_directory(tmp_path, monkeypatch):
    payload = b"CURRENT payload"
    current = tmp_path / "CURRENT"
    current.write_bytes(payload)
    contents = "data:application/octet-stream;base64," + base64.b64encode(
        payload
    ).decode("ascii")
    monkeypatch.chdir(tmp_path)

    saved = _save_upload(contents, "CURRENT", "volume")

    assert saved == str(current)
    assert not _is_managed_upload(saved)


def test_upload_copy_is_kept_outside_working_directory(tmp_path, monkeypatch):
    payload = b"uploaded elsewhere"
    contents = "data:application/octet-stream;base64," + base64.b64encode(
        payload
    ).decode("ascii")
    monkeypatch.chdir(tmp_path)

    saved = _save_upload(contents, "CURRENT", "volume")

    assert Path(saved).read_bytes() == payload
    assert _is_managed_upload(saved)
    assert Path(saved).parent != tmp_path


def test_synthetic_current_uses_bskan_surface_pipeline(tmp_path):
    current_path = tmp_path / "CURRENT"
    _write_synthetic_current(current_path)

    current = _load_current(str(current_path))
    surface = _surface_stm(current, "CONSTANT CURRENT", 4.0)

    assert tuple(current.grids) == (4, 3, 5)
    assert surface.shape == (3, 4)
    assert np.isfinite(surface).all()


def test_repeated_atom_overlay_is_cut_to_surface_bounds():
    structure = Atoms(
        "Cu2",
        scaled_positions=((0.1, 0.1, 0.4), (0.8, 0.7, 0.6)),
        cell=((2.0, 0.0, 0.0), (0.5, 1.5, 0.0), (0.0, 0.0, 8.0)),
        pbc=True,
    )
    x_min = 0.0
    y_min = 0.0
    x_max = float(structure.cell[0, 0]) * 2
    y_max = float(structure.cell[1, 1]) * 2

    traces = _overlay_traces(
        structure=structure,
        x_min=x_min,
        x_max=x_max,
        y_min=y_min,
        y_max=y_max,
        show_atoms=True,
        show_grid=False,
        layers=1,
        radius_type="atomic",
        radius_scale=10.0,
        grid_line_width=1.0,
        grid_line_color="#000000",
    )
    marker_traces = [
        trace for trace in traces if getattr(trace, "mode", None) == "markers"
    ]
    x_positions = np.concatenate(
        [np.asarray(trace.x, dtype=float) for trace in marker_traces]
    )
    y_positions = np.concatenate(
        [np.asarray(trace.y, dtype=float) for trace in marker_traces]
    )

    assert x_positions.size > 0
    assert np.all(x_positions >= x_min)
    assert np.all(x_positions < x_max)
    assert np.all(y_positions >= y_min)
    assert np.all(y_positions < y_max)


def test_dynamic_layout_is_fresh_without_disabling_asset_cache():
    app = build_app(debug_mode=False)
    client = app.server.test_client()

    layout_response = client.get("/_dash-layout")
    asset_response = client.get("/assets/autobskan.css")

    assert "no-store" in layout_response.headers["Cache-Control"]
    assert "no-store" not in asset_response.headers.get("Cache-Control", "")


def test_surface_can_be_recovered_from_plotly_lists():
    figure = {
        "data": [
            {
                "type": "heatmap",
                "x": [0.0, 1.0, 2.0],
                "y": [0.0, 2.0],
                "z": [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
            }
        ]
    }

    X, Y, Z, x_axis, y_axis = _surface_from_figure(figure)

    assert X.shape == (2, 3)
    assert Y.shape == (2, 3)
    np.testing.assert_allclose(Z, figure["data"][0]["z"])
    np.testing.assert_allclose(x_axis, [0.0, 1.0, 2.0])
    np.testing.assert_allclose(y_axis, [0.0, 2.0])


def test_surface_can_be_recovered_from_plotly_binary_arrays():
    x_axis = np.asarray([0.0, 1.0], dtype=np.float64)
    y_axis = np.asarray([0.0, 2.0], dtype=np.float64)
    z_data = np.asarray([[1.0, 2.0], [3.0, 4.0]], dtype=np.float64)

    def encoded(array, shape=None):
        value = {
            "dtype": array.dtype.str,
            "bdata": base64.b64encode(array.tobytes()).decode("ascii"),
        }
        if shape is not None:
            value["shape"] = ",".join(str(part) for part in shape)
        return value

    figure = {
        "data": [
            {
                "type": "heatmap",
                "x": encoded(x_axis),
                "y": encoded(y_axis),
                "z": encoded(z_data, z_data.shape),
            }
        ]
    }

    _X, _Y, recovered, recovered_x, recovered_y = _surface_from_figure(figure)

    np.testing.assert_allclose(recovered, z_data)
    np.testing.assert_allclose(recovered_x, x_axis)
    np.testing.assert_allclose(recovered_y, y_axis)


def test_apparent_barrier_fit_half_width_reaches_vasp_fitter():
    class Volumetric:
        cell = np.diag([2.0, 2.0, 10.0])

    volumetric = Volumetric()
    z_axis = np.arange(100, dtype=float) * 0.1
    decay = np.exp(-1.2 * z_axis - 0.08 * z_axis**3)
    volumetric.pot = np.broadcast_to(decay[:, None, None], (100, 2, 2)).copy()

    narrow = _surface_apparent_barrier_from_parchg(
        volumetric,
        height=4.0,
        fit_radius=0.2,
    )
    wide = _surface_apparent_barrier_from_parchg(
        volumetric,
        height=4.0,
        fit_radius=1.0,
    )

    assert np.all(np.isfinite(narrow))
    assert np.all(np.isfinite(wide))
    assert not np.allclose(narrow, wide)


def test_static_export_dependencies_are_declared():
    project_root = Path(__file__).resolve().parents[1]
    requirements = (project_root / "requirements.txt").read_text(encoding="utf-8")
    setup_source = (project_root / "setup.py").read_text(encoding="utf-8")

    assert "plotly>=6.1,<6.8" in requirements
    assert "kaleido>=1.0,<2" in requirements
    assert '"plotly>=6.1,<6.8"' in setup_source
    assert '"kaleido>=1.0,<2"' in setup_source


def test_cli_versions_match_package_version():
    assert autobskan.__version__ == "1.3.5"
    assert main.__version__ == autobskan.__version__

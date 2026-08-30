from __future__ import annotations

import _thread
import base64
import binascii
import csv
import datetime
import importlib.metadata as importlib_metadata
import io
import logging
import os
import re
import secrets
import tempfile
import threading
import time
from functools import lru_cache

from autobskan.gui.frontend import _prepare_matplotlib_cache

_prepare_matplotlib_cache()

import dash
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
import plotly.io as pio
import scipy.ndimage as ndimage
from dash import ClientsideFunction, Input, Output, State, dcc, no_update
from dash.exceptions import PreventUpdate
from flask import Response, jsonify, request, stream_with_context
from matplotlib import colormaps
from matplotlib import colors as mcolors
from scipy.interpolate import PchipInterpolator, RegularGridInterpolator, griddata

import autobskan
from autobskan.gui.cache import SurfaceCache
from autobskan.gui.layout import build_layout
from autobskan.image import AR, lwfplot, post_processing, stmplot

LOGGER = logging.getLogger("autobskan.gui")
if not LOGGER.handlers:
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s %(name)s: %(message)s",
    )
for noisy_logger in ("choreographer", "kaleido"):
    logging.getLogger(noisy_logger).setLevel(logging.WARNING)


COLORMAP_PRESETS = [
    "afmhot",
    "inferno",
    "magma",
    "plasma",
    "viridis",
    "cividis",
    "turbo",
    "RdBu_r",
    "coolwarm",
    "Greys",
]

DEFAULTS = {
    "simulation": "STM",
    "input_source": "BSKAN",
    "mode": "CONSTANT CURRENT",
    "fermi_level": None,
    "fit_radius": 0.5,
    "colormap_stm": "afmhot",
    "colormap_lwf": "RdBu_r",
    "color_range_mode": "AUTO",
    "brightness": 0.0,
    "contrast": 0.0,
    "layers": 1,
    "view_options": ["show_atoms", "show_repeated"],
    "repeat_x": 2,
    "repeat_y": 2,
    "scalebar": 5.0,
    "atom_radius_type": "ATOMIC",
    "atom_radius_scale": 10.0,
    "gaussian_sigma": 5.0,
    "xy_upsample": 1,
    "display_quality": "BALANCED",
    "grid_line_width": 1.5,
    "grid_line_color": "#d5dee8",
    "export_format": "png",
    "export_width": 4096,
    "export_height": 2048,
    "colorbar_width": 1600,
    "colorbar_height": 80,
    "line_profile_enabled": True,
    "line_magnet_enabled": False,
    "fft_enabled": True,
}

DISPLAY_PIXEL_BUDGETS = {
    "FAST": 200_000,
    "BALANCED": 350_000,
    "FULL": None,
}
_SURFACE_CACHE = SurfaceCache(
    max_entries=32,
    max_bytes=512 * 1024 * 1024,
    ttl_seconds=2 * 60 * 60,
)
MAX_EXPORT_DIMENSION = 12_000
MAX_EXPORT_PIXELS = 60_000_000
EXPORT_FORMATS = {
    "png": ("png", "image/png"),
    "jpeg": ("jpg", "image/jpeg"),
    "webp": ("webp", "image/webp"),
    "svg": ("svg", "image/svg+xml"),
    "pdf": ("pdf", "application/pdf"),
}


class _BrowserLifecycle:
    def __init__(self, enabled=False, shutdown_callback=None, grace_seconds=2.0):
        self.enabled = bool(enabled)
        self.shutdown_callback = shutdown_callback or _thread.interrupt_main
        self.grace_seconds = max(float(grace_seconds), 0.0)
        self._active = {}
        self._closing = set()
        self._ever_connected = False
        self._timer = None
        self._lock = threading.RLock()

    @property
    def active_count(self):
        with self._lock:
            return len(self._active)

    def register(self, client_id, connection_id=None):
        if not self.enabled:
            return 0
        connection_key = connection_id or "presence"
        with self._lock:
            self._ever_connected = True
            if client_id in self._closing:
                return len(self._active)
            self._active.setdefault(client_id, set()).add(connection_key)
            timer = self._timer
            self._timer = None
        if timer is not None:
            timer.cancel()
        return self.active_count

    def close(self, client_id, connection_id=None):
        if not self.enabled:
            return 0
        with self._lock:
            self._ever_connected = True
            if connection_id is None:
                self._closing.add(client_id)
                self._active.pop(client_id, None)
            else:
                connections = self._active.get(client_id)
                if connections is not None:
                    connections.discard(connection_id)
                    if not connections:
                        self._active.pop(client_id, None)
            remaining = len(self._active)
            if remaining == 0:
                self._schedule_shutdown_locked()
        return remaining

    def _schedule_shutdown_locked(self):
        if not self._ever_connected:
            return
        previous = self._timer
        timer = threading.Timer(self.grace_seconds, self._shutdown_if_unused)
        timer.daemon = True
        self._timer = timer
        if previous is not None:
            previous.cancel()
        timer.start()

    def _shutdown_if_unused(self):
        with self._lock:
            if self._active or not self._ever_connected:
                return
            self._timer = None
        _terminal_log("Last GUI window closed; stopping AutoBSKAN.")
        self.shutdown_callback()


def _gui_version_info():
    package_version = getattr(autobskan, "__version__", None)
    if not package_version:
        try:
            package_version = importlib_metadata.version("autobskan")
        except importlib_metadata.PackageNotFoundError:
            package_version = "editable"

    return {"package_version": package_version}


def _static_export_help():
    try:
        kaleido_version = importlib_metadata.version("kaleido")
    except importlib_metadata.PackageNotFoundError:
        return (
            "Kaleido is missing from this Python environment. Re-run "
            "`pip install -r requirements.txt` and `pip install -e .`."
        )
    return (
        f"Kaleido {kaleido_version} is installed. Kaleido 1 also requires "
        "Chrome or Chromium; run `plotly_get_chrome` if no browser is available."
    )


def _safe_path(path_value):
    if isinstance(path_value, str) and path_value.strip():
        stripped = path_value.strip().strip('"').strip("'")
        candidate = os.path.abspath(os.path.expanduser(stripped))
        if os.path.exists(candidate):
            return candidate
    return None


def _human_size(num_bytes):
    size = float(num_bytes)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024.0 or unit == "TB":
            return f"{size:.2f} {unit}"
        size /= 1024.0
    return f"{float(num_bytes):.2f} B"


def _path_with_size(path):
    try:
        return f"{path} ({_human_size(os.path.getsize(path))})"
    except OSError:
        return path


def _effective_path(primary_value, fallback_value=None):
    if isinstance(primary_value, str) and primary_value.strip():
        return _safe_path(primary_value)
    return _safe_path(fallback_value)


def _bskanin_export_text(
    *,
    simulation,
    input_source,
    volume_path,
    mode,
    iso_value,
    fit_radius,
    cmap_name,
    contrast,
    brightness,
    poscar_path,
    iteration,
    gaussian_sigma,
    gamma,
    display_atoms,
    layers,
    radius_type,
    size_ratio,
    display_cell,
    fermi_level=None,
):
    """Serialize current GUI settings as concise CLI input with absolute paths."""

    source = str(input_source).strip().upper()
    if source not in {"BSKAN", "VASP"}:
        raise ValueError(f"Unsupported input source: {input_source}")

    resolved_volume = _safe_path(volume_path)
    if resolved_volume is None or not os.path.isfile(resolved_volume):
        raise ValueError("Select an accessible volumetric file before exporting bskan.in.")

    resolved_poscar = None
    if isinstance(poscar_path, str) and poscar_path.strip():
        resolved_poscar = _safe_path(poscar_path)
        if resolved_poscar is None or not os.path.isfile(resolved_poscar):
            raise ValueError("The selected structure file is not accessible.")
    if (display_atoms or display_cell) and resolved_poscar is None:
        raise ValueError(
            "Atom or cell overlays require an accessible structure file in CLI exports."
        )

    nx, ny = (int(iteration[0]), int(iteration[1]))
    volume_keyword = "CURRENT" if source == "BSKAN" else "VOLUME"
    lines = [
        "TASK = IMAGE",
        f"SIMULATION = {simulation}",
        f"INPUT_SOURCE = {source}",
        f"{volume_keyword} = {resolved_volume}",
        "ISO_AUTO = FALSE",
        f"ISO = {float(iso_value):.10g}",
    ]
    if simulation == "PHI_APP":
        lines.append(f"FIT_RADIUS = {float(fit_radius):.10g}")
    if fermi_level is not None:
        lines.append(f"FERMI_LEVEL = {float(fermi_level):.10g}")

    lines.extend(
        [
            f"CMAP = {cmap_name}",
            f"CONTRAST = {float(contrast):.6g}",
            f"BRIGHTNESS = {float(brightness):.6g}",
            f"IMAGE_MODE = {mode}",
        ]
    )
    if resolved_poscar is not None:
        lines.append(f"POSCAR = {resolved_poscar}")
    else:
        lines.append(f"GAMMA = {float(gamma):.6f}")
    lines.extend(
        [
            f"ITERATION = {nx}, {ny}",
            f"BLUR_SIGMA = {float(gaussian_sigma):.6g}",
            f"DISPLAY_ATOMS = {'TRUE' if display_atoms else 'FALSE'}",
            f"DISPLAY_CELL = {'TRUE' if display_cell else 'FALSE'}",
        ]
    )
    if display_atoms:
        lines.extend(
            [
                f"LAYERS = {max(int(layers), 1)}",
                f"RADIUS_TYPE = {radius_type}",
                f"SIZE_RATIO = {int(round(float(size_ratio)))}",
            ]
        )
    return "\n".join(lines) + "\n"


_UPLOAD_SESSION = None
_CHUNK_UPLOADS = {}
_CHUNK_UPLOAD_LOCK = threading.RLock()
_CHUNK_UPLOAD_LIMIT = 8 * 1024 * 1024


def _managed_upload_dir():
    global _UPLOAD_SESSION
    if _UPLOAD_SESSION is None:
        _UPLOAD_SESSION = tempfile.TemporaryDirectory(prefix="autobskan-gui-")
    return _UPLOAD_SESSION.name


def _is_managed_upload(path):
    if _UPLOAD_SESSION is None or not isinstance(path, str):
        return False
    try:
        return os.path.commonpath(
            [os.path.realpath(path), os.path.realpath(_UPLOAD_SESSION.name)]
        ) == os.path.realpath(_UPLOAD_SESSION.name)
    except (OSError, ValueError):
        return False


def _file_matches_payload(path, payload, chunk_size=1024 * 1024):
    try:
        if not os.path.isfile(path) or os.path.getsize(path) != len(payload):
            return False
        offset = 0
        view = memoryview(payload)
        with open(path, "rb") as fileobj:
            while offset < len(payload):
                chunk = fileobj.read(min(chunk_size, len(payload) - offset))
                if not chunk or chunk != view[offset : offset + len(chunk)]:
                    return False
                offset += len(chunk)
        return True
    except OSError:
        return False


def _save_upload(contents, filename, prefix):
    if isinstance(contents, list):
        contents = contents[0] if contents else None
    if isinstance(filename, list):
        filename = filename[0] if filename else None
    if not isinstance(contents, str) or "," not in contents:
        return None
    _, encoded = contents.split(",", 1)
    try:
        payload = base64.b64decode(encoded, validate=True)
    except (binascii.Error, ValueError):
        return None

    raw_name = os.path.basename(filename or f"{prefix}.dat")
    safe_name = re.sub(r"[^A-Za-z0-9._-]", "_", raw_name)

    # Browsers do not expose the selected absolute path. When the exact same
    # file exists in the launch directory, reuse it instead of making a copy.
    local_candidate = os.path.abspath(raw_name)
    if _file_matches_payload(local_candidate, payload):
        return local_candidate

    upload_dir = _managed_upload_dir()
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    path = os.path.join(upload_dir, f"{prefix}_{stamp}_{safe_name}")
    with open(path, "wb") as fileobj:
        fileobj.write(payload)
    return path


@lru_cache(maxsize=16)
def _load_current_cached(path, mtime_ns, size_bytes):
    return stmplot.Current(path)


def _load_current(path):
    stat = os.stat(path)
    return _load_current_cached(path, stat.st_mtime_ns, stat.st_size)


_LOCPOT_LOAD_LOCK = threading.RLock()


@lru_cache(maxsize=2)
def _load_locpot_cached(path, mtime_ns, size_bytes):
    return AR.Chgcar(path)


def _load_locpot(path):
    stat = os.stat(path)
    cache_key = (path, stat.st_mtime_ns, stat.st_size)
    # Dash may fire parse, range, and render callbacks together. lru_cache
    # protects its dictionary but does not coalesce concurrent cache misses.
    with _LOCPOT_LOAD_LOCK:
        return _load_locpot_cached(*cache_key)


@lru_cache(maxsize=16)
def _load_fermi_energy_cached(path, mtime_ns, size_bytes):
    return lwfplot.fermi_energy_from_outcar(path)


def _resolve_lwf_fermi_level(volume_path, manual_value=None):
    if manual_value is not None and str(manual_value).strip().upper() not in {
        "",
        "AUTO",
        "NONE",
    }:
        value = float(manual_value)
        if not np.isfinite(value):
            raise ValueError("Fermi level should be a finite value in eV.")
        return value, "manual input"

    outcar_path = lwfplot.sibling_outcar_path(volume_path)
    if outcar_path is None:
        raise ValueError(
            "Phi_loc requires E_F. Enter the Fermi level in eV, or place OUTCAR "
            "in the same directory as LOCPOT."
        )
    stat = os.stat(outcar_path)
    value = _load_fermi_energy_cached(
        outcar_path,
        stat.st_mtime_ns,
        stat.st_size,
    )
    return value, outcar_path


@lru_cache(maxsize=16)
def _load_structure_cached(path, mtime_ns, size_bytes):
    return AR.to_new_cell(AR.read_vasp_robust(path))


def _load_structure(path):
    stat = os.stat(path)
    return _load_structure_cached(path, stat.st_mtime_ns, stat.st_size)


def _embedded_structure_from_volume(volumetric):
    atoms = getattr(volumetric, "atoms", None)
    if atoms is None:
        raise ValueError("The VASP volumetric header does not contain a structure.")
    return AR.to_new_cell(atoms.copy())


def _embedded_structure_metadata(volumetric, path):
    structure = _embedded_structure_from_volume(volumetric)
    a, b, c, alpha, beta, gamma = structure.cell.cellpar()
    return {
        "active": True,
        "volume_path": str(path),
        "volume_kind": _vasp_volume_kind(path),
        "atoms": int(len(structure)),
        "cell": [float(a), float(b), float(c)],
        "angles": [float(alpha), float(beta), float(gamma)],
    }


def _structure_override_is_confirmed(structure_path, decision, volume_path=None):
    if structure_path is None or not isinstance(decision, dict):
        return False
    return bool(
        decision.get("path") == structure_path
        and decision.get("decision") == "proceed"
        and (
            volume_path is None
            or decision.get("volume_path") == volume_path
        )
    )


def _parse_repeat(nx, ny):
    try:
        nx_i = max(int(nx), 1)
        ny_i = max(int(ny), 1)
    except Exception:
        nx_i, ny_i = 1, 1
    return nx_i, ny_i


def _brightness_limits(z_data, brightness):
    z = np.asarray(z_data, dtype=float)
    finite = z[np.isfinite(z)]
    if finite.size == 0:
        return 0.0, 1.0
    z_min = float(np.min(finite))
    z_max = float(np.max(finite))
    span = max(z_max - z_min, 1e-12)
    if brightness < 0:
        lower = z_min
        upper = z_max + (-brightness) * span
    else:
        lower = z_min - brightness * span
        upper = z_max
    if np.isclose(lower, upper):
        upper = lower + 1e-9
    return lower, upper


def _manual_color_limits(mode, vmin, vmax):
    if str(mode).upper() != "MANUAL":
        return None
    try:
        lower = float(vmin)
        upper = float(vmax)
    except (TypeError, ValueError):
        return None
    if not np.isfinite(lower) or not np.isfinite(upper) or lower >= upper:
        return None
    return lower, upper


def _effective_color_limits(
    z_data,
    brightness,
    mode="AUTO",
    vmin=None,
    vmax=None,
):
    manual = _manual_color_limits(mode, vmin, vmax)
    if manual is not None:
        return manual
    effective_brightness = 0.0 if str(mode).upper() == "MANUAL" else brightness
    return _brightness_limits(z_data, float(effective_brightness or 0.0))


def _colorscale(cmap_name, contrast, npoint=256):
    cmap_name = cmap_name if cmap_name in colormaps else "afmhot"
    cmap = colormaps[cmap_name]
    gamma = max(0.05, 1.0 + float(contrast))
    t = np.linspace(0.0, 1.0, npoint)
    rgba = cmap(np.power(t, gamma))
    scale = [
        [float(i) / float(npoint - 1), mcolors.to_hex(rgba[i, :3])]
        for i in range(npoint)
    ]
    return scale


def _grid_axes(X, Y, Z):
    x_min, x_max = float(np.nanmin(X)), float(np.nanmax(X))
    y_min, y_max = float(np.nanmin(Y)), float(np.nanmax(Y))
    x_axis = np.linspace(x_min, x_max, Z.shape[1])
    y_axis = np.linspace(y_min, y_max, Z.shape[0])
    return x_axis, y_axis, x_min, x_max, y_min, y_max


def _regularize_scalar_grid(X, Y, Z, x_axis, y_axis):
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    Z = np.asarray(Z, dtype=float)
    X_regular, Y_regular = np.meshgrid(
        np.asarray(x_axis, dtype=float),
        np.asarray(y_axis, dtype=float),
    )
    points = np.column_stack([X.ravel(), Y.ravel()])
    values = Z.ravel()
    finite = np.isfinite(points[:, 0]) & np.isfinite(points[:, 1]) & np.isfinite(values)
    if np.count_nonzero(finite) < 3:
        return np.full_like(X_regular, np.nan, dtype=float)

    Z_regular = griddata(
        points[finite],
        values[finite],
        (X_regular, Y_regular),
        method="linear",
    )
    if not np.isfinite(Z_regular).any():
        Z_regular = griddata(
            points[finite],
            values[finite],
            (X_regular, Y_regular),
            method="nearest",
        )
    return np.asarray(Z_regular, dtype=float)


def _rectilinear_surface(X, Y, Z):
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    Z = np.asarray(Z, dtype=float)
    x_axis, y_axis, *_bounds = _rectilinear_axes_from_grid(X, Y, Z)
    if _is_rectilinear_grid(X, Y):
        return np.asarray(x_axis), np.asarray(y_axis), Z
    return (
        np.asarray(x_axis),
        np.asarray(y_axis),
        _regularize_scalar_grid(X, Y, Z, x_axis, y_axis),
    )


def _display_surface(x_axis, y_axis, Z, quality="BALANCED"):
    x_axis = np.asarray(x_axis, dtype=float)
    y_axis = np.asarray(y_axis, dtype=float)
    Z = np.asarray(Z, dtype=float)
    budget = DISPLAY_PIXEL_BUDGETS.get(
        str(quality or "BALANCED").upper(),
        DISPLAY_PIXEL_BUDGETS["BALANCED"],
    )
    if budget is None or Z.size <= budget:
        return x_axis, y_axis, Z

    scale = np.sqrt(float(budget) / float(Z.size))
    rows = max(2, int(round(Z.shape[0] * scale)))
    columns = max(2, int(round(Z.shape[1] * scale)))
    if rows * columns > budget:
        columns = max(2, int(budget // rows))
    display_x = np.linspace(float(x_axis[0]), float(x_axis[-1]), columns)
    display_y = np.linspace(float(y_axis[0]), float(y_axis[-1]), rows)
    target_x, target_y = np.meshgrid(display_x, display_y)
    interpolator = RegularGridInterpolator(
        (y_axis, x_axis),
        Z,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )
    display_z = interpolator((target_y, target_x))
    return display_x, display_y, np.asarray(display_z, dtype=float)


def _store_surface_data(data):
    return _SURFACE_CACHE.put(data)


def _get_surface_data(key):
    return _SURFACE_CACHE.get(key)


def _surface_analysis_values(surface):
    return np.asarray(surface.get("analysis_z", surface["z"]), dtype=float)


def _is_rectilinear_grid(X, Y, tolerance=1e-8):
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    if X.ndim != 2 or Y.ndim != 2 or X.shape != Y.shape:
        return False
    x_ref = X[0, :]
    y_ref = Y[:, 0]
    return bool(
        np.allclose(X, x_ref.reshape(1, -1), atol=tolerance, rtol=0.0)
        and np.allclose(Y, y_ref.reshape(-1, 1), atol=tolerance, rtol=0.0)
    )


def _rectilinear_axes_from_grid(X, Y, Z):
    if _is_rectilinear_grid(X, Y):
        return (
            np.asarray(X[0, :], dtype=float),
            np.asarray(Y[:, 0], dtype=float),
            float(np.nanmin(X)),
            float(np.nanmax(X)),
            float(np.nanmin(Y)),
            float(np.nanmax(Y)),
        )
    return _grid_axes(X, Y, Z)


def _scalar_field_figure(
    X,
    Y,
    Z,
    colorscale,
    z_min,
    z_max,
    show_blur,
    gaussian_sigma,
):
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    Z = np.asarray(Z, dtype=float)
    x_axis, y_axis, x_min, x_max, y_min, y_max = _rectilinear_axes_from_grid(X, Y, Z)
    use_carpet = not _is_rectilinear_grid(X, Y)

    if use_carpet:
        Z_display = _regularize_scalar_grid(X, Y, Z, x_axis, y_axis)
    else:
        Z_display = Z

    blur_applied = bool(show_blur and float(gaussian_sigma) > 0)
    if blur_applied:
        Z_display = _gaussian_blur_surface(Z_display, gaussian_sigma)

    fig = go.Figure(
        go.Heatmap(
            x=x_axis,
            y=y_axis,
            z=Z_display,
            colorscale=colorscale,
            zmin=z_min,
            zmax=z_max,
            zsmooth=False,
            opacity=1.0,
            showscale=False,
            hovertemplate="x=%{x:.3f}<br>y=%{y:.3f}<br>value=%{z:.6g}<extra></extra>",
            meta={"autobskan_role": "scalar-field"},
        )
    )
    return (
        fig,
        x_axis,
        y_axis,
        x_min,
        x_max,
        y_min,
        y_max,
        blur_applied,
    )


def _lwf_cell_c_length(locpot):
    return lwfplot.cell_c_length(locpot)


def _lwf_topmost_z(structure=None, locpot=None):
    return lwfplot.topmost_z(structure, locpot)


def _vasp_volume_kind(path):
    if not isinstance(path, str):
        return "VASP volume"
    basename = os.path.basename(path).upper()
    if "PARCHG" in basename:
        return "PARCHG"
    if "LOCPOT" in basename:
        return "LOCPOT"
    if "CHGCAR" in basename:
        return "CHGCAR"
    return "VASP volume"


def _named_volume_source(path):
    if not isinstance(path, str):
        return None
    name = os.path.basename(path).upper()
    if any(kind in name for kind in ("PARCHG", "LOCPOT", "CHGCAR")):
        return "VASP"
    if re.search(r"(?:^|[_.-])CURRENT(?:$|[_.-])", name):
        return "BSKAN"
    return None


def _effective_input_source(input_source, _volume_path=None):
    source = str(input_source or "").strip().upper()
    if source not in {"BSKAN", "VASP"}:
        raise ValueError("Select either bSKAN or VASP as the data source.")
    return source


def _validate_volume_source(input_source, path):
    source = _effective_input_source(input_source, path)
    named_source = _named_volume_source(path)
    if named_source is not None and named_source != source:
        filename = os.path.basename(path)
        expected = "bSKAN" if named_source == "BSKAN" else "VASP"
        selected = "bSKAN" if source == "BSKAN" else "VASP"
        raise ValueError(
            f"{filename} is identified as {expected} data, but {selected} is selected."
        )
    return source


def _validate_simulation_input(simulation, input_source, path):
    simulation = str(simulation or "").strip().upper()
    source = _validate_volume_source(input_source, path)
    kind = _vasp_volume_kind(path) if source == "VASP" else "CURRENT"

    if simulation == "PHI_APP":
        if source == "VASP" and kind != "PARCHG":
            raise ValueError(
                "Apparent barrier height (Phi_app) requires PARCHG for VASP data. "
                "LOCPOT is reserved for Phi_loc."
            )
        return source

    if simulation == "LWF":
        if source != "VASP":
            raise ValueError("Phi_loc requires VASP LOCPOT data, not bSKAN CURRENT.")
        if kind != "LOCPOT":
            raise ValueError(
                "Phi_loc requires a LOCPOT file and is evaluated as "
                "V_LOCPOT(x,y,z0) - E_F at constant height."
            )
        return source

    if simulation != "STM":
        raise ValueError(f"Unsupported simulation type: {simulation}")
    return source


def _simulation_display_name(simulation):
    if simulation == "PHI_APP":
        return "Phi_app"
    if simulation == "LWF":
        return "Phi_loc"
    return str(simulation or "STM")


def _vasp_constant_height_label(path):
    kind = _vasp_volume_kind(path)
    if kind == "PARCHG":
        return "PARCHG density", "arb."
    if kind == "LOCPOT":
        return "Potential slice (eV)", "eV"
    if kind == "CHGCAR":
        return "Charge density", "arb."
    return "Volumetric value", "arb."


def _surface_stm(current, mode, iso_value):
    constant_current = mode == "CONSTANT CURRENT"
    z_surface = stmplot.get_surface_bskan(
        current=current,
        isosurface=iso_value,
        constant_current=constant_current,
        return_index=False,
    )
    z_surface = np.asarray(z_surface, dtype=float)
    finite_ratio = float(np.count_nonzero(np.isfinite(z_surface))) / float(
        z_surface.size
    )
    if finite_ratio < 0.95:
        _terminal_log(
            "Surface extraction produced non-finite values with exponential interpolation. "
            "Retrying with linear interpolation.",
            level="warning",
        )
        z_surface = stmplot.get_surface_bskan(
            current=current,
            isosurface=iso_value,
            constant_current=constant_current,
            return_index=False,
            interpolate="linear",
        )
        z_surface = np.asarray(z_surface, dtype=float)
        if not np.any(np.isfinite(z_surface)):
            raise ValueError(
                "Surface extraction failed: all values are non-finite "
                "(check CURRENT data and isosurface range)."
            )
    return z_surface


def _surface_lwf(locpot, height, fermi_level, topmost=0.0):
    return lwfplot.local_work_function_slice(
        locpot,
        height,
        fermi_level=fermi_level,
        topmost=topmost,
    )


def _surface_apparent_barrier_from_parchg(
    parchg,
    height,
    topmost=0.0,
    fit_radius=0.5,
):
    return lwfplot.apparent_barrier_from_vasp_density(
        parchg,
        height,
        topmost=topmost,
        fit_radius=fit_radius,
    )


def _prepare_stm_grid(current, structure, z_raw, show_repeat, nx, ny, gamma):
    if show_repeat:
        if structure is not None:
            result = post_processing.array_iter(
                z_raw,
                nx=nx,
                ny=ny,
                cell=structure.cell,
                input_filetype="BSKAN",
                diagonal_transform_for_hexagonal=True,
            )
        else:
            result = post_processing.array_iter(
                z_raw,
                nx=nx,
                ny=ny,
                cell=current.cell,
                input_filetype="BSKAN",
                gamma=gamma,
                l_a=float(current.cellpar[0]),
                l_b=float(current.cellpar[1]),
            )

        if isinstance(result, tuple) and len(result) == 3:
            X, Y, Z = result
            return np.asarray(X), np.asarray(Y), np.asarray(Z)

        # Fallback for branches of array_iter() that return only Z.
        Z = np.asarray(result)
        x = np.linspace(0.0, float(current.cellpar[0]) * nx, Z.shape[1])
        y = np.linspace(0.0, float(current.cellpar[1]) * ny, Z.shape[0])
        X, Y = np.meshgrid(x, y)
        return X, Y, Z

    x = np.linspace(0.0, float(current.cellpar[0]), current.grids[0])
    y = np.linspace(0.0, float(current.cellpar[1]), current.grids[1])
    X, Y = np.meshgrid(x, y)
    return X, Y, np.asarray(z_raw)


def _prepare_lwf_grid(locpot, structure, z_raw, show_repeat, nx, ny):
    if structure is not None:
        cell = AR.to_new_cell_onlycell(np.asarray(structure.cell))
    else:
        cell = AR.to_new_cell_onlycell(np.asarray(locpot.cell))

    z_base = np.asarray(z_raw, dtype=float)
    rep_x = max(int(nx), 1) if show_repeat else 1
    rep_y = max(int(ny), 1) if show_repeat else 1
    X, Y, Z = post_processing.array_iter(
        z_base,
        nx=rep_x,
        ny=rep_y,
        cell=cell,
        input_filetype="VASP",
        filled_array_for_contourf=False,
        verbose=False,
    )
    return np.asarray(X), np.asarray(Y), np.asarray(Z)


def _upsample_surface_xy(X, Y, Z, factor):
    try:
        n = max(int(round(float(factor))), 1)
    except Exception:
        n = 1
    if n <= 1:
        return (
            np.asarray(X, dtype=float),
            np.asarray(Y, dtype=float),
            np.asarray(Z, dtype=float),
            1,
        )
    X_up = ndimage.zoom(np.asarray(X, dtype=float), zoom=n, order=1)
    Y_up = ndimage.zoom(np.asarray(Y, dtype=float), zoom=n, order=1)
    Z_up = ndimage.zoom(np.asarray(Z, dtype=float), zoom=n, order=1)
    return X_up, Y_up, Z_up, n


def _overlay_traces(
    structure,
    x_min,
    x_max,
    y_min,
    y_max,
    show_atoms,
    show_grid,
    layers,
    radius_type,
    radius_scale,
    grid_line_width,
    grid_line_color,
):
    traces = []
    if structure is None:
        return traces

    if show_grid:
        fig, ax = plt.subplots()
        stmplot.plot_cell(
            structure,
            xmin=x_min,
            xmax=x_max,
            ymin=y_min,
            ymax=y_max,
            ax=ax,
            ls=":",
            c=grid_line_color,
            lw=float(grid_line_width),
        )
        x_line, y_line = [], []
        for line in ax.lines:
            x_data = np.asarray(line.get_xdata(), dtype=float).tolist()
            y_data = np.asarray(line.get_ydata(), dtype=float).tolist()
            x_line.extend(x_data + [None])
            y_line.extend(y_data + [None])
        plt.close(fig)
        if x_line:
            traces.append(
                go.Scatter(
                    x=x_line,
                    y=y_line,
                    mode="lines",
                    line={
                        "color": grid_line_color,
                        "width": float(grid_line_width),
                        "dash": "dot",
                    },
                    hoverinfo="skip",
                    cliponaxis=True,
                    showlegend=False,
                )
            )

    if show_atoms:
        fig, ax = plt.subplots()
        stmplot.plot_atoms(
            structure,
            xmin=x_min,
            xmax=x_max,
            ymin=y_min,
            ymax=y_max,
            ax=ax,
            top_n_layers=max(int(layers), 1),
            atoms_info="VESTA",
            radii_type=radius_type,
            radii_marker_scale=float(radius_scale),
        )
        grouped = {}
        for line in ax.lines:
            x_data = np.asarray(line.get_xdata(), dtype=float)
            y_data = np.asarray(line.get_ydata(), dtype=float)
            if x_data.size != 1 or y_data.size != 1:
                continue
            color = mcolors.to_hex(line.get_color())
            marker_size = float(line.get_markersize())
            key = (color, marker_size)
            if key not in grouped:
                grouped[key] = [[], []]
            grouped[key][0].append(float(x_data[0]))
            grouped[key][1].append(float(y_data[0]))
        plt.close(fig)
        for (color, marker_size), (xs, ys) in grouped.items():
            traces.append(
                go.Scatter(
                    x=xs,
                    y=ys,
                    mode="markers",
                    marker={
                        "size": marker_size,
                        "color": color,
                        "line": {"width": 0},
                    },
                    hoverinfo="skip",
                    cliponaxis=True,
                    showlegend=False,
                )
            )
    return traces


def _gaussian_blur_surface(Z, sigma):
    sigma = float(sigma)
    if sigma <= 0:
        return np.asarray(Z, dtype=float)

    z = np.asarray(Z, dtype=float)
    finite = np.isfinite(z)
    if not np.any(finite):
        return z

    weights = finite.astype(float)
    values = np.where(finite, z, 0.0)
    blurred_values = ndimage.gaussian_filter(
        values,
        sigma=sigma,
        order=0,
        mode="wrap",
    )
    blurred_weights = ndimage.gaussian_filter(
        weights,
        sigma=sigma,
        order=0,
        mode="wrap",
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        result = np.divide(
            blurred_values,
            blurred_weights,
            out=np.full_like(blurred_values, np.nan),
            where=blurred_weights > 1e-12,
        )
    return np.where(finite, result, np.nan)


def _colorbar_figure(colorscale, zmin, zmax, title):
    z_line = np.linspace(float(zmin), float(zmax), 400).reshape(1, -1)
    fig = go.Figure(
        go.Heatmap(
            z=z_line,
            colorscale=colorscale,
            zmin=float(zmin),
            zmax=float(zmax),
            opacity=0.0,
            colorbar={
                "orientation": "h",
                "x": 0.5,
                "xanchor": "center",
                "y": 0.5,
                "len": 0.94,
                "thickness": 20,
                "title": title,
            },
            hoverinfo="skip",
            showscale=True,
            meta={"autobskan_role": "colorbar-field"},
        )
    )
    fig.update_layout(
        template="plotly_white",
        title="Color scale",
        height=145,
        margin={"l": 24, "r": 24, "t": 34, "b": 18},
    )
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    return fig


def _export_colorbar_figure(colorscale, zmin, zmax):
    z_line = np.linspace(float(zmin), float(zmax), 1024)
    fig = go.Figure(
        go.Heatmap(
            z=np.vstack([z_line, z_line]),
            colorscale=colorscale,
            zmin=float(zmin),
            zmax=float(zmax),
            hoverinfo="skip",
            showscale=False,
            zsmooth=False,
        )
    )
    fig.update_layout(
        template="plotly_white",
        margin={"l": 0, "r": 0, "t": 0, "b": 0, "pad": 0},
        paper_bgcolor="#ffffff",
        plot_bgcolor="#ffffff",
        showlegend=False,
    )
    fig.update_xaxes(visible=False, fixedrange=True, domain=[0.0, 1.0])
    fig.update_yaxes(visible=False, fixedrange=True, domain=[0.0, 1.0])
    return fig


def _normalise_export_format(value):
    normalized = str(value or "png").strip().lower()
    if normalized == "jpg":
        normalized = "jpeg"
    if normalized not in EXPORT_FORMATS:
        raise ValueError(
            "Output format should be PNG, JPEG, WebP, SVG, or PDF."
        )
    return normalized


def _validate_export_dimensions(width, height):
    try:
        width_px = int(round(float(width)))
        height_px = int(round(float(height)))
    except (TypeError, ValueError) as exc:
        raise ValueError("Export width and height should be integer pixels.") from exc
    if not (256 <= width_px <= MAX_EXPORT_DIMENSION):
        raise ValueError(
            f"Export width should be between 256 and {MAX_EXPORT_DIMENSION} px."
        )
    if not (16 <= height_px <= MAX_EXPORT_DIMENSION):
        raise ValueError(
            f"Export height should be between 16 and {MAX_EXPORT_DIMENSION} px."
        )
    if width_px * height_px > MAX_EXPORT_PIXELS:
        raise ValueError(
            f"Export is limited to {MAX_EXPORT_PIXELS:,} pixels."
        )
    return width_px, height_px


def _export_file_details(value, stem):
    image_format = _normalise_export_format(value)
    extension, mime_type = EXPORT_FORMATS[image_format]
    return image_format, f"{stem}.{extension}", mime_type


def _render_export_figure(figure, image_format, width, height, opaque=False):
    width_px, height_px = _validate_export_dimensions(width, height)
    image_format = _normalise_export_format(image_format)
    if opaque or image_format in {"jpeg", "pdf"}:
        figure.update_layout(
            paper_bgcolor="#ffffff",
            plot_bgcolor="#ffffff",
        )
    return pio.to_image(
        figure,
        format=image_format,
        width=width_px,
        height=height_px,
        scale=1,
    )


def _cached_surface_figure(
    surface,
    colorscale,
    z_min,
    z_max,
    cmap_name,
    contrast,
    line_points=None,
):
    x_axis = np.asarray(surface["x_axis"], dtype=float)
    y_axis = np.asarray(surface["y_axis"], dtype=float)
    has_analysis_surface = "analysis_z" in surface
    Z = np.asarray(surface.get("analysis_z", surface["z"]), dtype=float)
    X, Y = np.meshgrid(x_axis, y_axis)
    figure, *_axes = _scalar_field_figure(
        X=X,
        Y=Y,
        Z=Z,
        colorscale=colorscale,
        z_min=float(z_min),
        z_max=float(z_max),
        show_blur=bool(surface.get("blur_applied")) and not has_analysis_surface,
        gaussian_sigma=float(surface.get("gaussian_sigma", 0.0)),
    )
    for trace in surface.get("overlay", []):
        figure.add_trace(trace)
    for shape in surface.get("shapes", []):
        figure.add_shape(shape)
    _add_line_selection_traces(
        figure,
        _normalise_line_points(line_points),
    )
    return figure


def _snap_point_to_local_optimum(X, Y, Z, point, window_radius=8):
    x0, y0 = float(point[0]), float(point[1])
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    Z = np.asarray(Z, dtype=float)
    rectilinear_axes = X.ndim == 1 and Y.ndim == 1 and Z.shape == (Y.size, X.size)
    if rectilinear_axes:
        j0 = int(np.nanargmin(np.abs(X - x0)))
        i0 = int(np.nanargmin(np.abs(Y - y0)))
    else:
        if X.shape != Y.shape or X.shape != Z.shape:
            return [x0, y0]
        d2 = (X - x0) ** 2 + (Y - y0) ** 2
        d2 = np.where(np.isfinite(d2), d2, np.inf)
        if not np.isfinite(np.min(d2)):
            return [x0, y0]
        i0, j0 = np.unravel_index(np.argmin(d2), d2.shape)
    i1 = max(i0 - int(window_radius), 0)
    i2 = min(i0 + int(window_radius) + 1, Z.shape[0])
    j1 = max(j0 - int(window_radius), 0)
    j2 = min(j0 + int(window_radius) + 1, Z.shape[1])

    sub = Z[i1:i2, j1:j2]
    if sub.size == 0:
        if rectilinear_axes:
            return [float(X[j0]), float(Y[i0])]
        return [float(X[i0, j0]), float(Y[i0, j0])]
    finite = np.isfinite(sub)
    if not np.any(finite):
        if rectilinear_axes:
            return [float(X[j0]), float(Y[i0])]
        return [float(X[i0, j0]), float(Y[i0, j0])]

    ci, cj = i0 - i1, j0 - j1
    center = (
        float(sub[ci, cj]) if np.isfinite(sub[ci, cj]) else float(np.nanmedian(sub))
    )
    median_val = float(np.nanmedian(sub))
    search_minimum = center <= median_val

    masked = np.where(finite, sub, np.nanmedian(sub))
    filt = (
        ndimage.minimum_filter(masked, size=3, mode="nearest")
        if search_minimum
        else ndimage.maximum_filter(masked, size=3, mode="nearest")
    )
    candidates_mask = np.isclose(masked, filt, atol=1e-12)
    cand = np.argwhere(candidates_mask)

    if cand.size == 0:
        if search_minimum:
            li, lj = np.unravel_index(np.nanargmin(masked), masked.shape)
        else:
            li, lj = np.unravel_index(np.nanargmax(masked), masked.shape)
    else:
        dist = np.sqrt((cand[:, 0] - ci) ** 2 + (cand[:, 1] - cj) ** 2)
        li, lj = cand[int(np.argmin(dist))]

    ig, jg = i1 + int(li), j1 + int(lj)
    if rectilinear_axes:
        return [float(X[jg]), float(Y[ig])]
    return [float(X[ig, jg]), float(Y[ig, jg])]


def _empty_figure(title, message, height=420):
    fig = go.Figure()
    fig.add_annotation(
        text=message,
        x=0.5,
        y=0.5,
        xref="paper",
        yref="paper",
        showarrow=False,
        font={"size": 14, "color": "#5b6976"},
    )
    fig.update_layout(
        template="plotly_white",
        title=title,
        height=height,
        margin={"l": 30, "r": 20, "t": 50, "b": 30},
    )
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    return fig


def _has_scalar_surface_trace(figure_dict):
    try:
        data = figure_dict.get("data") or []
    except AttributeError:
        return False
    for trace in data:
        if not isinstance(trace, dict):
            continue
        if trace.get("type") in {"heatmap", "contourcarpet"}:
            return True
    return False


def _plotly_array(values):
    if isinstance(values, dict) and "bdata" in values and "dtype" in values:
        array = np.frombuffer(
            base64.b64decode(values["bdata"]),
            dtype=np.dtype(values["dtype"]),
        )
        shape = values.get("shape")
        if isinstance(shape, str) and shape.strip():
            array = array.reshape(tuple(int(part) for part in shape.split(",")))
        elif isinstance(shape, (list, tuple)) and shape:
            array = array.reshape(tuple(int(part) for part in shape))
        return array
    return np.asarray(values)


def _surface_from_figure(figure_dict):
    """Return the rendered scalar field without re-running the volume analysis."""
    try:
        traces = figure_dict.get("data") or []
    except AttributeError as exc:
        raise ValueError("Rendered figure is unavailable.") from exc

    for trace in traces:
        if not isinstance(trace, dict) or trace.get("type") != "heatmap":
            continue
        z_data = trace.get("z")
        x_data = trace.get("x")
        y_data = trace.get("y")
        if z_data is None or x_data is None or y_data is None:
            continue

        Z = np.asarray(_plotly_array(z_data), dtype=float)
        x_values = np.asarray(_plotly_array(x_data), dtype=float)
        y_values = np.asarray(_plotly_array(y_data), dtype=float)
        if Z.ndim != 2:
            continue
        if x_values.ndim == 1 and y_values.ndim == 1:
            if x_values.size != Z.shape[1] or y_values.size != Z.shape[0]:
                continue
            X, Y = np.meshgrid(x_values, y_values)
            return X, Y, Z, x_values, y_values
        if x_values.shape == Z.shape and y_values.shape == Z.shape:
            x_axis, y_axis, *_ = _rectilinear_axes_from_grid(x_values, y_values, Z)
            return x_values, y_values, Z, x_axis, y_axis

    raise ValueError("Rendered scalar field is unavailable.")


def _normalise_line_points(points):
    normalised = []
    for point in list(points or [])[:2]:
        try:
            x_value = float(point[0])
            y_value = float(point[1])
        except (IndexError, TypeError, ValueError):
            continue
        if np.isfinite(x_value) and np.isfinite(y_value):
            normalised.append([x_value, y_value])
    return normalised


def _clip_line_points(points, x_axis, y_axis):
    x_values = np.asarray(x_axis, dtype=float)
    y_values = np.asarray(y_axis, dtype=float)
    if x_values.size == 0 or y_values.size == 0:
        return []
    x_min, x_max = float(np.nanmin(x_values)), float(np.nanmax(x_values))
    y_min, y_max = float(np.nanmin(y_values)), float(np.nanmax(y_values))
    return [
        [
            float(np.clip(point[0], x_min, x_max)),
            float(np.clip(point[1], y_min, y_max)),
        ]
        for point in _normalise_line_points(points)
    ]


def _add_line_selection_traces(figure, points):
    if len(points) == 2:
        figure.add_trace(
            go.Scatter(
                x=[points[0][0], points[1][0]],
                y=[points[0][1], points[1][1]],
                mode="lines",
                line={"color": "#0f766e", "width": 3},
                hoverinfo="skip",
                cliponaxis=True,
                showlegend=False,
                meta={"autobskan_role": "line-selection"},
            )
        )
    if points:
        figure.add_trace(
            go.Scatter(
                x=[point[0] for point in points],
                y=[point[1] for point in points],
                mode="markers+text",
                text=["P1", "P2"][: len(points)],
                textposition="top center",
                marker={
                    "size": 13,
                    "color": ["#0f766e", "#d97706"][: len(points)],
                    "line": {"color": "#ffffff", "width": 2},
                },
                textfont={"color": "#134e4a", "size": 13},
                hoverinfo="skip",
                cliponaxis=True,
                showlegend=False,
                meta={"autobskan_role": "line-selection"},
            )
        )


def _status_text(level, message):
    stamp = datetime.datetime.now().strftime("%H:%M:%S")
    return f"[{stamp}] {level}\n{message}"


def _status_wrap(show):
    return {"display": "block"} if show else {"display": "none"}


def _terminal_log(message, level="info"):
    line = f"[autobskan-gui] {message}"
    print(line, flush=True)
    if level == "error":
        LOGGER.error(message)
    elif level == "warning":
        LOGGER.warning(message)
    else:
        LOGGER.info(message)


def _line_profile(X, Y, Z, p1, p2, title, y_unit, y_range=None):
    if np.hypot(float(p2[0]) - float(p1[0]), float(p2[1]) - float(p1[1])) <= 1e-12:
        return _empty_figure(
            title,
            "P1 and P2 must be different points.",
            height=260,
        )

    xs = np.linspace(float(p1[0]), float(p2[0]), 600)
    ys = np.linspace(float(p1[1]), float(p2[1]), 600)
    dist = np.sqrt((xs - xs[0]) ** 2 + (ys - ys[0]) ** 2)

    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    Z = np.asarray(Z, dtype=float)
    if X.ndim == 1 and Y.ndim == 1 and Z.shape == (Y.size, X.size):
        interpolator = RegularGridInterpolator(
            (Y, X),
            Z,
            method="linear",
            bounds_error=False,
            fill_value=np.nan,
        )
        z_line = interpolator(np.column_stack([ys, xs]))
        if np.any(np.isnan(z_line)):
            nearest = RegularGridInterpolator(
                (Y, X),
                Z,
                method="nearest",
                bounds_error=False,
                fill_value=None,
            )
            z_line = np.where(
                np.isnan(z_line),
                nearest(np.column_stack([ys, xs])),
                z_line,
            )
    else:
        z_line = griddata(
            (X.ravel(), Y.ravel()),
            Z.ravel(),
            (xs, ys),
            method="linear",
        )
        if np.any(np.isnan(z_line)):
            z_near = griddata(
                (X.ravel(), Y.ravel()),
                Z.ravel(),
                (xs, ys),
                method="nearest",
            )
            z_line = np.where(np.isnan(z_line), z_near, z_line)

    if len(dist) > 3:
        spline = PchipInterpolator(dist, z_line)
        dist_fine = np.linspace(dist[0], dist[-1], 1200)
        z_fine = spline(dist_fine)
    else:
        dist_fine, z_fine = dist, z_line

    fig = go.Figure(
        go.Scatter(
            x=dist_fine,
            y=z_fine,
            mode="lines",
            line={"width": 2.2, "color": "#0b6f90"},
            name="Line profile",
        )
    )
    fig.update_layout(
        template="plotly_white",
        title=title,
        height=260,
        margin={"l": 60, "r": 20, "t": 50, "b": 40},
    )
    fig.update_xaxes(title="Distance along line (A)")
    yaxis_options = {"title": f"Value ({y_unit})"}
    if y_range is not None:
        lower, upper = y_range
        if np.isfinite(lower) and np.isfinite(upper) and lower < upper:
            yaxis_options.update(range=[float(lower), float(upper)], autorange=False)
    fig.update_yaxes(**yaxis_options)
    return fig


def _line_profile_csv(figure):
    try:
        data = go.Figure(figure).data
    except (TypeError, ValueError):
        data = ()
    profile = next(
        (
            trace
            for trace in data
            if trace.name == "Line profile"
            and "lines" in str(trace.mode or "")
        ),
        None,
    )
    if profile is None:
        raise ValueError("Select P1 and P2 before exporting a line profile.")

    def numeric_values(values):
        if isinstance(values, dict) and "bdata" in values:
            dtype = np.dtype(values.get("dtype", "f8"))
            return np.frombuffer(base64.b64decode(values["bdata"]), dtype=dtype)
        return np.asarray(values if values is not None else [], dtype=float)

    x_values = numeric_values(profile.x)
    y_values = numeric_values(profile.y)
    if x_values.size == 0 or x_values.shape != y_values.shape:
        raise ValueError("The plotted line profile does not contain valid x/y data.")

    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(["x", "y"])
    writer.writerows(zip(x_values, y_values))
    return buffer.getvalue()


def _fft_figure(Z, x_axis, y_axis):
    z = np.asarray(Z, dtype=float)
    z = np.nan_to_num(z - np.nanmean(z), nan=0.0)
    window = np.outer(np.hanning(z.shape[0]), np.hanning(z.shape[1]))
    fft2 = np.fft.fftshift(np.fft.fft2(z * window))
    magnitude = np.log1p(np.abs(fft2))

    dx = max(float(x_axis[-1] - x_axis[0]) / max(len(x_axis) - 1, 1), 1e-12)
    dy = max(float(y_axis[-1] - y_axis[0]) / max(len(y_axis) - 1, 1), 1e-12)
    qx = np.fft.fftshift(np.fft.fftfreq(z.shape[1], d=dx))
    qy = np.fft.fftshift(np.fft.fftfreq(z.shape[0], d=dy))

    fig = go.Figure(
        go.Heatmap(
            x=qx,
            y=qy,
            z=magnitude,
            colorscale="Magma",
            colorbar={"title": "log(1+|FFT|)"},
        )
    )
    fig.update_layout(
        template="plotly_white",
        title="2D FFT",
        height=320,
        margin={"l": 60, "r": 20, "t": 50, "b": 40},
    )
    fig.update_xaxes(
        title="qx (1/A)",
        showgrid=False,
        constrain="domain",
        scaleanchor="y",
        scaleratio=1.0,
    )
    fig.update_yaxes(
        title="qy (1/A)",
        showgrid=False,
    )
    return fig


def build_app(
    debug_mode=False,
    shutdown_on_last_client=False,
    shutdown_callback=None,
    shutdown_grace=2.0,
):
    app = dash.Dash(__name__, title="AutoBSKAN")
    browser_lifecycle = _BrowserLifecycle(
        enabled=shutdown_on_last_client,
        shutdown_callback=shutdown_callback,
        grace_seconds=shutdown_grace,
    )
    app.server.extensions["autobskan_browser_lifecycle"] = browser_lifecycle

    @app.server.before_request
    def _log_dash_update_requests():
        if debug_mode and request.path == "/_dash-update-component":
            _terminal_log(f"Dash update request: method={request.method}")

    @app.server.after_request
    def _cache_policy(response):
        # Layout and callback state must stay fresh. Dash assets already carry a
        # version or mtime in their URL and can use the browser cache safely.
        dynamic_paths = {"/", "/_dash-layout", "/_dash-dependencies"}
        if request.path in dynamic_paths or request.path == "/_dash-update-component":
            response.headers["Cache-Control"] = (
                "no-store, no-cache, must-revalidate, max-age=0"
            )
            response.headers["Pragma"] = "no-cache"
            response.headers["Expires"] = "0"
        return response

    def _browser_client_id(raw_client_id):
        value = str(raw_client_id or "")
        if not re.fullmatch(r"[A-Za-z0-9_-]{12,128}", value):
            return None
        return value

    @app.server.post("/_autobskan-client/open/<client_id>")
    def _register_browser_client(client_id):
        valid_id = _browser_client_id(client_id)
        if valid_id is None:
            return jsonify(error="Invalid browser client id."), 400
        active = browser_lifecycle.register(valid_id)
        return jsonify(active=active, shutdown_enabled=browser_lifecycle.enabled)

    @app.server.get("/_autobskan-client/events/<client_id>")
    def _stream_browser_client(client_id):
        valid_id = _browser_client_id(client_id)
        if valid_id is None:
            return jsonify(error="Invalid browser client id."), 400
        if not browser_lifecycle.enabled:
            return Response(status=204)

        connection_id = secrets.token_urlsafe(12)
        browser_lifecycle.register(valid_id, connection_id)

        @stream_with_context
        def _event_stream():
            try:
                while True:
                    yield ": keepalive\n\n"
                    time.sleep(1.0)
            except GeneratorExit:
                pass
            finally:
                browser_lifecycle.close(valid_id, connection_id)

        return Response(
            _event_stream(),
            mimetype="text/event-stream",
            headers={
                "Cache-Control": "no-store",
                "X-Accel-Buffering": "no",
            },
        )

    @app.server.post("/_autobskan-client/close/<client_id>")
    def _close_browser_client(client_id):
        valid_id = _browser_client_id(client_id)
        if valid_id is None:
            return jsonify(error="Invalid browser client id."), 400
        active = browser_lifecycle.close(valid_id)
        return jsonify(active=active, shutdown_enabled=browser_lifecycle.enabled)

    @app.server.post("/_autobskan-upload/start")
    def _start_chunk_upload():
        payload = request.get_json(silent=True) or {}
        raw_name = os.path.basename(str(payload.get("filename") or "LOCPOT"))
        safe_name = re.sub(r"[^A-Za-z0-9._-]", "_", raw_name)
        try:
            expected_size = int(payload.get("size"))
        except (TypeError, ValueError):
            expected_size = 0
        if expected_size <= 0:
            return jsonify(error="Upload size should be positive."), 400

        upload_id = secrets.token_urlsafe(24)
        stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        path = os.path.join(
            _managed_upload_dir(),
            f"volume_{stamp}_{upload_id[:8]}_{safe_name}",
        )
        with open(path, "xb"):
            pass
        with _CHUNK_UPLOAD_LOCK:
            _CHUNK_UPLOADS[upload_id] = {
                "path": path,
                "expected": expected_size,
                "received": 0,
                "filename": raw_name,
            }
        return jsonify(upload_id=upload_id, chunk_size=_CHUNK_UPLOAD_LIMIT)

    @app.server.post("/_autobskan-upload/chunk/<upload_id>")
    def _write_chunk_upload(upload_id):
        try:
            offset = int(request.headers.get("X-Upload-Offset", "-1"))
        except ValueError:
            offset = -1
        chunk = request.get_data(cache=False)
        if not chunk or len(chunk) > _CHUNK_UPLOAD_LIMIT:
            return jsonify(error="Invalid upload chunk size."), 400

        with _CHUNK_UPLOAD_LOCK:
            state = _CHUNK_UPLOADS.get(upload_id)
            if state is None:
                return jsonify(error="Upload session was not found."), 404
            if offset != state["received"]:
                return jsonify(error="Upload chunk offset is out of sequence."), 409
            if offset + len(chunk) > state["expected"]:
                return jsonify(error="Upload exceeds the declared file size."), 400
            with open(state["path"], "ab") as fileobj:
                fileobj.write(chunk)
            state["received"] += len(chunk)
            received = state["received"]
        return jsonify(received=received)

    @app.server.post("/_autobskan-upload/finish/<upload_id>")
    def _finish_chunk_upload(upload_id):
        with _CHUNK_UPLOAD_LOCK:
            state = _CHUNK_UPLOADS.get(upload_id)
            if state is None:
                return jsonify(error="Upload session was not found."), 404
            if (
                state["received"] != state["expected"]
                or os.path.getsize(state["path"]) != state["expected"]
            ):
                return jsonify(error="Upload is incomplete."), 409
            _CHUNK_UPLOADS.pop(upload_id, None)
        _terminal_log(
            "Uploaded volumetric file in chunks: "
            f"{state['filename']} -> {state['path']}"
        )
        return jsonify(path=state["path"], filename=state["filename"])

    @app.server.delete("/_autobskan-upload/<upload_id>")
    def _cancel_chunk_upload(upload_id):
        with _CHUNK_UPLOAD_LOCK:
            state = _CHUNK_UPLOADS.pop(upload_id, None)
        if state is not None:
            try:
                os.unlink(state["path"])
            except OSError:
                pass
        return jsonify(cancelled=state is not None)

    def _status_ui(show):
        return _status_wrap(bool(show) and debug_mode)

    version_info = _gui_version_info()
    app.layout = build_layout(
        defaults=DEFAULTS,
        colormaps=COLORMAP_PRESETS,
        version_label=f"v{version_info['package_version']}",
        debug_mode=debug_mode,
        empty_figure=_empty_figure,
        status_text=_status_text,
        status_wrap=_status_wrap,
    )

    app.clientside_callback(
        ClientsideFunction(
            namespace="autobskan",
            function_name="surfaceErrors",
        ),
        Output("user-error-banner", "children"),
        Output("user-error-banner", "className"),
        Input("render-status", "children"),
        Input("render-status", "color"),
        Input("volume-status", "children"),
        Input("volume-status", "color"),
        Input("structure-status", "children"),
        Input("structure-status", "color"),
    )

    @app.callback(
        Output("startup-log-store", "data"),
        Input("startup-ping", "n_intervals"),
    )
    def _startup_ping(n_intervals):
        _terminal_log(f"Client startup ping received: n_intervals={n_intervals}")
        return {
            "n_intervals": int(n_intervals or 0),
            "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
        }

    @app.callback(
        Output("structure-path", "value", allow_duplicate=True),
        Output("structure-file-store", "data", allow_duplicate=True),
        Output("structure-upload-ready-store", "data", allow_duplicate=True),
        Input("upload-structure", "contents"),
        State("upload-structure", "filename"),
        prevent_initial_call=True,
    )
    def _on_structure_upload(contents, filename):
        if not contents:
            raise PreventUpdate
        path = _save_upload(contents, filename, "structure")
        if path is None:
            _terminal_log(
                "Upload structure file failed: empty upload payload.", level="error"
            )
            raise PreventUpdate
        _terminal_log(f"Uploaded structure file: {filename or '(unnamed)'} -> {path}")
        visible_path = path if not _is_managed_upload(path) else ""
        return (
            visible_path,
            path,
            datetime.datetime.now().isoformat(timespec="microseconds"),
        )

    @app.callback(
        Output("structure-upload-filename-store", "data"),
        Input("upload-structure", "filename"),
        prevent_initial_call=True,
    )
    def _remember_structure_upload_filename(filename):
        if isinstance(filename, list):
            filename = filename[0] if filename else ""
        text = str(filename or "")
        if text:
            _terminal_log(f"Structure file selected in browser: {text}")
        return text

    @app.callback(
        Output("volume-path-feedback", "children"),
        Output("volume-path-feedback", "c"),
        Input("volume-path", "value"),
        Input("volume-file-store", "data"),
        Input("volume-upload-filename-store", "data"),
    )
    def _volume_path_feedback(volume_path, volume_store_path, selected_filename):
        typed = _safe_path(volume_path)
        stored = _safe_path(volume_store_path)

        if typed is not None:
            return f"Resolved path: {_path_with_size(typed)}", "teal"

        if isinstance(volume_path, str) and volume_path.strip():
            return f"Path not found: {volume_path}", "red"

        if stored is not None:
            display_name = os.path.basename(str(selected_filename or ""))
            if not display_name:
                display_name = os.path.basename(stored)
            return (
                f"Uploaded: {display_name} ({_human_size(os.path.getsize(stored))})",
                "teal",
            )

        if isinstance(selected_filename, str) and selected_filename.strip():
            return (
                "Selected in browser: "
                f"{selected_filename}\nWaiting for content transfer to server.",
                "yellow",
            )

        return "No volumetric file selected.", "dimmed"

    @app.callback(
        Output("structure-path-feedback", "children"),
        Output("structure-path-feedback", "c"),
        Input("structure-path", "value"),
        Input("structure-file-store", "data"),
        Input("structure-upload-filename-store", "data"),
        Input("embedded-structure-store", "data"),
        Input("structure-override-decision-store", "data"),
    )
    def _structure_path_feedback(
        structure_path,
        structure_store_path,
        selected_filename,
        embedded_structure,
        override_decision,
    ):
        typed = _safe_path(structure_path)
        stored = _safe_path(structure_store_path)
        explicit_path = typed or stored
        embedded_active = bool(
            isinstance(embedded_structure, dict)
            and embedded_structure.get("active")
        )

        if explicit_path is not None and embedded_active:
            display_name = os.path.basename(str(selected_filename or ""))
            if not display_name:
                display_name = os.path.basename(explicit_path)
            decision = (
                override_decision.get("decision")
                if isinstance(override_decision, dict)
                and override_decision.get("path") == explicit_path
                and override_decision.get("volume_path")
                == embedded_structure.get("volume_path")
                else None
            )
            if decision == "proceed":
                return (
                    f"Explicit structure active: {display_name}\n"
                    "The embedded VASP structure is overridden.",
                    "teal",
                )
            if decision == "cancel":
                return (
                    f"Selected structure ignored: {display_name}\n"
                    "The embedded VASP structure remains active.",
                    "yellow",
                )
            return (
                f"Selected structure pending confirmation: {display_name}",
                "yellow",
            )

        if typed is not None:
            return f"Resolved path: {_path_with_size(typed)}", "teal"

        if isinstance(structure_path, str) and structure_path.strip():
            return f"Path not found: {structure_path}", "red"

        if stored is not None:
            display_name = os.path.basename(str(selected_filename or ""))
            if not display_name:
                display_name = os.path.basename(stored)
            return (
                f"Uploaded: {display_name} ({_human_size(os.path.getsize(stored))})",
                "teal",
            )

        if isinstance(selected_filename, str) and selected_filename.strip():
            return (
                "Selected in browser: "
                f"{selected_filename}\nWaiting for content transfer to server.",
                "yellow",
            )

        if embedded_active:
            return (
                f"Loaded automatically from {embedded_structure['volume_kind']} "
                f"header: {embedded_structure['atoms']} atoms.",
                "teal",
            )

        return "No structure file selected.", "dimmed"

    @app.callback(
        Output("embedded-structure-message", "children"),
        Output("embedded-structure-banner", "style"),
        Input("embedded-structure-store", "data"),
        Input("structure-path", "value"),
        Input("structure-file-store", "data"),
        Input("structure-override-decision-store", "data"),
    )
    def _embedded_structure_banner(
        embedded_structure,
        structure_path,
        structure_store_path,
        override_decision,
    ):
        if not isinstance(embedded_structure, dict) or not embedded_structure.get(
            "active"
        ):
            return "", {"display": "none"}
        explicit_path = _effective_path(structure_path, structure_store_path)
        explicit_active = _structure_override_is_confirmed(
            explicit_path,
            override_decision,
            embedded_structure.get("volume_path"),
        )
        a, b, c = embedded_structure["cell"]
        state = (
            "Available; an explicitly confirmed structure is currently active."
            if explicit_active
            else "Loaded automatically and active."
        )
        message = (
            f"{embedded_structure['volume_kind']} header | "
            f"{embedded_structure['atoms']} atoms | "
            f"cell {a:.3f}, {b:.3f}, {c:.3f} A\n{state}"
        )
        return message, {"display": "grid"}

    @app.callback(
        Output("structure-override-modal", "opened"),
        Output("structure-override-pending-store", "data"),
        Input("structure-path", "value"),
        Input("structure-file-store", "data"),
        Input("embedded-structure-store", "data"),
        Input("structure-override-decision-store", "data"),
    )
    def _request_structure_override(
        structure_path,
        structure_store_path,
        embedded_structure,
        override_decision,
    ):
        explicit_path = _effective_path(structure_path, structure_store_path)
        if (
            explicit_path is None
            or not isinstance(embedded_structure, dict)
            or not embedded_structure.get("active")
        ):
            return False, {}
        if (
            isinstance(override_decision, dict)
            and override_decision.get("path") == explicit_path
            and override_decision.get("volume_path")
            == embedded_structure.get("volume_path")
            and override_decision.get("decision") in {"proceed", "cancel"}
        ):
            return False, {
                "path": explicit_path,
                "volume_path": embedded_structure.get("volume_path"),
            }
        return True, {
            "path": explicit_path,
            "name": os.path.basename(explicit_path),
            "volume_path": embedded_structure.get("volume_path"),
        }

    @app.callback(
        Output("structure-override-decision-store", "data"),
        Input("structure-override-proceed-button", "n_clicks"),
        Input("structure-override-cancel-button", "n_clicks"),
        State("structure-override-pending-store", "data"),
        prevent_initial_call=True,
    )
    def _resolve_structure_override(proceed_clicks, cancel_clicks, pending):
        path = pending.get("path") if isinstance(pending, dict) else None
        if path is None:
            raise PreventUpdate
        trigger = dash.callback_context.triggered_id
        if trigger == "structure-override-proceed-button" and proceed_clicks:
            decision = "proceed"
        elif trigger == "structure-override-cancel-button" and cancel_clicks:
            decision = "cancel"
        else:
            raise PreventUpdate
        _terminal_log(
            f"Explicit structure override decision: {decision} | file={path}"
        )
        return {
            "path": path,
            "volume_path": pending.get("volume_path"),
            "decision": decision,
        }

    @app.callback(
        Output("volume-status", "children"),
        Output("volume-status", "color"),
        Output("volume-status-wrap", "style"),
        Output("embedded-structure-store", "data"),
        Input("volume-path", "value"),
        Input("volume-file-store", "data"),
        Input("simulation-type", "value"),
        Input("input-source", "value"),
        running=[
            (
                Output("render-overlay", "className", allow_duplicate=True),
                "render-overlay is-active",
                "render-overlay",
            ),
        ],
    )
    def _volume_parse_status(volume_path, volume_store_path, simulation, input_source):
        path = _effective_path(volume_path, volume_store_path)
        if path is None:
            if isinstance(volume_path, str) and volume_path.strip():
                msg = _status_text(
                    "ERROR",
                    f"File not found: {volume_path}\nCheck path and try again.",
                )
                _terminal_log(f"Volume path not found: {volume_path}", level="warning")
                return msg, "red", _status_ui(True), {}
            return (
                _status_text("WAIT", "No volumetric file parsed yet."),
                "gray",
                _status_ui(True),
                {},
            )
        try:
            effective_source = _validate_simulation_input(
                simulation,
                input_source,
                path,
            )
            if effective_source == "BSKAN":
                current = _load_current(path)
                message = (
                    f"Parsed CURRENT: {os.path.basename(path)}\n"
                    f"Path: {path}\n"
                    f"Grid (nx, ny, nz): {tuple(int(v) for v in current.grids)}\n"
                    f"Valid isosurface: {current.iso_min:.6g} < iso < {current.iso_max:.6g}"
                )
                _terminal_log(
                    "Parsed CURRENT: "
                    f"{path} | grids={tuple(int(v) for v in current.grids)} "
                    f"| iso=({current.iso_min},{current.iso_max})"
                )
                embedded_structure = {}
            else:
                locpot = _load_locpot(path)
                kind = _vasp_volume_kind(path)
                embedded_structure = _embedded_structure_metadata(locpot, path)
                message = (
                    f"Parsed {kind}: {os.path.basename(path)}\n"
                    f"Path: {path}\n"
                    f"Grid (nx, ny, nz): {tuple(int(v) for v in locpot.grids)}\n"
                    f"Cell (a, b, c): {locpot.atoms.cell.cellpar()[0]:.3f}, "
                    f"{locpot.atoms.cell.cellpar()[1]:.3f}, {locpot.atoms.cell.cellpar()[2]:.3f} A"
                )
                _terminal_log(
                    f"Parsed {kind}: "
                    f"{path} | grids={tuple(int(v) for v in locpot.grids)}"
                )
            return (
                _status_text("OK", message),
                "green",
                _status_ui(True),
                embedded_structure,
            )
        except Exception as exc:
            message = (
                f"Failed to parse volumetric file: {path}\n"
                f"{type(exc).__name__}: {exc}\n"
                "Use a valid file for the selected style: CURRENT for bSKAN, "
                "PARCHG/LOCPOT/CHGCAR for VASP."
            )
            LOGGER.exception("Failed to parse volumetric file: %s", path)
            _terminal_log(
                f"Failed to parse volumetric file: {path} | {type(exc).__name__}: {exc}",
                level="error",
            )
            return _status_text("ERROR", message), "red", _status_ui(True), {}

    @app.callback(
        Output("structure-status", "children"),
        Output("structure-status", "color"),
        Output("structure-status-wrap", "style"),
        Input("structure-path", "value"),
        Input("structure-file-store", "data"),
        running=[
            (
                Output("render-overlay", "className", allow_duplicate=True),
                "render-overlay is-active",
                "render-overlay",
            ),
        ],
    )
    def _structure_parse_status(structure_path, structure_store_path):
        path = _effective_path(structure_path, structure_store_path)
        if path is None:
            if isinstance(structure_path, str) and structure_path.strip():
                msg = _status_text(
                    "ERROR",
                    f"File not found: {structure_path}\nCheck path and try again.",
                )
                _terminal_log(
                    f"Structure path not found: {structure_path}", level="warning"
                )
                return msg, "red", _status_ui(True)
            return (
                _status_text("WAIT", "No structure file parsed yet."),
                "gray",
                _status_ui(True),
            )
        try:
            structure = _load_structure(path)
            a, b, c, alpha, beta, gamma = structure.cell.cellpar()
            message = (
                f"Parsed structure: {os.path.basename(path)}\n"
                f"Path: {path}\n"
                f"Atoms: {len(structure)}\n"
                f"Cell: a={a:.3f}, b={b:.3f}, c={c:.3f} A, gamma={gamma:.3f} deg"
            )
            _terminal_log(
                f"Parsed structure: {path} | atoms={len(structure)} | gamma={gamma:.3f}"
            )
            return _status_text("OK", message), "green", _status_ui(True)
        except Exception as exc:
            message = (
                f"Failed to parse structure file: {path}\n"
                f"{type(exc).__name__}: {exc}\n"
                "Use a valid POSCAR/vasp file."
            )
            LOGGER.exception("Failed to parse structure file: %s", path)
            _terminal_log(
                f"Failed to parse structure file: {path} | {type(exc).__name__}: {exc}",
                level="error",
            )
            return _status_text("ERROR", message), "red", _status_ui(True)

    @app.callback(
        Output("mode-container", "style"),
        Output("fermi-level-container", "style"),
        Output("fit-radius-container", "style"),
        Output("simulation-definition", "children"),
        Output("volume-source-hint", "children"),
        Output("volume-file-label", "children"),
        Output("stm-mode", "value"),
        Output("colormap-dropdown", "value"),
        Input("simulation-type", "value"),
        Input("input-source", "value"),
        State("colormap-dropdown", "value"),
    )
    def _sync_simulation_mode(simulation, input_source, current_cmap):
        if simulation == "PHI_APP":
            file_label = "CURRENT file" if input_source == "BSKAN" else "PARCHG file"
            source_hint = (
                "bSKAN CURRENT: current/LDOS decay"
                if input_source == "BSKAN"
                else "VASP PARCHG: partial-density decay"
            )
            cmap_value = (
                DEFAULTS["colormap_lwf"]
                if current_cmap in [None, DEFAULTS["colormap_stm"]]
                else no_update
            )
            return (
                {"display": "none"},
                {"display": "none"},
                {"display": "grid"},
                r"**Apparent barrier height** $\Phi_{\rm app}=0.952495(\partial_z\ln Q)^2$; constant-height fit.",
                source_hint,
                file_label,
                "CONSTANT HEIGHT",
                cmap_value,
            )

        if simulation == "LWF":
            cmap_value = (
                DEFAULTS["colormap_lwf"]
                if current_cmap in [None, DEFAULTS["colormap_stm"]]
                else no_update
            )
            return (
                {"display": "none"},
                {"display": "grid"},
                {"display": "none"},
                r"**Local work function** $\Phi_{\rm loc}(x,y;z_0)=V_{\rm LOCPOT}(x,y,z_0)-E_{\rm F}$; constant height.",
                "VASP LOCPOT only; E_F is read from OUTCAR or entered below",
                "LOCPOT file",
                "CONSTANT HEIGHT",
                cmap_value,
            )

        cmap_value = (
            DEFAULTS["colormap_stm"]
            if current_cmap in [None, DEFAULTS["colormap_lwf"]]
            else no_update
        )
        file_label = (
            "CURRENT file" if input_source == "BSKAN" else "VASP volumetric file"
        )
        source_hint = (
            "bSKAN CURRENT"
            if input_source == "BSKAN"
            else "VASP PARCHG, LOCPOT, or CHGCAR"
        )
        return (
            {"display": "block"},
            {"display": "none"},
            {"display": "none"},
            "**STM simulation** surface analysis.",
            source_hint,
            file_label,
            no_update,
            cmap_value,
        )

    @app.callback(
        Output("fermi-level-help", "children"),
        Output("fermi-level-help", "c"),
        Input("simulation-type", "value"),
        Input("volume-path", "value"),
        Input("volume-file-store", "data"),
        Input("fermi-level-input", "value"),
    )
    def _fermi_level_feedback(
        simulation,
        volume_path,
        volume_store_path,
        manual_value,
    ):
        if simulation != "LWF":
            return "", "dimmed"
        if manual_value is not None:
            try:
                value, _source = _resolve_lwf_fermi_level("", manual_value)
                return f"Using manual E_F = {value:.8g} eV.", "teal"
            except (TypeError, ValueError) as exc:
                return f"Invalid E_F: {exc}", "red"

        volume_real = _effective_path(volume_path, volume_store_path)
        if volume_real is None:
            return "Enter E_F or select a LOCPOT with a sibling OUTCAR.", "dimmed"
        outcar_path = lwfplot.sibling_outcar_path(volume_real)
        if outcar_path is None:
            return "No sibling OUTCAR found. Enter E_F manually.", "yellow"
        try:
            value, _source = _resolve_lwf_fermi_level(volume_real)
            return f"Auto-detected E_F = {value:.8g} eV from {outcar_path}.", "teal"
        except (OSError, ValueError) as exc:
            return f"Could not read E_F from OUTCAR: {exc}", "red"

    @app.callback(
        Output("gamma-input", "value"),
        Output("gamma-input", "disabled"),
        Input("structure-path", "value"),
        Input("structure-file-store", "data"),
        Input("volume-path", "value"),
        Input("volume-file-store", "data"),
        Input("input-source", "value"),
        Input("structure-override-decision-store", "data"),
        State("gamma-input", "value"),
    )
    def _sync_gamma_from_structure(
        structure_path,
        structure_store_path,
        volume_path,
        volume_store_path,
        input_source,
        override_decision,
        gamma_value,
    ):
        structure_real = _effective_path(structure_path, structure_store_path)
        volume_real = _effective_path(volume_path, volume_store_path)
        try:
            if input_source == "VASP" and volume_real is not None:
                if _structure_override_is_confirmed(
                    structure_real,
                    override_decision,
                    volume_real,
                ):
                    structure = _load_structure(structure_real)
                else:
                    structure = _embedded_structure_from_volume(
                        _load_locpot(volume_real)
                    )
                return float(np.round(structure.cell.cellpar()[-1], 6)), True
            if structure_real is not None:
                structure = _load_structure(structure_real)
                return float(np.round(structure.cell.cellpar()[-1], 6)), True
        except Exception:
            pass
        return gamma_value, False

    @app.callback(
        Output("isosurface-label", "children"),
        Output("isosurface-slider", "min"),
        Output("isosurface-slider", "max"),
        Output("isosurface-slider", "step"),
        Output("isosurface-range-default-store", "data"),
        Output("isosurface-help", "children"),
        Input("simulation-type", "value"),
        Input("input-source", "value"),
        Input("stm-mode", "value"),
        Input("volume-path", "value"),
        Input("structure-path", "value"),
        Input("volume-file-store", "data"),
        Input("structure-file-store", "data"),
        Input("structure-override-decision-store", "data"),
        Input("fit-radius-slider", "value"),
    )
    def _update_isosurface_range(
        simulation,
        input_source,
        stm_mode,
        volume_path,
        structure_path,
        volume_store_path,
        structure_store_path,
        override_decision,
        fit_radius_value,
    ):
        try:
            fit_radius = float(fit_radius_value)
        except (TypeError, ValueError):
            fit_radius = DEFAULTS["fit_radius"]
        if not np.isfinite(fit_radius) or fit_radius <= 0.0:
            fit_radius = DEFAULTS["fit_radius"]

        def range_value(default_value, minimum, maximum):
            return float(np.clip(default_value, minimum, maximum))

        volume_real = _effective_path(volume_path, volume_store_path)
        if volume_real is None:
            if simulation == "STM" and stm_mode == "CONSTANT CURRENT":
                return (
                    "log10(Isosurface)",
                    -5.0,
                    0.0,
                    0.01,
                    -2.5,
                    "Range set after file load.",
                )
            label = (
                "Fit height (A)"
                if simulation == "PHI_APP"
                else "Height above surface (A)"
            )
            if simulation == "PHI_APP":
                h_min = fit_radius
                h_max = max(6.0 - fit_radius, h_min + 1e-3)
                value = range_value(2.5, h_min, h_max)
                return (
                    label,
                    h_min,
                    h_max,
                    0.001,
                    value,
                    f"Range reserves the full +/-{fit_radius:.3f} A fit window.",
                )
            return label, 0.0, 6.0, 0.001, 2.5, "Range set after file load."

        try:
            effective_source = _validate_simulation_input(
                simulation,
                input_source,
                volume_real,
            )
            if (
                effective_source == "BSKAN"
                and simulation == "STM"
                and stm_mode == "CONSTANT CURRENT"
            ):
                current = _load_current(volume_real)
                iso_min = float(np.ceil(np.log10(current.iso_min)))
                iso_max = float(np.floor(np.log10(current.iso_max)))
                if iso_max <= iso_min:
                    center = float(np.log10((current.iso_min + current.iso_max) * 0.5))
                    iso_min = center - 0.5
                    iso_max = center + 0.5
                value = float(np.clip(-2.5, iso_min, iso_max))
                msg = f"Valid range: {iso_min:.2f} to {iso_max:.2f} (log10 scale)"
                return "log10(Isosurface)", iso_min, iso_max, 0.01, value, msg

            if effective_source == "BSKAN":
                current = _load_current(volume_real)
                h_min, h_max = lwfplot.current_height_window(current)
                if simulation == "PHI_APP":
                    h_min += fit_radius
                    h_max -= fit_radius
                    if h_max <= h_min:
                        raise ValueError(
                            "Selected fit half-width does not fit inside the CURRENT z-window."
                        )
                value = range_value(min(2.5, h_max), h_min, h_max)
                msg = (
                    f"Valid range: {h_min:.3f} to {h_max:.3f} A "
                    f"(from CURRENT c={float(current.cell[2, 2]):.3f} A, no extrapolation)"
                )
                if simulation == "PHI_APP":
                    msg += f"; reserves +/-{fit_radius:.3f} A for fitting"
                label = (
                    "Fit height (A)"
                    if simulation == "PHI_APP"
                    else "Height (A)"
                )
                return label, h_min, h_max, 0.001, value, msg

            locpot = _load_locpot(volume_real)
            c_len = _lwf_cell_c_length(locpot)
            structure_real = _effective_path(structure_path, structure_store_path)
            structure = _embedded_structure_from_volume(locpot)
            if _structure_override_is_confirmed(
                structure_real,
                override_decision,
                volume_real,
            ):
                try:
                    structure = _load_structure(structure_real)
                except Exception:
                    structure = _embedded_structure_from_volume(locpot)
            top = _lwf_topmost_z(structure, locpot)

            if simulation == "STM" and stm_mode == "CONSTANT CURRENT":
                iso_min_raw, iso_max_raw = stmplot.get_vasp_density_isorange(locpot)
                iso_min = float(np.log10(max(iso_min_raw, np.nextafter(0.0, 1.0))))
                iso_max = float(np.log10(iso_max_raw))
                if iso_max <= iso_min:
                    center = float(np.log10((iso_min_raw + iso_max_raw) * 0.5))
                    iso_min = center - 0.5
                    iso_max = center + 0.5
                value = float(np.clip((iso_min + iso_max) * 0.5, iso_min, iso_max))
                msg = f"Valid VASP density range: {iso_min:.2f} to {iso_max:.2f} (log10 scale)"
                return "log10(Isosurface)", iso_min, iso_max, 0.01, value, msg

            fit_margin = fit_radius if simulation == "PHI_APP" else 0.0
            h_min = fit_margin
            z_ceiling = (
                float(lwfplot.vasp_z_axis(locpot)[-1])
                if simulation == "PHI_APP"
                else c_len
            )
            h_max = max(z_ceiling - top - fit_margin, 1e-3)
            if h_max <= h_min:
                raise ValueError(
                    "Selected fit half-width does not fit inside the vacuum z-window."
                )
            value = range_value(min(2.5, h_max), h_min, h_max)
            msg = f"Valid range: {h_min:.3f} to {h_max:.3f} A"
            if simulation == "PHI_APP":
                msg += f"; reserves +/-{fit_radius:.3f} A for fitting"
            label = (
                "Fit height (A)"
                if simulation == "PHI_APP"
                else "Height above surface (A)"
            )
            return label, h_min, h_max, 0.001, value, msg
        except Exception as exc:
            label = (
                "Fit height (A)"
                if simulation == "PHI_APP"
                else "Height above surface (A)"
            )
            return label, 0.0, 6.0, 0.001, 2.5, f"Range fallback: {exc}"

    @app.callback(
        Output("volume-path", "value", allow_duplicate=True),
        Output("structure-path", "value", allow_duplicate=True),
        Output("volume-file-store", "data", allow_duplicate=True),
        Output("structure-file-store", "data", allow_duplicate=True),
        Output("volume-upload-filename-store", "data", allow_duplicate=True),
        Output("structure-upload-filename-store", "data", allow_duplicate=True),
        Output("simulation-type", "value"),
        Output("input-source", "value"),
        Output("stm-mode", "value", allow_duplicate=True),
        Output("fermi-level-input", "value"),
        Output("fit-radius-slider", "value"),
        Output("colormap-dropdown", "value", allow_duplicate=True),
        Output("brightness-slider", "value"),
        Output("color-range-mode", "value"),
        Output("contrast-slider", "value"),
        Output("layers-input", "value"),
        Output("view-options", "value"),
        Output("repeat-x", "value"),
        Output("repeat-y", "value"),
        Output("scalebar-length", "value"),
        Output("atom-radius-type", "value"),
        Output("atom-radius-scale", "value"),
        Output("gaussian-sigma", "value"),
        Output("xy-upsample", "value"),
        Output("display-quality", "value"),
        Output("grid-line-width", "value"),
        Output("grid-line-color", "value"),
        Output("export-format", "value"),
        Output("export-width", "value", allow_duplicate=True),
        Output("export-height", "value", allow_duplicate=True),
        Output("colorbar-width", "value"),
        Output("colorbar-height", "value"),
        Output("line-profile-enabled", "checked"),
        Output("line-magnet-enabled", "checked"),
        Output("fft-enabled", "checked"),
        Output("line-points-store", "data", allow_duplicate=True),
        Output("render-status", "children", allow_duplicate=True),
        Output("render-status", "color", allow_duplicate=True),
        Output("render-status-wrap", "style", allow_duplicate=True),
        Input("set-default-button", "n_clicks"),
        State("volume-path", "value"),
        State("structure-path", "value"),
        State("volume-file-store", "data"),
        State("structure-file-store", "data"),
        State("volume-upload-filename-store", "data"),
        State("structure-upload-filename-store", "data"),
        prevent_initial_call=True,
    )
    def _set_defaults(
        _n_clicks,
        volume_path_current,
        structure_path_current,
        volume_store_current,
        structure_store_current,
        volume_upload_filename_current,
        structure_upload_filename_current,
    ):
        _terminal_log("Reset to defaults triggered. Restoring built-in GUI defaults.")
        reset_msg = _status_text(
            "OK",
            "Defaults restored.\nPlot options were reset. Selected file paths were kept.",
        )
        return (
            volume_path_current,
            structure_path_current,
            volume_store_current,
            structure_store_current,
            volume_upload_filename_current,
            structure_upload_filename_current,
            DEFAULTS["simulation"],
            DEFAULTS["input_source"],
            DEFAULTS["mode"],
            DEFAULTS["fermi_level"],
            DEFAULTS["fit_radius"],
            DEFAULTS["colormap_stm"],
            DEFAULTS["brightness"],
            DEFAULTS["color_range_mode"],
            DEFAULTS["contrast"],
            DEFAULTS["layers"],
            DEFAULTS["view_options"],
            DEFAULTS["repeat_x"],
            DEFAULTS["repeat_y"],
            DEFAULTS["scalebar"],
            DEFAULTS["atom_radius_type"],
            DEFAULTS["atom_radius_scale"],
            DEFAULTS["gaussian_sigma"],
            DEFAULTS["xy_upsample"],
            DEFAULTS["display_quality"],
            DEFAULTS["grid_line_width"],
            DEFAULTS["grid_line_color"],
            DEFAULTS["export_format"],
            DEFAULTS["export_width"],
            DEFAULTS["export_height"],
            DEFAULTS["colorbar_width"],
            DEFAULTS["colorbar_height"],
            True,
            DEFAULTS["line_magnet_enabled"],
            True,
            [],
            reset_msg,
            "blue",
            _status_ui(True),
        )

    @app.callback(
        Output("appearance-store", "data"),
        Input("colormap-dropdown", "value"),
        Input("contrast-slider", "value"),
    )
    def _update_appearance_store(cmap_name, contrast):
        return {
            "cmap_name": cmap_name,
            "contrast": float(contrast or 0.0),
            "colorscale": _colorscale(cmap_name, float(contrast or 0.0)),
        }

    @app.callback(
        Output("color-range-mode", "value", allow_duplicate=True),
        Input("color-range-auto-button", "n_clicks"),
        prevent_initial_call=True,
    )
    def _restore_auto_color_range(n_clicks):
        if not n_clicks:
            raise PreventUpdate
        return "AUTO"

    app.clientside_callback(
        ClientsideFunction(
            namespace="autobskan",
            function_name="syncColorRangeControls",
        ),
        Output("color-vmin-input", "value"),
        Output("color-vmax-input", "value"),
        Output("color-vmin-input", "disabled"),
        Output("color-vmax-input", "disabled"),
        Output("brightness-slider", "disabled"),
        Output("brightness-input", "disabled"),
        Input("color-range-mode", "value"),
        Input("plot-meta-store", "data"),
        State("color-vmin-input", "value"),
        State("color-vmax-input", "value"),
    )

    app.clientside_callback(
        ClientsideFunction(
            namespace="autobskan",
            function_name="colorRangeStatus",
        ),
        Output("color-range-help", "children"),
        Output("color-range-help", "className"),
        Input("color-range-mode", "value"),
        Input("color-vmin-input", "value"),
        Input("color-vmax-input", "value"),
        Input("plot-meta-store", "data"),
    )

    app.clientside_callback(
        ClientsideFunction(
            namespace="autobskan",
            function_name="syncOverlayControls",
        ),
        Output("repeat-x", "disabled"),
        Output("repeat-y", "disabled"),
        Output("layers-input", "disabled"),
        Output("atom-radius-type", "disabled"),
        Output("atom-radius-scale", "disabled"),
        Output("atom-radius-scale-input", "disabled"),
        Output("gaussian-sigma", "disabled"),
        Output("gaussian-sigma-input", "disabled"),
        Output("grid-line-width", "disabled"),
        Output("grid-line-color", "disabled"),
        Input("view-options", "value"),
    )

    for slider_id, input_id in [
        ("isosurface-slider", "isosurface-input"),
        ("fit-radius-slider", "fit-radius-input"),
        ("brightness-slider", "brightness-input"),
        ("contrast-slider", "contrast-input"),
        ("gaussian-sigma", "gaussian-sigma-input"),
        ("atom-radius-scale", "atom-radius-scale-input"),
    ]:
        app.clientside_callback(
            ClientsideFunction(
                namespace="autobskan",
                function_name="syncSliderNumber",
            ),
            Output(slider_id, "value", allow_duplicate=True),
            Output(input_id, "value"),
            Input(slider_id, "value"),
            Input(input_id, "value"),
            State(slider_id, "min"),
            State(slider_id, "max"),
            prevent_initial_call=True,
        )

    app.clientside_callback(
        ClientsideFunction(
            namespace="autobskan",
            function_name="syncSliderBounds",
        ),
        Output("isosurface-input", "min"),
        Output("isosurface-input", "max"),
        Output("isosurface-input", "step"),
        Input("isosurface-slider", "min"),
        Input("isosurface-slider", "max"),
        Input("isosurface-slider", "step"),
    )

    app.clientside_callback(
        ClientsideFunction(
            namespace="autobskan",
            function_name="applySliderRange",
        ),
        Output("isosurface-slider", "value", allow_duplicate=True),
        Input("isosurface-slider", "min"),
        Input("isosurface-slider", "max"),
        Input("isosurface-range-default-store", "data"),
        State("isosurface-slider", "value"),
        prevent_initial_call=True,
    )

    app.clientside_callback(
        ClientsideFunction(
            namespace="autobskan",
            function_name="scheduleRenderControls",
        ),
        Output("isosurface-render-store", "data"),
        Input("isosurface-slider", "value"),
        Input("gaussian-sigma", "value"),
        Input("atom-radius-scale", "value"),
        Input("fit-radius-slider", "value"),
        State("view-options", "value"),
        State("simulation-type", "value"),
        prevent_initial_call=True,
    )

    app.clientside_callback(
        ClientsideFunction(
            namespace="autobskan",
            function_name="scheduleContextRender",
        ),
        Output("render-request-store", "data"),
        Input("gamma-input", "value"),
        Input("layers-input", "value"),
        Input("repeat-x", "value"),
        Input("repeat-y", "value"),
        Input("atom-radius-type", "value"),
        Input("grid-line-width", "value"),
        Input("grid-line-color", "value"),
        State("view-options", "value"),
        State("input-source", "value"),
        prevent_initial_call=True,
    )

    app.clientside_callback(
        ClientsideFunction(
            namespace="autobskan",
            function_name="beginRender",
        ),
        Output("render-overlay", "className", allow_duplicate=True),
        Input("render-button", "n_clicks"),
        Input("volume-path", "value"),
        Input("volume-upload-ready-store", "data"),
        Input("structure-path", "value"),
        Input("structure-upload-ready-store", "data"),
        Input("structure-override-decision-store", "data"),
        Input("simulation-type", "value"),
        Input("input-source", "value"),
        Input("stm-mode", "value"),
        Input("fermi-level-input", "value"),
        Input("isosurface-render-store", "data"),
        Input("view-options", "value"),
        Input("scalebar-length", "value"),
        Input("xy-upsample", "value"),
        Input("display-quality", "value"),
        Input("render-request-store", "data"),
        prevent_initial_call=True,
    )

    app.clientside_callback(
        ClientsideFunction(
            namespace="autobskan",
            function_name="updateAppearance",
        ),
        Output("main-image", "figure", allow_duplicate=True),
        Output("colorbar-graph", "figure", allow_duplicate=True),
        Output("line-profile-graph", "figure", allow_duplicate=True),
        Input("appearance-store", "data"),
        Input("brightness-slider", "value"),
        Input("color-range-mode", "value"),
        Input("color-vmin-input", "value"),
        Input("color-vmax-input", "value"),
        State("main-image", "figure"),
        State("colorbar-graph", "figure"),
        State("line-profile-graph", "figure"),
        State("plot-meta-store", "data"),
        prevent_initial_call=True,
    )

    app.clientside_callback(
        ClientsideFunction(
            namespace="autobskan",
            function_name="syncExportDimensions",
        ),
        Output("export-width", "value"),
        Output("export-height", "value"),
        Output("export-size-status", "children"),
        Input("export-width", "value"),
        Input("export-height", "value"),
        Input("plot-meta-store", "data"),
    )

    app.clientside_callback(
        ClientsideFunction(
            namespace="autobskan",
            function_name="toggleControls",
        ),
        Output("controls-open-store", "data"),
        Output("control-rail", "className"),
        Output("controls-backdrop", "className"),
        Output("app-interactive-wrap", "className"),
        Input("controls-toggle-button", "n_clicks"),
        Input("controls-close-button", "n_clicks"),
        Input("controls-backdrop", "n_clicks"),
        Input("render-button", "n_clicks"),
        Input("controls-keyboard-store", "data"),
        State("controls-open-store", "data"),
        prevent_initial_call=True,
    )

    app.clientside_callback(
        ClientsideFunction(
            namespace="autobskan",
            function_name="updateLineSelection",
        ),
        Output("line-points-store", "data"),
        Output("main-image", "figure", allow_duplicate=True),
        Input("main-image", "clickData"),
        Input("clear-line-button", "n_clicks"),
        Input("line-profile-enabled", "checked"),
        State("line-points-store", "data"),
        State("main-image", "figure"),
        prevent_initial_call=True,
    )

    app.clientside_callback(
        ClientsideFunction(
            namespace="autobskan",
            function_name="applyLineDisplay",
        ),
        Output("main-image", "figure", allow_duplicate=True),
        Input("line-display-points-store", "data"),
        State("main-image", "figure"),
        prevent_initial_call=True,
    )

    @app.callback(
        Output("main-image", "figure"),
        Output("colorbar-graph", "figure"),
        Output("line-profile-graph", "figure"),
        Output("fft-graph", "figure"),
        Output("plot-meta-store", "data"),
        Output("surface-key-store", "data"),
        Output("render-status", "children"),
        Output("render-status", "color"),
        Output("render-status-wrap", "style"),
        Input("render-button", "n_clicks"),
        Input("volume-path", "value"),
        Input("volume-upload-ready-store", "data"),
        Input("structure-path", "value"),
        Input("structure-upload-ready-store", "data"),
        Input("structure-override-decision-store", "data"),
        Input("simulation-type", "value"),
        Input("input-source", "value"),
        Input("stm-mode", "value"),
        Input("fermi-level-input", "value"),
        Input("isosurface-render-store", "data"),
        Input("view-options", "value"),
        Input("scalebar-length", "value"),
        Input("xy-upsample", "value"),
        Input("display-quality", "value"),
        Input("render-request-store", "data"),
        State("volume-file-store", "data"),
        State("structure-file-store", "data"),
        State("line-profile-enabled", "checked"),
        State("line-magnet-enabled", "checked"),
        State("fft-enabled", "checked"),
        State("line-points-store", "data"),
        State("colormap-dropdown", "value"),
        State("brightness-slider", "value"),
        State("contrast-slider", "value"),
        State("isosurface-slider", "value"),
        State("atom-radius-scale", "value"),
        State("gaussian-sigma", "value"),
        State("fit-radius-slider", "value"),
        State("gamma-input", "value"),
        State("layers-input", "value"),
        State("repeat-x", "value"),
        State("repeat-y", "value"),
        State("atom-radius-type", "value"),
        State("grid-line-width", "value"),
        State("grid-line-color", "value"),
        State("color-range-mode", "value"),
        State("color-vmin-input", "value"),
        State("color-vmax-input", "value"),
        running=[
            (
                Output("processing-indicator", "style"),
                {
                    "display": "inline-flex",
                    "alignItems": "center",
                    "gap": "0.35rem",
                },
                {
                    "display": "none",
                },
            ),
            (
                Output("render-overlay", "className"),
                "render-overlay is-active",
                "render-overlay",
            ),
        ],
    )
    def _render(
        _render_clicks,
        volume_path,
        _volume_upload_ready,
        structure_path,
        _structure_upload_ready,
        structure_override_decision,
        simulation,
        input_source,
        stm_mode,
        fermi_level_input,
        iso_render_value,
        view_options,
        scalebar_length,
        xy_upsample,
        display_quality,
        _render_request,
        volume_store_path,
        structure_store_path,
        line_enabled,
        line_magnet_enabled,
        fft_enabled,
        line_points,
        cmap_name,
        brightness,
        contrast,
        iso_slider_current,
        atom_radius_scale_current,
        gaussian_sigma_current,
        fit_radius_current,
        gamma_value,
        layers,
        repeat_x,
        repeat_y,
        atom_radius_type,
        grid_line_width,
        grid_line_color,
        color_range_mode,
        color_vmin,
        color_vmax,
    ):
        trigger = dash.callback_context.triggered_id
        render_controls = iso_render_value if isinstance(iso_render_value, dict) else {}
        if trigger == "isosurface-render-store":
            iso_slider = render_controls.get("isosurface", iso_slider_current)
            atom_radius_scale = render_controls.get(
                "atom_radius_scale", atom_radius_scale_current
            )
            gaussian_sigma = render_controls.get(
                "gaussian_sigma", gaussian_sigma_current
            )
            fit_radius = render_controls.get("fit_radius", fit_radius_current)
        else:
            iso_slider = iso_slider_current
            atom_radius_scale = atom_radius_scale_current
            gaussian_sigma = gaussian_sigma_current
            fit_radius = fit_radius_current
        try:
            fit_radius = float(fit_radius)
        except (TypeError, ValueError):
            fit_radius = float("nan")
        _terminal_log(
            "Render callback: "
            f"simulation={simulation}, source={input_source}, mode={stm_mode}, "
            f"volume={volume_path}, volume_store={volume_store_path}, "
            f"structure={structure_path}, structure_store={structure_store_path}"
        )
        volume_real = _effective_path(volume_path, volume_store_path)
        if volume_real is None:
            empty_main = _empty_figure(
                "Surface map",
                "Set CURRENT/PARCHG/LOCPOT file path or upload a file.",
                height=580,
            )
            empty_cbar = _empty_figure(
                "Color scale", "Render an image first.", height=145
            )
            empty_line = _empty_figure(
                "Line profile", "Select two points after rendering.", height=260
            )
            empty_fft = _empty_figure("2D FFT", "Render an image first.", height=260)
            if isinstance(volume_path, str) and volume_path.strip():
                info = _status_text(
                    "ERROR",
                    f"Cannot access file path: {volume_path}\n"
                    "Check path and permissions, then retry.",
                )
                _terminal_log(
                    f"Render skipped, invalid volumetric path: {volume_path}",
                    level="warning",
                )
                return (
                    empty_main,
                    empty_cbar,
                    empty_line,
                    empty_fft,
                    {},
                    "",
                    info,
                    "red",
                    _status_ui(True),
                )

            info = _status_text(
                "WAIT",
                "No volumetric file selected.\nProvide CURRENT/PARCHG/LOCPOT path or upload a file.",
            )
            _terminal_log("Render waiting for volumetric file.", level="warning")
            return (
                empty_main,
                empty_cbar,
                empty_line,
                empty_fft,
                {},
                "",
                info,
                "gray",
                _status_ui(True),
            )

        explicit_structure = None
        structure_real = _effective_path(structure_path, structure_store_path)
        if structure_real is not None:
            try:
                explicit_structure = _load_structure(structure_real)
            except Exception:
                explicit_structure = None
        structure = explicit_structure
        structure_source = (
            f"explicit file {os.path.basename(structure_real)}"
            if explicit_structure is not None
            else None
        )
        if input_source == "VASP" and structure_real is not None:
            override_confirmed = _structure_override_is_confirmed(
                structure_real,
                structure_override_decision,
                volume_real,
            )
            if trigger in {
                "structure-path",
                "structure-upload-ready-store",
            } and not override_confirmed:
                raise PreventUpdate
            if (
                trigger == "structure-override-decision-store"
                and not override_confirmed
            ):
                raise PreventUpdate

        mode = (
            "CONSTANT HEIGHT"
            if simulation in {"PHI_APP", "LWF"}
            else stm_mode
        )
        nx, ny = _parse_repeat(repeat_x, repeat_y)
        opts = set(view_options or [])
        show_repeat = "show_repeated" in opts
        show_atoms = "show_atoms" in opts
        show_blur = "show_blurred" in opts
        show_grid = "show_grids" in opts
        stm_height_window = None
        fermi_level_used = None
        fermi_level_source = None

        try:
            if simulation == "PHI_APP" and (
                not np.isfinite(fit_radius) or fit_radius <= 0.0
            ):
                raise ValueError("Fit half-width should be a positive number in A.")
            effective_source = _validate_simulation_input(
                simulation,
                input_source,
                volume_real,
            )
            if effective_source == "BSKAN":
                current = _load_current(volume_real)
                stm_height_window = lwfplot.current_height_window(current)
                if simulation == "STM":
                    if mode == "CONSTANT CURRENT":
                        iso_candidate = 10.0 ** float(iso_slider)
                        iso_lower = np.nextafter(float(current.iso_min), np.inf)
                        iso_upper = np.nextafter(float(current.iso_max), -np.inf)
                        if iso_lower >= iso_upper:
                            iso_lower = float(current.iso_min)
                            iso_upper = float(current.iso_max)
                        iso_value = float(np.clip(iso_candidate, iso_lower, iso_upper))
                    else:
                        h_min = np.nextafter(float(stm_height_window[0]), np.inf)
                        h_max = np.nextafter(float(stm_height_window[1]), -np.inf)
                        if h_min >= h_max:
                            h_min, h_max = stm_height_window
                        iso_value = float(np.clip(float(iso_slider), h_min, h_max))
                    z_raw = _surface_stm(current, mode, iso_value)
                    value_unit = "A" if mode == "CONSTANT CURRENT" else "nA"
                    cbar_title = (
                        "Height (A)" if mode == "CONSTANT CURRENT" else "Current (nA)"
                    )
                elif simulation == "PHI_APP":
                    h_min = np.nextafter(float(stm_height_window[0]), np.inf)
                    h_max = np.nextafter(float(stm_height_window[1]), -np.inf)
                    if h_min >= h_max:
                        h_min, h_max = stm_height_window
                    iso_value = float(np.clip(float(iso_slider), h_min, h_max))
                    if (
                        iso_value - fit_radius < float(stm_height_window[0])
                        or iso_value + fit_radius > float(stm_height_window[1])
                    ):
                        raise ValueError(
                            "Fit height +/- half-width must stay inside the "
                            "CURRENT z-window. Reduce the fit half-width or "
                            "move the fit height."
                        )
                    z_raw = lwfplot.apparent_barrier_from_current(
                        current,
                        height=iso_value,
                        fit_radius=fit_radius,
                    )
                    value_unit = "eV"
                    cbar_title = r"$\Phi_{\rm app}\; (\mathrm{eV})$"
                else:
                    raise ValueError("Phi_loc requires VASP LOCPOT data.")
                X, Y, Z = _prepare_stm_grid(
                    current=current,
                    structure=structure,
                    z_raw=z_raw,
                    show_repeat=show_repeat,
                    nx=nx,
                    ny=ny,
                    gamma=float(gamma_value if gamma_value is not None else 90.0),
                )
            else:
                locpot = _load_locpot(volume_real)
                structure = _embedded_structure_from_volume(locpot)
                structure_source = f"embedded {_vasp_volume_kind(volume_real)} header"
                if (
                    explicit_structure is not None
                    and _structure_override_is_confirmed(
                        structure_real,
                        structure_override_decision,
                        volume_real,
                    )
                ):
                    structure = explicit_structure
                    structure_source = (
                        f"explicit file {os.path.basename(structure_real)}"
                    )
                top = _lwf_topmost_z(structure, locpot)
                if simulation == "STM":
                    if (
                        mode == "CONSTANT CURRENT"
                        and _vasp_volume_kind(volume_real) == "LOCPOT"
                    ):
                        mode = "CONSTANT HEIGHT"
                    if mode == "CONSTANT CURRENT":
                        iso_min_raw, iso_max_raw = stmplot.get_vasp_density_isorange(
                            locpot
                        )
                        iso_value = float(
                            np.clip(
                                10.0 ** float(iso_slider),
                                np.nextafter(float(iso_min_raw), np.inf),
                                np.nextafter(float(iso_max_raw), -np.inf),
                            )
                        )
                    else:
                        c_len = _lwf_cell_c_length(locpot)
                        h_min = 0.0
                        h_max = max(c_len - top, 1e-3)
                        iso_value = float(np.clip(float(iso_slider), h_min, h_max))
                    z_raw = stmplot.get_surface_vasp(
                        locpot,
                        iso_value,
                        constant_current=mode == "CONSTANT CURRENT",
                        interpolate=(
                            "linear"
                            if mode == "CONSTANT HEIGHT"
                            and _vasp_volume_kind(volume_real) == "LOCPOT"
                            else "exponential"
                        ),
                    )
                    if mode == "CONSTANT CURRENT":
                        cbar_title = "Height (A)"
                        value_unit = "A"
                    else:
                        cbar_title, value_unit = _vasp_constant_height_label(
                            volume_real
                        )
                elif simulation == "PHI_APP":
                    fit_height = float(iso_slider)
                    fit_center_abs = top + fit_height
                    z_ceiling = float(lwfplot.vasp_z_axis(locpot)[-1])
                    if (
                        fit_center_abs - fit_radius < top
                        or fit_center_abs + fit_radius > z_ceiling
                    ):
                        raise ValueError(
                            "Fit height +/- half-width must stay inside the "
                            "PARCHG vacuum z-window. Reduce the fit half-width "
                            "or move the fit height."
                        )
                    z_raw = _surface_apparent_barrier_from_parchg(
                        locpot,
                        fit_height,
                        topmost=top,
                        fit_radius=fit_radius,
                    )
                    cbar_title = r"$\Phi_{\rm app}\; (\mathrm{eV})$"
                    value_unit = "eV"
                else:
                    fermi_level_used, fermi_level_source = _resolve_lwf_fermi_level(
                        volume_real,
                        fermi_level_input,
                    )
                    z_raw = _surface_lwf(
                        locpot,
                        float(iso_slider),
                        fermi_level=fermi_level_used,
                        topmost=top,
                    )
                    cbar_title = r"$\Phi_{\rm loc}\; (\mathrm{eV})$"
                    value_unit = "eV"
                X, Y, Z = _prepare_lwf_grid(
                    locpot=locpot,
                    structure=structure,
                    z_raw=z_raw,
                    show_repeat=show_repeat,
                    nx=nx,
                    ny=ny,
                )
        except Exception as exc:
            empty_main = _empty_figure(
                "Surface map", f"Rendering failed: {exc}", height=580
            )
            empty_cbar = _empty_figure(
                "Color scale", "Fix rendering first.", height=145
            )
            empty_line = _empty_figure(
                "Line profile", "Fix rendering first.", height=260
            )
            empty_fft = _empty_figure("2D FFT", "Fix rendering first.", height=260)
            message = _status_text(
                "ERROR",
                f"Rendering failed for {os.path.basename(volume_real)}\n"
                f"{type(exc).__name__}: {exc}",
            )
            LOGGER.exception("Rendering failed for %s", volume_real)
            _terminal_log(
                f"Rendering failed for {volume_real} | {type(exc).__name__}: {exc}",
                level="error",
            )
            return (
                empty_main,
                empty_cbar,
                empty_line,
                empty_fft,
                {},
                "",
                message,
                "red",
                _status_ui(True),
            )

        X, Y, Z, upsample_factor = _upsample_surface_xy(X, Y, Z, xy_upsample)
        full_x_axis, full_y_axis, Z_full = _rectilinear_surface(X, Y, Z)
        blur_applied = bool(show_blur and float(gaussian_sigma) > 0)
        Z_visual_full = (
            _gaussian_blur_surface(Z_full, gaussian_sigma)
            if blur_applied
            else Z_full
        )
        display_x_axis, display_y_axis, Z_plot = _display_surface(
            full_x_axis,
            full_y_axis,
            Z_visual_full,
            display_quality,
        )
        X_plot, Y_plot = np.meshgrid(display_x_axis, display_y_axis)
        color_scale = _colorscale(cmap_name, contrast)
        z_min, z_max = _effective_color_limits(
            Z_visual_full,
            brightness,
            color_range_mode,
            color_vmin,
            color_vmax,
        )
        (
            main_fig,
            _figure_x_axis,
            _figure_y_axis,
            x_min,
            x_max,
            y_min,
            y_max,
            _,
        ) = _scalar_field_figure(
            X=X_plot,
            Y=Y_plot,
            Z=Z_plot,
            colorscale=color_scale,
            z_min=z_min,
            z_max=z_max,
            show_blur=False,
            gaussian_sigma=gaussian_sigma,
        )

        overlay_warning = None
        try:
            overlay = _overlay_traces(
                structure=structure,
                x_min=x_min,
                x_max=x_max,
                y_min=y_min,
                y_max=y_max,
                show_atoms=show_atoms,
                show_grid=show_grid,
                layers=layers,
                radius_type=atom_radius_type,
                radius_scale=atom_radius_scale,
                grid_line_width=grid_line_width,
                grid_line_color=grid_line_color,
            )
        except Exception as exc:
            LOGGER.exception("Overlay rendering failed.")
            overlay = []
            overlay_warning = (
                f"Atom/grid overlay failed: {type(exc).__name__}: {exc}. "
                "Main image was rendered without overlays."
            )
        for trace in overlay:
            main_fig.add_trace(trace)

        try:
            scale_val = float(scalebar_length)
        except Exception:
            scale_val = 0.0
        if scale_val > 0.0 and scale_val < (x_max - x_min):
            x_end = x_max - 0.06 * (x_max - x_min)
            x_start = x_end - scale_val
            y_pos = y_min + 0.08 * (y_max - y_min)
            # Black outline for visibility.
            main_fig.add_shape(
                type="line",
                x0=x_start,
                y0=y_pos,
                x1=x_end,
                y1=y_pos,
                line={"color": "black", "width": 8},
            )
            main_fig.add_shape(
                type="line",
                x0=x_start,
                y0=y_pos,
                x1=x_end,
                y1=y_pos,
                line={"color": "white", "width": 4},
            )

        line_on = bool(line_enabled)
        raw_points = _clip_line_points(line_points, full_x_axis, full_y_axis)
        if line_on and line_magnet_enabled:
            points = [
                _snap_point_to_local_optimum(
                    full_x_axis,
                    full_y_axis,
                    Z_visual_full,
                    p,
                )
                for p in raw_points
            ]
        elif line_on:
            points = raw_points
        else:
            points = []

        _add_line_selection_traces(main_fig, points)

        surface_key = _store_surface_data(
            {
                "x_axis": full_x_axis,
                "y_axis": full_y_axis,
                "z": Z_full,
                "analysis_z": Z_visual_full,
                "display_x_axis": display_x_axis,
                "display_y_axis": display_y_axis,
                "display_z": Z_plot,
                "overlay": [trace.to_plotly_json() for trace in overlay],
                "shapes": [
                    shape.to_plotly_json() for shape in (main_fig.layout.shapes or ())
                ],
                "blur_applied": bool(blur_applied),
                "gaussian_sigma": float(gaussian_sigma or 0.0),
                "x_min": float(x_min),
                "x_max": float(x_max),
                "y_min": float(y_min),
                "y_max": float(y_max),
            }
        )

        simulation_title = {
            "STM": "STM surface map",
            "PHI_APP": r"$\Phi_{\rm app}\; \mathrm{map}$",
            "LWF": r"$\Phi_{\rm loc}\; \mathrm{map}$",
        }.get(simulation, "Surface map")
        main_fig.update_layout(
            template="plotly_white",
            title=simulation_title,
            height=580,
            margin={"l": 60, "r": 20, "t": 48, "b": 36},
            uirevision=f"{volume_real}:{simulation}:{effective_source}:{mode}",
        )
        main_fig.update_xaxes(
            title="x (A)",
            showgrid=False,
            constrain="domain",
            range=[x_min, x_max],
            autorange=False,
        )
        main_fig.update_yaxes(
            title="y (A)",
            showgrid=False,
            scaleanchor="x",
            scaleratio=1.0,
            range=[y_min, y_max],
            autorange=False,
        )

        line_fig = no_update
        fft_fig = no_update

        warnings = []
        if structure is None and (show_atoms or show_grid):
            warnings.append(
                "Structure file is missing or invalid, so atom/grid overlay is skipped."
            )
        if structure is None and show_repeat:
            warnings.append(
                "Structure-based iteration correction is unavailable; fallback iteration is used."
            )
        if overlay_warning:
            warnings.append(overlay_warning)
        finite_vals = np.asarray(Z_visual_full, dtype=float)
        finite_vals = finite_vals[np.isfinite(finite_vals)]
        if finite_vals.size == 0:
            finite_min = float("nan")
            finite_max = float("nan")
            warnings.append(
                "Rendered values are non-finite. "
                "Try another mode/isovalue or check CURRENT/PARCHG/LOCPOT data."
            )
        else:
            finite_min = float(np.min(finite_vals))
            finite_max = float(np.max(finite_vals))
        status_lines = [
            f"Rendered {_simulation_display_name(simulation)} | {effective_source} | {mode}",
            f"File: {os.path.basename(volume_real)}",
            f"Physical grid: {Z_full.shape[1]} x {Z_full.shape[0]}",
            f"Display grid: {Z_plot.shape[1]} x {Z_plot.shape[0]} ({display_quality})",
            f"XY interpolation factor: x{upsample_factor}",
            f"x: {x_min:.3f} .. {x_max:.3f} A",
            f"y: {y_min:.3f} .. {y_max:.3f} A",
            f"value: {finite_min:.6g} .. {finite_max:.6g} {value_unit}",
        ]
        if simulation == "PHI_APP":
            fit_center = float(iso_slider)
            status_lines.append(
                "Fit z-window: "
                f"{fit_center - fit_radius:.3f} .. "
                f"{fit_center + fit_radius:.3f} A "
                f"(center {fit_center:.3f} A, half-width {fit_radius:.3f} A)"
            )
        if stm_height_window is not None:
            status_lines.append(
                "CURRENT z-window: "
                f"{stm_height_window[0]:.3f} .. {stm_height_window[1]:.3f} A "
                "(no extrapolation)"
            )
        if fermi_level_used is not None:
            status_lines.append(
                f"E_F: {fermi_level_used:.8g} eV ({fermi_level_source})"
            )
        if structure_source is not None:
            status_lines.append(f"Structure: {structure_source}")
        status_color = "green"
        if warnings:
            status_lines.extend(warnings)
            status_color = "yellow"
        status_message = _status_text("OK", "\n".join(status_lines))
        _terminal_log(
            "Rendered "
            f"{_simulation_display_name(simulation)} | {effective_source} | {mode} "
            f"| file={volume_real} | shape={Z.shape}"
        )

        meta = {
            "colorscale": color_scale,
            "zmin": float(z_min),
            "zmax": float(z_max),
            "unit": value_unit,
            "title": cbar_title,
            "x_min": float(x_min),
            "x_max": float(x_max),
            "y_min": float(y_min),
            "y_max": float(y_max),
            "data_min": float(finite_min),
            "data_max": float(finite_max),
            "surface_key": surface_key,
            "blur_applied": bool(blur_applied),
            "physical_shape": [int(Z_full.shape[0]), int(Z_full.shape[1])],
            "display_shape": [int(Z_plot.shape[0]), int(Z_plot.shape[1])],
        }
        cbar_fig = _colorbar_figure(color_scale, z_min, z_max, cbar_title)
        return (
            main_fig,
            cbar_fig,
            line_fig,
            fft_fig,
            meta,
            surface_key,
            status_message,
            status_color,
            _status_ui(True),
        )

    @app.callback(
        Output("line-display-points-store", "data"),
        Output("line-profile-graph", "figure", allow_duplicate=True),
        Input("line-points-store", "data"),
        Input("line-profile-enabled", "checked"),
        Input("line-magnet-enabled", "checked"),
        Input("surface-key-store", "data"),
        State("plot-meta-store", "data"),
        State("brightness-slider", "value"),
        State("color-range-mode", "value"),
        State("color-vmin-input", "value"),
        State("color-vmax-input", "value"),
        prevent_initial_call=True,
    )
    def _refresh_line_analysis(
        line_points,
        line_enabled,
        line_magnet_enabled,
        surface_key,
        meta,
        brightness,
        color_range_mode,
        color_vmin,
        color_vmax,
    ):
        surface = _get_surface_data(surface_key)
        if surface is None:
            return [], _empty_figure(
                "Line profile",
                "Render a surface map first.",
                height=280,
            )

        raw_points = (
            _clip_line_points(
                line_points,
                surface["x_axis"],
                surface["y_axis"],
            )
            if line_enabled
            else []
        )
        if not line_enabled:
            return [], _empty_figure(
                "Line profile",
                "Line profile is disabled.",
                height=280,
            )
        if not raw_points:
            return [], _empty_figure(
                "Line profile",
                "Select P1 on the surface map.",
                height=280,
            )

        if line_magnet_enabled:
            points = [
                _snap_point_to_local_optimum(
                    surface["x_axis"],
                    surface["y_axis"],
                    _surface_analysis_values(surface),
                    point,
                )
                for point in raw_points
            ]
        else:
            points = raw_points

        if len(points) == 1:
            return points, _empty_figure(
                "Line profile",
                "P1 selected. Select P2 to draw the line.",
                height=280,
            )

        if len(points) == 2:
            y_range = _effective_color_limits(
                _surface_analysis_values(surface),
                brightness,
                color_range_mode,
                color_vmin,
                color_vmax,
            )
            line_figure = _line_profile(
                X=surface["x_axis"],
                Y=surface["y_axis"],
                Z=_surface_analysis_values(surface),
                p1=points[0],
                p2=points[1],
                title="Line profile",
                y_unit=(meta or {}).get("unit", "a.u."),
                y_range=y_range,
            )
        else:
            line_figure = _empty_figure("Line profile", "Select P1.", height=280)
        return points, line_figure

    @app.callback(
        Output("save-line-profile-button", "disabled"),
        Input("line-profile-graph", "figure"),
    )
    def _line_profile_export_disabled(figure):
        try:
            _line_profile_csv(figure or {})
        except ValueError:
            return True
        return False

    @app.callback(
        Output("download-line-profile", "data"),
        Input("save-line-profile-button", "n_clicks"),
        State("line-profile-graph", "figure"),
        prevent_initial_call=True,
    )
    def _save_line_profile_csv(n_clicks, figure):
        if not n_clicks:
            raise PreventUpdate
        try:
            payload = _line_profile_csv(figure or {})
        except ValueError:
            raise PreventUpdate
        _terminal_log("Line-profile CSV prepared for download.")
        return dcc.send_string(payload, "line_profile.csv")

    @app.callback(
        Output("fft-graph", "figure", allow_duplicate=True),
        Input("fft-enabled", "checked"),
        Input("surface-key-store", "data"),
        Input("analysis-tabs", "value"),
        State("plot-meta-store", "data"),
        prevent_initial_call=True,
    )
    def _refresh_fft_analysis(fft_enabled, surface_key, active_tab, meta):
        if not fft_enabled:
            return _empty_figure("2D FFT", "FFT is disabled.", height=280)
        if active_tab != "fft":
            raise PreventUpdate
        surface = _get_surface_data(surface_key)
        if not meta or surface is None:
            return _empty_figure("2D FFT", "Render a surface map first.", height=280)
        try:
            return _fft_figure(
                _surface_analysis_values(surface),
                surface["x_axis"],
                surface["y_axis"],
            )
        except Exception:
            return _empty_figure(
                "2D FFT",
                "Rendered scalar data is unavailable.",
                height=280,
            )

    @app.callback(
        Output("main-image", "figure", allow_duplicate=True),
        Input("main-image", "relayoutData"),
        State("main-image", "figure"),
        State("plot-meta-store", "data"),
        prevent_initial_call=True,
    )
    def _enforce_main_axes_bounds(relayout_data, figure_dict, meta):
        if not relayout_data or not figure_dict or not meta:
            raise PreventUpdate
        if not _has_scalar_surface_trace(figure_dict):
            raise PreventUpdate
        if not isinstance(relayout_data, dict):
            raise PreventUpdate

        reset_requested = (
            relayout_data.get("xaxis.autorange") is True
            or relayout_data.get("yaxis.autorange") is True
            or relayout_data.get("autosize") is True
        )
        if not reset_requested:
            raise PreventUpdate

        try:
            x_min = float(meta.get("x_min"))
            x_max = float(meta.get("x_max"))
            y_min = float(meta.get("y_min"))
            y_max = float(meta.get("y_max"))
        except Exception:
            raise PreventUpdate

        fig = go.Figure(figure_dict)
        fig.update_xaxes(
            range=[x_min, x_max],
            autorange=False,
            constrain="domain",
        )
        fig.update_yaxes(
            range=[y_min, y_max],
            autorange=False,
            scaleanchor="x",
            scaleratio=1.0,
        )
        return fig

    @app.callback(
        Output("main-image", "figure", allow_duplicate=True),
        Input("set-default-button", "n_clicks"),
        State("main-image", "figure"),
        State("plot-meta-store", "data"),
        prevent_initial_call=True,
    )
    def _enforce_main_axes_on_default_reset(n_clicks, figure_dict, meta):
        if not n_clicks or not figure_dict or not meta:
            raise PreventUpdate
        if not _has_scalar_surface_trace(figure_dict):
            raise PreventUpdate
        try:
            x_min = float(meta.get("x_min"))
            x_max = float(meta.get("x_max"))
            y_min = float(meta.get("y_min"))
            y_max = float(meta.get("y_max"))
        except Exception:
            raise PreventUpdate

        fig = go.Figure(figure_dict)
        fig.update_xaxes(
            range=[x_min, x_max],
            autorange=False,
            constrain="domain",
        )
        fig.update_yaxes(
            range=[y_min, y_max],
            autorange=False,
            scaleanchor="x",
            scaleratio=1.0,
        )
        return fig

    @app.callback(
        Output("download-image", "data"),
        Output("render-status", "children", allow_duplicate=True),
        Output("render-status", "color", allow_duplicate=True),
        Output("render-status-wrap", "style", allow_duplicate=True),
        Input("save-image-button", "n_clicks"),
        State("surface-key-store", "data"),
        State("appearance-store", "data"),
        State("brightness-slider", "value"),
        State("color-range-mode", "value"),
        State("color-vmin-input", "value"),
        State("color-vmax-input", "value"),
        State("export-format", "value"),
        State("export-width", "value"),
        State("export-height", "value"),
        State("line-display-points-store", "data"),
        running=[
            (Output("save-image-button", "loading"), True, False),
        ],
        prevent_initial_call=True,
    )
    def _save_image(
        n_clicks,
        surface_key,
        appearance,
        brightness,
        color_range_mode,
        color_vmin,
        color_vmax,
        export_format,
        export_width,
        export_height,
        line_points,
    ):
        _terminal_log(f"Save image clicked: n_clicks={n_clicks}")
        surface = _get_surface_data(surface_key)
        if not n_clicks or surface is None:
            _terminal_log(
                "Save image skipped: no rendered figure is available.", level="warning"
            )
            msg = _status_text("WAIT", "Save image skipped.\nRender an image first.")
            return no_update, msg, "yellow", _status_ui(True)
        try:
            appearance = appearance or {}
            colorscale = appearance.get("colorscale") or _colorscale("afmhot", 0.0)
            cmap_name = appearance.get("cmap_name", "afmhot")
            contrast = float(appearance.get("contrast", 0.0))
            z_min, z_max = _effective_color_limits(
                _surface_analysis_values(surface),
                brightness,
                color_range_mode,
                color_vmin,
                color_vmax,
            )
            figure = _cached_surface_figure(
                surface=surface,
                colorscale=colorscale,
                z_min=z_min,
                z_max=z_max,
                cmap_name=cmap_name,
                contrast=contrast,
                line_points=line_points,
            )
            x_range = [float(surface["x_min"]), float(surface["x_max"])]
            y_range = [float(surface["y_min"]), float(surface["y_max"])]
            width_px, height_px = _validate_export_dimensions(
                export_width,
                export_height,
            )
            image_format, filename, _mime_type = _export_file_details(
                export_format,
                f"autobskan_image_{width_px}x{height_px}",
            )

            figure.update_layout(
                margin={"l": 0, "r": 0, "t": 0, "b": 0, "pad": 0},
                title=None,
                width=width_px,
                height=height_px,
                paper_bgcolor="rgba(0,0,0,0)",
                plot_bgcolor="rgba(0,0,0,0)",
                showlegend=False,
            )
            figure.update_xaxes(
                title=None,
                visible=False,
                showticklabels=False,
                showgrid=False,
                showline=False,
                zeroline=False,
                range=x_range,
                autorange=False,
                fixedrange=True,
                constrain="domain",
            )
            figure.update_yaxes(
                title=None,
                visible=False,
                showticklabels=False,
                showgrid=False,
                showline=False,
                zeroline=False,
                range=y_range,
                autorange=False,
                fixedrange=True,
                scaleanchor="x",
                scaleratio=1.0,
                constrain="domain",
            )
            image_bytes = _render_export_figure(
                figure,
                image_format,
                width_px,
                height_px,
            )
        except Exception as exc:
            _terminal_log(
                f"Save image failed: {type(exc).__name__}: {exc}",
                level="error",
            )
            msg = _status_text(
                "ERROR",
                "Save image failed.\n"
                f"{type(exc).__name__}: {exc}\n"
                f"{_static_export_help()}",
            )
            return no_update, msg, "red", _status_ui(True)
        _terminal_log(
            f"Save image prepared: {width_px}x{height_px} {image_format.upper()}."
        )
        msg = _status_text(
            "OK",
            f"Image export ready ({width_px} x {height_px} px, "
            f"{image_format.upper()}).\nDownloading {filename}",
        )
        return (
            dcc.send_bytes(
                lambda fileobj: fileobj.write(image_bytes), filename
            ),
            msg,
            "green",
            _status_ui(True),
        )

    @app.callback(
        Output("download-colorbar", "data"),
        Output("render-status", "children", allow_duplicate=True),
        Output("render-status", "color", allow_duplicate=True),
        Output("render-status-wrap", "style", allow_duplicate=True),
        Input("save-cbar-button", "n_clicks"),
        State("plot-meta-store", "data"),
        State("appearance-store", "data"),
        State("brightness-slider", "value"),
        State("color-range-mode", "value"),
        State("color-vmin-input", "value"),
        State("color-vmax-input", "value"),
        State("export-format", "value"),
        State("colorbar-width", "value"),
        State("colorbar-height", "value"),
        running=[
            (Output("save-cbar-button", "loading"), True, False),
        ],
        prevent_initial_call=True,
    )
    def _save_colorbar(
        n_clicks,
        meta,
        appearance,
        brightness,
        color_range_mode,
        color_vmin,
        color_vmax,
        export_format,
        colorbar_width,
        colorbar_height,
    ):
        _terminal_log(f"Save colorbar clicked: n_clicks={n_clicks}")
        if not n_clicks or not meta:
            _terminal_log(
                "Save colorbar skipped: no plot meta data is available.",
                level="warning",
            )
            msg = _status_text("WAIT", "Save colorbar skipped.\nRender an image first.")
            return no_update, msg, "yellow", _status_ui(True)
        appearance = appearance or {}
        colorscale = appearance.get("colorscale") or meta.get("colorscale")
        data_min = meta.get("data_min")
        data_max = meta.get("data_max")
        if data_min is not None and data_max is not None:
            zmin, zmax = _effective_color_limits(
                np.asarray([data_min, data_max], dtype=float),
                brightness,
                color_range_mode,
                color_vmin,
                color_vmax,
            )
        else:
            zmin = meta.get("zmin")
            zmax = meta.get("zmax")
        if colorscale is None or zmin is None or zmax is None:
            _terminal_log(
                "Save colorbar skipped: incomplete colorbar metadata.", level="warning"
            )
            msg = _status_text(
                "ERROR", "Save colorbar failed.\nIncomplete colorbar metadata."
            )
            return no_update, msg, "red", _status_ui(True)

        try:
            width_px, height_px = _validate_export_dimensions(
                colorbar_width,
                colorbar_height,
            )
            image_format, filename, _mime_type = _export_file_details(
                export_format,
                f"autobskan_colorbar_{width_px}x{height_px}",
            )
            fig = _export_colorbar_figure(colorscale, zmin, zmax)
            image_bytes = _render_export_figure(
                fig,
                image_format,
                width_px,
                height_px,
                opaque=True,
            )
        except Exception as exc:
            _terminal_log(
                f"Save colorbar failed: {type(exc).__name__}: {exc}",
                level="error",
            )
            msg = _status_text(
                "ERROR",
                "Save colorbar failed.\n"
                f"{type(exc).__name__}: {exc}\n"
                f"{_static_export_help()}",
            )
            return no_update, msg, "red", _status_ui(True)
        _terminal_log(
            f"Save colorbar prepared: {width_px}x{height_px} "
            f"{image_format.upper()}."
        )
        msg = _status_text(
            "OK",
            f"Colorbar export ready ({width_px} x {height_px} px, "
            f"{image_format.upper()}).\nDownloading {filename}",
        )
        return (
            dcc.send_bytes(
                lambda fileobj: fileobj.write(image_bytes),
                filename,
            ),
            msg,
            "green",
            _status_ui(True),
        )

    @app.callback(
        Output("download-fft-image", "data"),
        Output("render-status", "children", allow_duplicate=True),
        Output("render-status", "color", allow_duplicate=True),
        Output("render-status-wrap", "style", allow_duplicate=True),
        Input("save-fft-button", "n_clicks"),
        State("fft-graph", "figure"),
        prevent_initial_call=True,
    )
    def _save_fft_image(n_clicks, figure_dict):
        _terminal_log(f"Save FFT clicked: n_clicks={n_clicks}")
        if not n_clicks or not figure_dict:
            msg = _status_text("WAIT", "Save FFT skipped.\nRender FFT first.")
            return no_update, msg, "yellow", _status_ui(True)
        try:
            figure = go.Figure(figure_dict)
            png_bytes = pio.to_image(
                figure,
                format="png",
                width=1200,
                height=900,
                scale=2,
            )
        except Exception as exc:
            msg = _status_text(
                "ERROR",
                "Save FFT failed.\n"
                f"{type(exc).__name__}: {exc}\n"
                f"{_static_export_help()}",
            )
            return no_update, msg, "red", _status_ui(True)
        msg = _status_text("OK", "FFT export ready.\nDownloading autobskan_fft.png")
        return (
            dcc.send_bytes(
                lambda fileobj: fileobj.write(png_bytes), "autobskan_fft.png"
            ),
            msg,
            "green",
            _status_ui(True),
        )

    @app.callback(
        Output("download-bskanin", "data"),
        Output("render-status", "children", allow_duplicate=True),
        Output("render-status", "color", allow_duplicate=True),
        Output("render-status-wrap", "style", allow_duplicate=True),
        Input("export-bskanin-button", "n_clicks"),
        State("simulation-type", "value"),
        State("input-source", "value"),
        State("stm-mode", "value"),
        State("fermi-level-input", "value"),
        State("volume-path", "value"),
        State("volume-file-store", "data"),
        State("structure-path", "value"),
        State("structure-file-store", "data"),
        State("gamma-input", "value"),
        State("colormap-dropdown", "value"),
        State("brightness-slider", "value"),
        State("contrast-slider", "value"),
        State("isosurface-slider", "value"),
        State("fit-radius-slider", "value"),
        State("layers-input", "value"),
        State("view-options", "value"),
        State("repeat-x", "value"),
        State("repeat-y", "value"),
        State("atom-radius-type", "value"),
        State("atom-radius-scale", "value"),
        State("gaussian-sigma", "value"),
        prevent_initial_call=True,
    )
    def _export_bskanin(
        n_clicks,
        simulation,
        input_source,
        stm_mode,
        fermi_level_input,
        volume_path,
        volume_store_path,
        structure_path,
        structure_store_path,
        gamma_value,
        cmap_name,
        brightness,
        contrast,
        iso_slider,
        fit_radius,
        layers,
        view_options,
        repeat_x,
        repeat_y,
        atom_radius_type,
        atom_radius_scale,
        gaussian_sigma,
    ):
        if not n_clicks:
            raise PreventUpdate

        mode = (
            "CONSTANT HEIGHT"
            if simulation in {"PHI_APP", "LWF"}
            else stm_mode
        )
        opts = set(view_options or [])
        show_blur = "show_blurred" in opts
        nx, ny = _parse_repeat(repeat_x, repeat_y)

        current_target = _effective_path(volume_path, volume_store_path)
        poscar_target = _effective_path(structure_path, structure_store_path)
        gamma = float(gamma_value if gamma_value is not None else 90.0)
        iso_value = (
            10.0 ** float(iso_slider)
            if simulation == "STM" and mode == "CONSTANT CURRENT"
            else float(iso_slider)
        )

        try:
            resolved_target = _safe_path(current_target)
            if resolved_target is None:
                raise ValueError(
                    "Select an accessible CURRENT or VASP volumetric file before exporting."
                )
            _validate_simulation_input(
                simulation,
                input_source,
                resolved_target,
            )
            fermi_level_export = None
            if simulation == "LWF":
                fermi_level_export, _fermi_source = _resolve_lwf_fermi_level(
                    resolved_target,
                    fermi_level_input,
                )
            content = _bskanin_export_text(
                simulation=simulation,
                input_source=input_source,
                volume_path=resolved_target,
                mode=mode,
                iso_value=iso_value,
                fit_radius=fit_radius,
                cmap_name=cmap_name,
                contrast=contrast,
                brightness=brightness,
                poscar_path=poscar_target,
                iteration=(nx, ny),
                gaussian_sigma=float(gaussian_sigma) if show_blur else 0.0,
                gamma=gamma,
                display_atoms="show_atoms" in opts,
                layers=layers,
                radius_type=atom_radius_type,
                size_ratio=atom_radius_scale,
                display_cell="show_grids" in opts,
                fermi_level=fermi_level_export,
            )
        except (OSError, TypeError, ValueError) as exc:
            msg = _status_text("ERROR", f"Cannot export settings.\n{exc}")
            return no_update, msg, "red", _status_ui(True)
        status_color = "green"
        msg = _status_text(
            "OK",
            "Exported GUI settings to bskan.in format.\nDownloading bskan_gui_export.in",
        )
        return (
            dcc.send_string(content, "bskan_gui_export.in"),
            msg,
            status_color,
            _status_ui(True),
        )

    return app

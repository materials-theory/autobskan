import os
import re

import numpy as np

from autobskan.image import AR, stmplot


APPARENT_BARRIER_COEFFICIENT = 0.952495
_FLOAT_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
_OUTCAR_FERMI_PATTERNS = (
    re.compile(rf"E-fermi\s*:\s*({_FLOAT_PATTERN})", re.IGNORECASE),
    re.compile(rf"Fermi\s+energy\s*:?\s*({_FLOAT_PATTERN})", re.IGNORECASE),
)


def cell_c_length(volumetric):
    cell = AR.to_new_cell_onlycell(np.asarray(volumetric.cell, dtype=float))
    c_length = float(cell[2][2])
    if c_length <= 0:
        raise ValueError("Volumetric c lattice length should be positive.")
    return c_length


def topmost_z(structure=None, volumetric=None):
    if structure is not None:
        atoms = structure
    elif volumetric is not None and getattr(volumetric, "atoms", None) is not None:
        atoms = AR.to_new_cell(volumetric.atoms)
    else:
        return 0.0
    return float(np.max(atoms.get_positions()[:, 2]))


def current_height_window(current):
    dz_shift = -float(current.topmost_atom)
    return dz_shift, float(current.cell[2, 2] + dz_shift)


def current_z_axis(current):
    nz = int(current.grids[-1])
    dz = float(current.cell[2, 2]) / float(nz)
    z_min, _ = current_height_window(current)
    return np.arange(nz, dtype=float) * dz + z_min


def vasp_z_axis(volumetric):
    nz = int(np.asarray(volumetric.pot).shape[0])
    c_length = cell_c_length(volumetric)
    return np.arange(nz, dtype=float) * c_length / float(nz)


def surface_potential_slice(volumetric, height, topmost=0.0):
    potential = np.asarray(volumetric.pot, dtype=float)
    nz = potential.shape[0]
    c_length = cell_c_length(volumetric)
    z_abs = float(topmost) + float(height)
    z_abs = min(max(z_abs, 0.0), c_length * (nz - 1) / nz)
    dz = c_length / nz
    z_index = z_abs / dz
    i_low = int(np.floor(z_index))
    if i_low >= nz - 1:
        i_low = nz - 2
        frac = 1.0
    else:
        frac = z_index - i_low
    return (1.0 - frac) * potential[i_low, :, :] + frac * potential[i_low + 1, :, :]


def local_work_function_slice(
    volumetric,
    height,
    fermi_level,
    topmost=0.0,
):
    """Return a finite-height LOCPOT map referenced to the Fermi level."""
    fermi = float(fermi_level)
    if not np.isfinite(fermi):
        raise ValueError("Fermi level should be a finite value in eV.")
    return surface_potential_slice(
        volumetric,
        height,
        topmost=topmost,
    ) - fermi


def fermi_energy_from_outcar(path):
    """Read the last reported VASP Fermi energy from an OUTCAR file."""
    fermi = None
    with open(path, "r", encoding="utf-8", errors="replace") as fileobj:
        for line in fileobj:
            for pattern in _OUTCAR_FERMI_PATTERNS:
                match = pattern.search(line)
                if match is not None:
                    fermi = float(match.group(1))
                    break
    if fermi is None:
        raise ValueError(f"No Fermi energy was found in OUTCAR: {path}")
    return fermi


def sibling_outcar_path(volume_path):
    """Return a sibling OUTCAR path for a local VASP volumetric file, if present."""
    if not isinstance(volume_path, str) or not volume_path.strip():
        return None
    candidate = os.path.join(
        os.path.dirname(os.path.realpath(volume_path)),
        "OUTCAR",
    )
    return candidate if os.path.isfile(candidate) else None


def resolve_fermi_level(volume_path, fermi_level=None):
    """Resolve an explicit Fermi level or read it from a sibling OUTCAR."""
    if fermi_level is not None and str(fermi_level).strip().upper() not in {
        "",
        "AUTO",
        "NONE",
    }:
        value = float(fermi_level)
        if not np.isfinite(value):
            raise ValueError("Fermi level should be a finite value in eV.")
        return value

    outcar_path = sibling_outcar_path(volume_path)
    if outcar_path is None:
        raise ValueError(
            "LWF requires E_F: enter the Fermi level in eV or place OUTCAR "
            "in the same directory as LOCPOT."
        )
    return fermi_energy_from_outcar(outcar_path)


def _fit_log_decay_uniform(values, z_axis, center, fit_radius):
    values = np.asarray(values, dtype=float)
    z_axis = np.asarray(z_axis, dtype=float)
    lower = float(center) - float(fit_radius)
    upper = float(center) + float(fit_radius)
    fit_mask = (z_axis >= lower) & (z_axis <= upper)
    if np.count_nonzero(fit_mask) < 3:
        nearest = int(np.clip(np.round(np.searchsorted(z_axis, center)), 0, len(z_axis) - 1))
        start = max(nearest - 2, 0)
        stop = min(nearest + 3, len(z_axis))
        fit_mask = np.zeros(len(z_axis), dtype=bool)
        fit_mask[start:stop] = True
    if np.count_nonzero(fit_mask) < 3:
        raise ValueError("Apparent-barrier fit needs at least three z grid points.")

    z_fit = z_axis[fit_mask]
    data_fit = values[fit_mask, :, :]
    log_data = np.full_like(data_fit, np.nan, dtype=float)
    positive = data_fit > 0.0
    log_data[positive] = np.log(data_fit[positive])
    valid = np.isfinite(log_data)
    n_valid = np.sum(valid, axis=0).astype(float)
    z3 = z_fit[:, None, None]
    with np.errstate(invalid="ignore", divide="ignore"):
        sum_z = np.sum(np.where(valid, z3, 0.0), axis=0)
        sum_y = np.nansum(log_data, axis=0)
        sum_zz = np.sum(np.where(valid, z3 * z3, 0.0), axis=0)
        sum_zy = np.nansum(log_data * z3, axis=0)
        denom = n_valid * sum_zz - sum_z * sum_z
        slope = (n_valid * sum_zy - sum_z * sum_y) / denom
        barrier = APPARENT_BARRIER_COEFFICIENT * slope * slope
    barrier[(n_valid < 3.0) | ~np.isfinite(barrier)] = np.nan
    if not np.any(np.isfinite(barrier)):
        raise ValueError("Apparent-barrier fit failed: no positive values in the z window.")
    return barrier


def _fit_log_decay_variable(values, z_axis, centers, fit_radius):
    values = np.asarray(values, dtype=float)
    z_axis = np.asarray(z_axis, dtype=float)
    centers = np.asarray(centers, dtype=float)
    result = np.full(values.shape[1:], np.nan, dtype=float)
    ny, nx = result.shape
    for iy in range(ny):
        for ix in range(nx):
            center = centers[iy, ix]
            if not np.isfinite(center):
                continue
            fit_mask = (z_axis >= center - fit_radius) & (z_axis <= center + fit_radius)
            if np.count_nonzero(fit_mask) < 3:
                nearest = int(np.argmin(np.abs(z_axis - center)))
                start = max(nearest - 2, 0)
                stop = min(nearest + 3, len(z_axis))
                fit_mask = np.zeros(len(z_axis), dtype=bool)
                fit_mask[start:stop] = True
            y = values[fit_mask, iy, ix]
            x = z_axis[fit_mask]
            positive = y > 0.0
            if np.count_nonzero(positive) < 3:
                continue
            x = x[positive]
            log_y = np.log(y[positive])
            slope, _intercept = np.polyfit(x, log_y, 1)
            result[iy, ix] = APPARENT_BARRIER_COEFFICIENT * slope * slope
    if not np.any(np.isfinite(result)):
        raise ValueError("Apparent-barrier fit failed: no positive values in the z window.")
    return result


def apparent_barrier_from_vasp_density(volumetric, height, topmost=0.0, fit_radius=0.5):
    density = np.asarray(volumetric.pot, dtype=float)
    z_axis = vasp_z_axis(volumetric)
    z_center = float(topmost) + float(height)
    z_center = min(max(z_center, 0.0), float(z_axis[-1]))
    return _fit_log_decay_uniform(density, z_axis, z_center, fit_radius)


def apparent_barrier_from_current(
    current,
    height=None,
    fit_radius=0.5,
    isosurface=None,
    reference="CONSTANT HEIGHT",
):
    values = np.asarray(current.cur_3d, dtype=float)
    z_axis = current_z_axis(current)
    reference = str(reference).upper()
    if reference.startswith("CONSTANT CURRENT"):
        if isosurface is None:
            raise ValueError(
                "CURRENT apparent-barrier analysis with a constant-current "
                "reference needs an isosurface."
            )
        centers = stmplot.get_surface_bskan(
            current=current,
            isosurface=float(isosurface),
            constant_current=True,
            return_index=False,
        )
        return _fit_log_decay_variable(values, z_axis, centers, fit_radius)

    if height is None:
        raise ValueError(
            "CURRENT apparent-barrier analysis with a constant-height "
            "reference needs a height."
        )
    z_min, z_max = current_height_window(current)
    z_center = float(np.clip(float(height), z_min, z_max))
    return _fit_log_decay_uniform(values, z_axis, z_center, fit_radius)

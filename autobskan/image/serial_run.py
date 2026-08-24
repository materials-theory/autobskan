import copy
import os

# from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndimage
from matplotlib.figure import Figure

from autobskan.image import AR, lwfplot, stmplot
from autobskan.image.post_processing import array_iter
from autobskan.image.stmplot import Current, get_surface_bskan, plot_atoms, plot_cell
from autobskan.io.input import Options


def _prepare_output_dir(image_dir):
    before_enter_image_dir = os.getcwd()
    if not os.path.exists(image_dir):
        os.mkdir(image_dir)
    os.chdir(image_dir)
    return before_enter_image_dir


def _save_scalar_map(X, Y, Z, filename, options, unit_label=None):
    brightness = options["BRIGHTNESS"]
    resolution = options["CONTOUR_RESOLUTION"]
    brighter, darker = (0, -brightness) if brightness < 0 else (brightness, 0)
    cmap = copy.deepcopy(plt.get_cmap(options["CMAP"]))
    cmap.set_gamma(1 + options["CONTRAST"])
    z_span = np.nanmax(Z) - np.nanmin(Z)
    if not np.isfinite(z_span) or z_span == 0:
        z_span = 1e-12
    levels = np.linspace(
        np.nanmin(Z) - brighter * z_span,
        np.nanmax(Z) + darker * z_span,
        int(resolution * (1 + abs(brightness))),
    )
    fig, ax = plt.subplots()
    con = ax.contourf(X, Y, Z, cmap=cmap, levels=levels)
    if options["DISPLAY_CBAR"]:
        cbar = plt.colorbar(con, ax=ax)
        if unit_label:
            cbar.set_label(unit_label)
    ax.axis("off")
    ax.set_aspect("equal")
    ax.set_xlim(np.min(X), np.max(X))
    ax.set_ylim(np.min(Y), np.max(Y))
    fig.savefig(f"{filename}.png", dpi=300, bbox_inches="tight", pad_inches=0)
    plt.close(fig)


def _height_values(options, h_min, h_max):
    if options["ISO_AUTO"] == "LOGSCALE":
        raise IOError("ISO_AUTO=LOGSCALE is not available for height mode")
    if options["ISO_AUTO"]:
        return list(np.linspace(h_min, h_max, options["ISO"]))
    return [options["ISO"]] if type(options["ISO"]) in [int, float] else options["ISO"]


def _constant_current_iso_values(options, iso_min, iso_max):
    """Return valid constant-current isosurfaces strictly inside the data range."""

    lower = float(iso_min)
    upper = float(iso_max)
    if not np.isfinite(lower) or not np.isfinite(upper):
        raise ValueError(
            f"Constant-current isosurface bounds must be finite, got [{lower}, {upper}]."
        )
    if lower <= 0.0 or upper <= 0.0:
        raise ValueError(
            "Constant-current automatic isosurfaces require positive CURRENT values; "
            f"got [{lower}, {upper}]."
        )
    if lower >= upper:
        raise ValueError(
            f"Constant-current isosurface range is empty: [{lower}, {upper}]."
        )

    iso_auto = options["ISO_AUTO"]
    if iso_auto is True or iso_auto == "LOGSCALE":
        first_decade = int(np.ceil(np.log10(lower)))
        last_decade = int(np.floor(np.log10(upper)))
        if first_decade <= last_decade:
            values = np.power(
                10.0,
                np.arange(first_decade, last_decade + 1, dtype=float),
            )
            values = values[(values > lower) & (values < upper)]
        else:
            values = np.empty(0, dtype=float)

        if values.size == 0:
            raise ValueError(
                "ISO_AUTO=LOGSCALE found no power of ten strictly inside "
                f"[{lower}, {upper}]. Use ISO_AUTO=TRUE for evenly spaced "
                "values or ISO_AUTO=FALSE for an explicit isosurface."
            )
    elif iso_auto == "LINEAR":
        count = int(options["ISO"])
        if count < 1:
            raise ValueError(f"ISO must be a positive image count, got {count}.")
        values = np.linspace(lower, upper, count + 2)[1:-1]
    elif iso_auto is False:
        explicit = options["ISO"]
        values = (
            np.asarray([explicit], dtype=float)
            if np.isscalar(explicit)
            else np.asarray(explicit, dtype=float).reshape(-1)
        )
    else:
        raise ValueError(
            "ISO_AUTO must be TRUE/LOGSCALE, LINEAR, or FALSE. "
            f"Received: {iso_auto}"
        )

    return [float(value) for value in values]


def _repeat_current_grid(current, Z_orig, options, poscar=None):
    nx, ny = options["ITERATION"]
    if nx * ny != 1:
        if poscar is not None:
            return array_iter(
                Z_orig,
                nx=nx,
                ny=ny,
                cell=poscar.cell,
                diagonal_transform_for_hexagonal=options["DIAGONAL_TRANSFORM_FOR_HEXAGONAL"],
            )
        return array_iter(
            Z_orig,
            nx=nx,
            ny=ny,
            gamma=options["GAMMA"],
            l_a=float(current.cellpar[0]),
            l_b=float(current.cellpar[1]),
        )
    x = np.linspace(0, float(current.cellpar[0]), current.grids[0])
    y = np.linspace(0, float(current.cellpar[1]), current.grids[1])
    X, Y = np.meshgrid(x, y)
    return X, Y, Z_orig


def _repeat_vasp_grid(volumetric, Z_orig, options, poscar=None):
    nx, ny = options["ITERATION"]
    cell = AR.to_new_cell_onlycell(np.asarray(poscar.cell if poscar is not None else volumetric.cell))
    return array_iter(
        Z_orig,
        nx=nx,
        ny=ny,
        cell=cell,
        input_filetype="VASP",
        verbose=False,
    )


def single_current(current:Current,
                   bskan_input:Options,
                   image_dir     = ".",
                   fig:Figure = None,
                   ):
    """
    * current : Current instance
    * bskan_input : Bskan_input instance (in io.py)
    * image_dir : where to store the images
    * fit : If you want to extract the axes, you can put plt.figure as an input -> subplots will be returned
    * blur : If True, gaussian filter will be applied as bskan_input.blur_sigma
    """

    _options = bskan_input.option_dict
    constant_current = _options["IMAGE_MODE"]=="CONSTANT CURRENT"
    poscar = None
    poscar_path = _options["POSCAR"]
    if isinstance(poscar_path, str) and os.path.exists(poscar_path):
        poscar = AR.read_vasp_robust(poscar_path)
    elif _options["DISPLAY_ATOMS"] or _options["DISPLAY_CELL"]:
        print(f"[WARN] POSCAR not found: {poscar_path}. Atom/cell overlays are disabled.")

    diagonal_transform_for_hexagonal = _options["DIAGONAL_TRANSFORM_FOR_HEXAGONAL"]
    display_atoms = _options["DISPLAY_ATOMS"]
    display_cell = _options["DISPLAY_CELL"]
    display_cbar = _options["DISPLAY_CBAR"]
    display_cbar_on_separate_plot = _options["DISPLAY_CBAR_ON_SEPARATE_PLOT"]
    blur = True if _options["BLUR_SIGMA"] != 0 else False

    nx, ny = _options["ITERATION"]
    plot_repeat = True if nx * ny != 1 else False

    atoms_info = _options["ATOMS_INFO"]

    return_subplots = True if fig is not None else False

    output_dir = None
    if not return_subplots:
        output_dir = os.path.abspath(os.path.expanduser(image_dir))
        os.makedirs(output_dir, exist_ok=True)

    if constant_current:
        iso_min = float(current.iso_min)
        iso_max = float(current.iso_max)
        iso_data = _constant_current_iso_values(_options, iso_min, iso_max)

    else:
        _zmin = current.topmost_atom
        _dz = current.cell[2,2]
        dz_shift = -_zmin
        iso_min, iso_max = dz_shift, _dz + dz_shift
        if _options["ISO_AUTO"] == "LOGSCALE":
            raise IOError("ISO_AUTO=LOGSCALE is not available for constant height mode")
        elif _options["ISO_AUTO"]:
            iso_data = list(np.linspace(iso_min, iso_max, _options["ISO"]))
        else:
            assert not _options["ISO_AUTO"], f"Invalid input of ISO_AUTO: {_options['ISO_AUTO']}"
            iso_data = [_options["ISO"]] if type(_options["ISO"]) in [int, float] else _options["ISO"]

    real_x, real_y = current.cellpar[:2]
    if return_subplots:
        if len(iso_data) == 1:
            axes = [fig.subplots(1, 1)]
        else:
            axes = fig.subplots(1, len(iso_data))
    generated_files = []
    for _i, iso in enumerate(iso_data):
        _fname = f"{iso}"
        try:
            Z_orig = get_surface_bskan(current = current, isosurface = iso,
                                       constant_current = constant_current, return_index = False)
        except IOError as exc:
            print(
                f"[WARN] Skipped ISO={iso:.10g} for "
                f"{getattr(current, 'filename', 'CURRENT')}: {exc}"
            )
            continue

        brightness = _options["BRIGHTNESS"]
        resolution = _options["CONTOUR_RESOLUTION"]

        ax_stm = plt.figure(figsize=(real_x, real_y)).gca() if not return_subplots else axes[_i]

        (brighter, darker) = (0, -brightness) if brightness < 0 else (brightness, 0)

        cmap = copy.deepcopy(plt.get_cmap(_options["CMAP"]))
        cmap.set_gamma(1+_options["CONTRAST"])

        levels = np.linspace(np.min(Z_orig) - brighter * (np.max(Z_orig)-np.min(Z_orig)),
                             np.max(Z_orig) + darker * (np.max(Z_orig)-np.min(Z_orig)),
                             int(resolution * (1 + abs(brightness))))

        if plot_repeat:
            if poscar is not None:
                X, Y, Z = array_iter(
                    Z_orig,
                    nx=nx,
                    ny=ny,
                    cell=poscar.cell,
                    diagonal_transform_for_hexagonal=diagonal_transform_for_hexagonal,
                )
            else:
                X, Y, Z = array_iter(
                    Z_orig,
                    nx=nx,
                    ny=ny,
                    gamma=_options["GAMMA"],
                    l_a=real_x,
                    l_b=real_y,
                )
            _fname += f"_{nx:d}x{ny:d}"
        else:
            x = np.linspace(0, real_x, current.grids[0])
            y = np.linspace(0, real_y, current.grids[1])
            X, Y = np.meshgrid(x, y)
            Z = Z_orig
            nx, ny = 1, 1

        if blur:
            Z = ndimage.gaussian_filter(Z, sigma=_options["BLUR_SIGMA"], order=0)

        con = ax_stm.contourf(X, Y, Z, cmap = cmap, levels=levels)

        if display_cbar:
            _fname += "_cbar"
            if display_cbar_on_separate_plot:
                _ax_for_cbar = plt.figure().gca() # TODO: Do we need to set figsize? To make it consistent?
                plt.colorbar(con, ax = _ax_for_cbar)
                _ax_for_cbar.remove()
                colorbar_path = (
                    os.path.join(output_dir, _fname + "_only.png")
                    if output_dir is not None
                    else _fname + "_only.png"
                )
                plt.savefig(
                    colorbar_path,
                    dpi=300,
                    bbox_inches="tight",
                    pad_inches=0,
                )
                plt.close()
                del _ax_for_cbar
            else:
                plt.colorbar(con, ax=ax_stm)

        ax_stm.axis("off")
        ax_stm.set_aspect("equal")

        xmin, xmax = np.min(X), np.max(X)
        ymin, ymax = np.min(Y), np.max(Y)

        ax_stm.set_xlim(xmin, xmax)
        ax_stm.set_ylim(ymin, ymax)

        if display_cell and poscar is not None:
            plot_cell(poscar, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, ls = ":", c = "k", ax=ax_stm,
                      diagonal_transform_for_hexagonal = diagonal_transform_for_hexagonal)
            _fname += "_cell"

        if display_atoms and poscar is not None:
            plot_atoms(poscar, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                       ax = ax_stm, top_n_layers = _options["LAYERS"],
                       diagonal_transform_for_hexagonal = diagonal_transform_for_hexagonal,
                       atoms_info = atoms_info,
                       radii_type=_options["RADIUS_TYPE"],
                       radii_marker_scale=_options["SIZE_RATIO"])
            _fname += "_atoms"

        if not return_subplots:
            image_path = os.path.join(output_dir, f"{_fname}.png")
            plt.savefig(image_path, dpi = 300, bbox_inches='tight', pad_inches=0)
            plt.close()
            generated_files.append(image_path)

    if return_subplots:
        return fig

    if not generated_files:
        raise RuntimeError(
            "No STM image was generated for "
            f"{getattr(current, 'filename', 'CURRENT')} using isosurface range "
            f"[{iso_min}, {iso_max}]."
        )
    return generated_files


def single_apparent_barrier_current(
    current: Current,
    bskan_input: Options,
    image_dir=".",
):
    options = bskan_input.option_dict
    poscar = None
    poscar_path = options["POSCAR"]
    if isinstance(poscar_path, str) and os.path.exists(poscar_path):
        poscar = AR.read_vasp_robust(poscar_path)
    h_min, h_max = lwfplot.current_height_window(current)
    iso_data = _height_values(options, h_min, h_max)
    before_enter_image_dir = _prepare_output_dir(image_dir)
    try:
        for height in iso_data:
            Z_orig = lwfplot.apparent_barrier_from_current(
                current,
                height=float(height),
                fit_radius=float(options["FIT_RADIUS"]),
            )
            X, Y, Z = _repeat_current_grid(current, Z_orig, options, poscar=poscar)
            if options["BLUR_SIGMA"] != 0:
                Z = ndimage.gaussian_filter(Z, sigma=options["BLUR_SIGMA"], order=0)
            nx, ny = options["ITERATION"]
            fname = f"PHI_APP_CURRENT_{float(height):.6g}A_{nx}x{ny}"
            _save_scalar_map(
                X,
                Y,
                Z,
                fname,
                options,
                unit_label=r"$\Phi_{\rm app}$ (eV)",
            )
    finally:
        os.chdir(before_enter_image_dir)


def single_lwf_current(current: Current, bskan_input: Options, image_dir="."):
    """Backward-compatible alias for the former, misnamed CURRENT LWF path."""
    return single_apparent_barrier_current(current, bskan_input, image_dir=image_dir)


def single_vasp_volume(volumetric, bskan_input: Options, volume_path="VASP_VOLUME", image_dir="."):
    options = bskan_input.option_dict
    poscar = None
    poscar_path = options["POSCAR"]
    if isinstance(poscar_path, str) and os.path.exists(poscar_path):
        poscar = AR.read_vasp_robust(poscar_path)
    structure = AR.to_new_cell(poscar) if poscar is not None else None
    top = lwfplot.topmost_z(structure, volumetric)
    c_len = lwfplot.cell_c_length(volumetric)
    simulation = options["SIMULATION"]
    mode = options["IMAGE_MODE"]
    volume_name = os.path.basename(str(volume_path)).upper()
    is_parchg = "PARCHG" in volume_name
    is_locpot = "LOCPOT" in volume_name
    is_chgcar = "CHGCAR" in volume_name

    if simulation == "STM" and mode == "CONSTANT CURRENT":
        iso_min, iso_max = stmplot.get_vasp_density_isorange(volumetric)
        iso_data = _constant_current_iso_values(options, iso_min, iso_max)
    else:
        fit_margin = (
            float(options["FIT_RADIUS"])
            if simulation == "PHI_APP" and is_parchg
            else 0.0
        )
        iso_data = _height_values(options, 0.0, max(c_len - top - fit_margin, 1e-3))

    if simulation == "PHI_APP" and not is_parchg:
        raise ValueError(
            "VASP apparent barrier height (Phi_app) requires PARCHG. "
            "Use LWF for LOCPOT."
        )
    if simulation == "LWF" and not is_locpot:
        raise ValueError("LWF requires LOCPOT and cannot be evaluated from PARCHG/CHGCAR.")
    fermi_level = (
        lwfplot.resolve_fermi_level(volume_path, options.get("FERMI_LEVEL"))
        if simulation == "LWF"
        else None
    )

    before_enter_image_dir = _prepare_output_dir(image_dir)
    try:
        for iso in iso_data:
            if simulation == "STM":
                Z_orig = stmplot.get_surface_vasp(
                    volumetric,
                    float(iso),
                    constant_current=mode == "CONSTANT CURRENT",
                    interpolate="linear"
                    if mode == "CONSTANT HEIGHT" and is_locpot
                    else "exponential",
                )
                if mode == "CONSTANT CURRENT":
                    unit = "Height (A)"
                    prefix = "STM_VASP_CC"
                elif is_parchg:
                    unit = "PARCHG density"
                    prefix = "STM_VASP_CH_PARCHG"
                elif is_locpot:
                    unit = "Potential slice (eV)"
                    prefix = "STM_VASP_CH_LOCPOT"
                elif is_chgcar:
                    unit = "Charge density"
                    prefix = "STM_VASP_CH_CHGCAR"
                else:
                    unit = "Volumetric value"
                    prefix = "STM_VASP_CH"
            elif simulation == "PHI_APP":
                Z_orig = lwfplot.apparent_barrier_from_vasp_density(
                    volumetric,
                    float(iso),
                    topmost=top,
                    fit_radius=float(options["FIT_RADIUS"]),
                )
                unit = r"$\Phi_{\rm app}$ (eV)"
                prefix = "PHI_APP_PARCHG"
            else:
                Z_orig = lwfplot.local_work_function_slice(
                    volumetric,
                    float(iso),
                    fermi_level=fermi_level,
                    topmost=top,
                )
                unit = r"LWF, $\Phi_{\rm loc}$ (eV)"
                prefix = "LWF_LOCPOT"
            X, Y, Z = _repeat_vasp_grid(volumetric, Z_orig, options, poscar=poscar)
            if options["BLUR_SIGMA"] != 0:
                Z = ndimage.gaussian_filter(Z, sigma=options["BLUR_SIGMA"], order=0)
            nx, ny = options["ITERATION"]
            fname = f"{prefix}_{float(iso):.6g}_{nx}x{ny}"
            _save_scalar_map(X, Y, Z, fname, options, unit_label=unit)
    finally:
        os.chdir(before_enter_image_dir)


# Example!
def ex_run_serial_in_python_env(current, bskan_input):
    figs = plt.figure(figsize=(15, 20))
    bskan_input.option_dict["BRIGHTNESS"] = 0.1
    bskan_input.option_dict["ITERATION"] = 3, 2
    bskan_input.option_dict["CONTRAST"] = 0
    bskan_input.option_dict["DISPLAY_CBAR"] = True
    bskan_input.option_dict["DISPLAY_ATOMS"] = True
    bskan_input.option_dict["DISPLAY_CELL"] = True
    bskan_input.option_dict["LAYERS"] = 1
    bskan_input.option_dict["IMAGE_MODE"] = "CONSTANT HEIGHT"
    bskan_input.option_dict["ISO_AUTO"] = False
    bskan_input.option_dict["ISO"] = 4
    bskan_input.option_dict["CONTOUR_RESOLUTION"] = 1000
    single_current(current, bskan_input, image_dir="yaho", fig=figs)
    plt.tight_layout()
    plt.show()

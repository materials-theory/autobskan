##################################################################
## -------------------------About code------------------------- ##
## bSKAN automation ( IMAGE GENERATION )                        ##
## [ Giyeok Lee, Ethan / Inu Kim ] of MTG                       ##
## ------------------------------------------------------------ ##
##################################################################

import os

import matplotlib.pyplot as plt
import numpy as np
from ase.build import make_supercell
from ase.data import covalent_radii, vdw_radii
from ase.geometry.cell import cell_to_cellpar as ctc
from scipy.interpolate import CubicSpline, interp1d

from autobskan.image import AR, vesta_ar, vesta_colors, vesta_ir, vesta_vr

try:
    from ase.data.colors import cpk_colors, jmol_colors
except ImportError:
    from ase.data import cpk_colors, jmol_colors


class Current:
    def __init__(self, filename="CURRENT"):
        with open(filename, 'r') as f:
            ################################### Structural Information
            system_name = f.readline().strip()
            lattice = float(f.readline().strip())
            cell = []
            for i in range(3):
                cell.append(f.readline().split())
            cell = np.array(cell, dtype='d')
            natoms = np.array(f.readline().split(), dtype=int)
            tot_natoms = int(np.sum(natoms))

            # Selective Dynamics & Coord_Type
            what1 = f.readline().strip()
            if what1.upper().startswith("S"):
                selective = True
                what2 = f.readline().strip()
                if what2.upper().startswith("D"):
                    cartesian = False
                else:
                    cartesian = True
            else:
                selective = False
                if what1.upper().startswith("D"):
                    cartesian = False
                else:
                    cartesian = True
            if cartesian:
                coord_type = "CARTESIAN"
            else:
                coord_type = "DIRECT"

            # Coordination
            coord = []
            selective_tags = []
            if not (selective):
                for i in range(tot_natoms):
                    coord.append(f.readline().split())
                coord = np.array(coord, dtype='d')
            else:
                for i in range(tot_natoms):
                    line = f.readline().split()
                    coord.append(line[0:3])
                    selective_tags.append(line[3:])
                coord = np.array(coord, dtype='d')
                selective_tags = np.array(selective_tags, dtype=str)

            ################################### CURRENT Information
            if f.readline().strip() != '':
                raise IOError("something wrong in structural data. The number of positions seems wrong")
            grids = np.array(f.readline().split(), dtype=int)
            nx, ny, nz = grids

            if len(grids) != 3:
                raise IOError("number of grids seems wrong. is it 3D grids? (nx, ny, nz)")

            cur = []
            for i in range(int(np.ceil(nx * ny * nz / 5))):
                cur += f.readline().split()
            #                 cur = np.hstack((cur, np.array(f.readline().split(), dtype = "d"))) <------ takes too long time
            cur = np.array(cur, dtype="d")
            cur_3d = cur.reshape(nz, ny, nx)

            self.system_name = system_name
            self.lattice = lattice
            self.cell = cell
            # self.cellpar = Cell(cell).cellpar()
            self.cellpar = ctc(cell)
            self.natoms = natoms
            self.coord_type = coord_type
            self.coord = coord
            if selective:
                self.selective = selective
                self.selective_tags = selective_tags

            self.grids = grids
            self.cur = cur
            self.cur_3d = cur_3d
            self.filename = os.path.abspath(filename)
            # self.iso_max = np.min(cur_3d[0, :, :])
            # self.iso_min = np.max(cur_3d[nz - 1, :, :])
            self.iso_min = np.max(cur_3d[-1])
            self.iso_max = np.min(cur_3d[0])

            if np.max(np.abs(self.cur_3d)) == 0:
                raise IOError(f"Only 0 values in this current file: {filename}")
            if self.iso_max < self.iso_min:
                raise IOError("This current has a hole. Some grid has a larger minimum value than other grids' maximum")

            self.topmost_atom = np.max(coord @ cell[:,2])

# 2022.02.06 changed (input current obj instead of current.cur_3d matrix)
def get_surface_bskan(current:Current,
                      isosurface,
                      constant_current=True,
                      return_index=False,
                      interpolate="exponential",
                      numeric_precision = 1e-4):
    cur_3d = current.cur_3d
    # cur_3d = current.cur_3d.copy() # Do we need to copy this? This will increase the memory consumption

    _zmin = current.topmost_atom # should be minus value
    _dz = current.cell[2,2]
    if np.any(np.abs(current.cell[2,:2]) > numeric_precision):
        raise IOError("When c lattice vector is not parallel to the cartesian z axis, it is not supported in BSKAN")
        # Even bSKAN doesn't support this case.

    dz = _dz/current.grids[-1] # increment
    if not getattr(current, "_geometry_warning_emitted", False):
        geometry_notes = []
        if np.abs(_zmin + 0.5) > numeric_precision:
            geometry_notes.append(
                f"topmost atom z={_zmin:.9g} A (usual bSKAN STM reference: -0.5 A)"
            )
        if np.abs(dz-0.052918) > numeric_precision:
            geometry_notes.append(
                f"z-grid spacing={dz:.9g} A (usual bSKAN STM spacing: 0.052918 A)"
            )
        if geometry_notes:
            source = os.path.basename(getattr(current, "filename", "CURRENT"))
            print(
                f"[WARN] CURRENT header check for {source}: "
                + "; ".join(geometry_notes)
                + ". Stored header values are used; verify the non-SCF STM tag "
                "if these values are unintended."
            )
        current._geometry_warning_emitted = True

    dz_shift = -_zmin
    # ZGRID = np.linspace(_zmin, _zmin+_dz, current.grids[0])
    # ZGRID -= current.topmost_atom # Make Topmost atom as zero


    if constant_current:
        # Return the matrix of height (Angstrom if real, else index)
        if isosurface >= current.iso_max: # when isosurface >= np.min(cur_3d[0]) # Minimum of maximum isosurfaces
            raise IOError(f"Isosurface should be less than {current.iso_max}")
        elif isosurface <= current.iso_min: # when isosurface <= np.max(cur_3d[-1]) # Maximum of minimum isosurfaces
            # raise IOError(f"Isosurface should be larger than {np.min(cur_3d)}")
            raise IOError(f"Isosurface should be larger than {current.iso_min}")
        else:
            pass

        # i_low = np.argmin(cur_3d > isosurface, axis=0) - 1 # If there are two points, it should return upper one
        i_low = cur_3d.shape[0] - np.argmax(cur_3d[::-1,:,:] > isosurface, axis=0) - 1
        indices = np.indices(i_low.shape)

        if interpolate == "exponential":
            lower = np.log(cur_3d[i_low, indices[0], indices[1]])
            upper = np.log(cur_3d[i_low + 1, indices[0], indices[1]])
            b = lower - upper  # Always posivite
            z_index = (i_low + 1) - (np.log(isosurface) - upper) / b

        elif interpolate == "linear":
            lower = cur_3d[i_low, indices[0], indices[1]]
            upper = cur_3d[i_low + 1, indices[0], indices[1]]
            b = lower - upper
            z_index = (i_low + 1) - (np.abs(upper - isosurface)) / (np.abs(b))

        elif interpolate == "cubic":
            cs = CubicSpline(range(len(cur_3d)), cur_3d, axis=0, extrapolate=False)
            z_index = cs.solve(isosurface)  # But it takes quite a long time..!
            z_index = np.array(z_index)
            z_index = np.array([ele[-1] for ele in z_index.flatten()]).reshape(z_index.shape)

        else:
            raise NotImplementedError(f"{interpolate} method is not implemented. Available: exponential, linear, cubic")

        z_index = np.asarray(z_index, dtype=float)
        z_index = np.clip(z_index, 0, len(cur_3d) - 1)

        if return_index:
            result = z_index
            return result

        else:
    #         result = z_index * 0.0529177 + 0.5  # distance from topmost z (Angstrom)
            result = z_index * dz + dz_shift
            return result

    else:
        # constant height mode: return the matrix of current
    #     if isosurface > len(cur_3d) * 0.0529177 + 0.5:
    #         raise IOError(f"The height should be less than {len(cur_3d) * 0.0529177 + 0.5}")
        if isosurface > _dz + dz_shift:
            raise IOError(f"The height should be less than {_dz + dz_shift} Angstrom")
            # Usually it should be 5.2918 + 0.5 or 5.344718 + 0.5 (depends on the number of grid, nz)
        elif isosurface < dz_shift:
            raise IOError(f"The height should be higher than {dz_shift} Angstrom")
            # Usually it should be 0.5
        else:
            pass

        z_index = (isosurface - dz_shift) / dz
        i_low = np.array(np.floor(z_index), dtype=int)
        i_low = np.clip(i_low, 0, len(cur_3d) - 2)
        t = z_index - i_low

        if interpolate == "exponential":
            lower = np.log(cur_3d[i_low])
            upper = np.log(cur_3d[i_low + 1])
            # Exponential interpolation in log-current domain.
            log_result = (1.0 - t) * lower + t * upper
            result = np.exp(log_result)
            return result

        elif interpolate == "cubic":
            cs = CubicSpline(range(len(cur_3d)), cur_3d, axis=0, extrapolate=False)
            result = cs(z_index)
            return result

        elif interpolate == "linear":
            lin = interp1d([i_low, i_low + 1], (cur_3d[i_low], cur_3d[i_low + 1]), axis=0)
            result = lin(z_index)
            return result
        else:
            raise NotImplementedError(f"{interpolate} method is not implemented. Available: exponential, linear, cubic")


def _vasp_c_length(chgcar):
    cell = AR.to_new_cell_onlycell(np.asarray(chgcar.cell, dtype=float))
    c_length = float(cell[2][2])
    if c_length <= 0:
        raise ValueError("VASP volumetric c lattice length should be positive.")
    return c_length


def _vasp_z_axis(chgcar):
    nz = int(np.asarray(chgcar.pot).shape[0])
    c_length = _vasp_c_length(chgcar)
    return np.arange(nz, dtype=float) * c_length / float(nz)


def _vasp_topmost_atom(chgcar, topmost=None):
    if topmost is not None:
        return float(topmost)
    if getattr(chgcar, "atoms", None) is None:
        return 0.0
    atoms = AR.to_new_cell(chgcar.atoms)
    return float(np.max(atoms.get_positions()[:, 2]))


def get_vasp_density_isorange(chgcar, topmost=None):
    density_3d = np.asarray(chgcar.pot, dtype=float)
    z_axis = _vasp_z_axis(chgcar)
    topmost_atom = _vasp_topmost_atom(chgcar, topmost=topmost)
    above = z_axis >= topmost_atom
    if np.count_nonzero(above) < 2:
        above = np.ones_like(z_axis, dtype=bool)
    density_above = density_3d[above, :, :]
    positive = density_above > 0.0
    if not np.any(positive):
        raise ValueError("VASP density isosurface range needs positive volumetric values.")
    vacuum_floor = float(np.nanmax(np.where(positive[-1], density_above[-1], np.nan)))
    reachable = np.nanmax(np.where(positive, density_above, np.nan), axis=0)
    iso_ceiling = float(np.nanmin(reachable))
    if not np.isfinite(vacuum_floor):
        vacuum_floor = float(np.nanmin(density_above[positive]))
    if not np.isfinite(iso_ceiling) or iso_ceiling <= vacuum_floor:
        positive_values = density_above[positive]
        vacuum_floor = float(np.nanpercentile(positive_values, 1.0))
        iso_ceiling = float(np.nanpercentile(positive_values, 80.0))
    if iso_ceiling <= vacuum_floor:
        raise ValueError("Cannot determine a valid VASP density isosurface range.")
    return vacuum_floor, iso_ceiling


def get_surface_vasp(
    chgcar,
    isosurface,
    constant_current=True,
    return_index=False,
    interpolate="exponential",
    topmost=None,
    numeric_precision=1e-4,
):
    density_3d = np.asarray(chgcar.pot, dtype=float)
    if density_3d.ndim != 3:
        raise IOError("VASP volumetric data should be a 3D grid.")

    if np.any(np.abs(np.asarray(chgcar.cell, dtype=float)[2, :2]) > numeric_precision):
        raise IOError("When c lattice vector is not parallel to the cartesian z axis, it is not supported in VASP")

    topmost_atom = _vasp_topmost_atom(chgcar, topmost=topmost)
    z_axis = _vasp_z_axis(chgcar)
    dz = z_axis[1] - z_axis[0]

    if constant_current:
        iso_min, iso_max = get_vasp_density_isorange(chgcar, topmost=topmost_atom)
        if isosurface >= iso_max:
            raise IOError(f"Isosurface should be less than {iso_max}")
        elif isosurface <= iso_min:
            raise IOError(f"Isosurface should be larger than {iso_min}")
        else:
            pass

        above = z_axis >= topmost_atom
        if np.count_nonzero(above) < 2:
            above = np.ones_like(z_axis, dtype=bool)

        density_above = density_3d[above, :, :]
        z_above = z_axis[above]
        crossing = (density_above[:-1] >= isosurface) & (density_above[1:] < isosurface)
        has_crossing = np.any(crossing, axis=0)
        i_low = crossing.shape[0] - np.argmax(crossing[::-1, :, :], axis=0) - 1
        indices = np.indices(i_low.shape)

        lower = density_above[i_low, indices[0], indices[1]]
        upper = density_above[i_low + 1, indices[0], indices[1]]
        z_lower = z_above[i_low]
        z_upper = z_above[i_low + 1]

        with np.errstate(invalid="ignore", divide="ignore"):
            if interpolate == "exponential":
                log_lower = np.log(lower)
                log_upper = np.log(upper)
                b = log_lower - log_upper
                z_result = z_upper - (np.log(isosurface) - log_upper) / b * (z_upper - z_lower)

            elif interpolate == "linear":
                b = lower - upper
                z_result = z_upper - (np.abs(upper - isosurface)) / (np.abs(b)) * (z_upper - z_lower)

            elif interpolate == "cubic":
                cs = CubicSpline(range(len(density_above)), density_above, axis=0, extrapolate=False)
                z_index = cs.solve(isosurface)
                z_index = np.array(z_index)
                z_index = np.array([ele[-1] for ele in z_index.flatten()]).reshape(i_low.shape)
                z_result = z_above[0] + z_index * dz

            else:
                raise NotImplementedError(f"{interpolate} method is not implemented. Available: exponential, linear, cubic")

        z_result = np.asarray(z_result, dtype=float)
        z_result[~has_crossing] = np.nan

        if return_index:
            result = z_result / dz
            return result

        else:
            result = z_result - topmost_atom
            if not np.any(np.isfinite(result)):
                raise ValueError("No VASP density isosurface crossing was found above the topmost atom.")
            return result

    else:
        z_abs = topmost_atom + float(isosurface)
        z_min = topmost_atom
        z_max = float(z_axis[-1])
        if z_abs > z_max:
            raise IOError(f"The height should be less than {z_max - topmost_atom} Angstrom")
        elif z_abs < z_min:
            raise IOError("The height should be higher than 0.0 Angstrom")
        else:
            pass

        z_abs = min(max(z_abs, 0.0), float(z_axis[-1]))
        z_index = z_abs / dz
        i_low = int(np.floor(z_index))
        i_low = int(np.clip(i_low, 0, len(z_axis) - 2))
        t = z_index - i_low

        if return_index:
            result = z_index
            return result

        if interpolate == "exponential":
            lower = density_3d[i_low]
            upper = density_3d[i_low + 1]
            result = np.full_like(lower, np.nan, dtype=float)
            positive = (lower > 0.0) & (upper > 0.0)
            result[positive] = np.exp(
                (1.0 - t) * np.log(lower[positive]) + t * np.log(upper[positive])
            )
            result[~positive] = (1.0 - t) * lower[~positive] + t * upper[~positive]
            return result

        elif interpolate == "cubic":
            cs = CubicSpline(range(len(density_3d)), density_3d, axis=0, extrapolate=False)
            result = cs(z_index)
            return result

        elif interpolate == "linear":
            lin = interp1d([i_low, i_low + 1], (density_3d[i_low], density_3d[i_low + 1]), axis=0)
            result = lin(z_index)
            return result

        else:
            raise NotImplementedError(f"{interpolate} method is not implemented. Available: exponential, linear, cubic")


get_surface_parchg = get_surface_vasp


def _match_cell_to_image(poscar, xmin, xmax, ymin, ymax, numeric_precision=1e-4):
    h, w = ymax - ymin, xmax - xmin
    lx, ly = np.diag(poscar.cell)[:2]
    times_x, times_y = int(np.ceil(w / lx)), int(np.ceil(h / ly))
    e_xmin = poscar.cell[1, 0] * times_y
    shift_y = 0

    if np.abs(e_xmin) < numeric_precision:
        shift_x = 0
    elif e_xmin > 0:
        shift_x = -int(np.ceil(e_xmin / lx))
        times_x += -shift_x
    else:
        shift_x = 0
        times_x += int(np.ceil(abs(e_xmin) / lx))

    if np.abs(xmin) > numeric_precision:
        shift_x += int(np.floor(xmin / lx))

    if np.abs(ymin) > numeric_precision:
        shift_y += int(np.floor(ymin / ly))

    translation = shift_x * poscar.cell[0] + shift_y * poscar.cell[1]

    return times_x, times_y, translation


def _cut_position_to_image_bounds(
    position,
    xmin,
    xmax,
    ymin,
    ymax,
    numeric_precision=1e-4,
):
    """Return an XY point inside the half-open image bounds, or ``None``."""

    x, y = np.asarray(position, dtype=float)[:2]
    if x < xmin - numeric_precision or y < ymin - numeric_precision:
        return None
    if x > xmax + numeric_precision or y > ymax + numeric_precision:
        return None

    # xmin/ymin belong to the image. xmax/ymax are periodic duplicates of
    # those edges and must be excluded so repeated atoms do not leak in.
    if np.isclose(x, xmax, atol=numeric_precision, rtol=0.0):
        return None
    if np.isclose(y, ymax, atol=numeric_precision, rtol=0.0):
        return None
    if np.isclose(x, xmin, atol=numeric_precision, rtol=0.0):
        x = float(xmin)
    if np.isclose(y, ymin, atol=numeric_precision, rtol=0.0):
        y = float(ymin)
    return float(x), float(y)


def plot_atoms(poscar, xmax, ymax,
               xmin=0, ymin=0,
               ax=None,
               top_n_layers=1,
               al_tol=0.5,
               diagonal_transform_for_hexagonal=True,
               numeric_precision=1e-6,
               atoms_info="ASE-jmol",
               radii_type="default",
               radii_marker_scale=10
               ):
    if atoms_info.strip().upper().startswith("ASE"):
        if radii_type == "default":
            radii_type = "covalent"
        if atoms_info.strip().upper() == "ASE-CPK":
            colors = cpk_colors
            #        [O] vdw_radii / covalent_radii read from ase.data
            #        [X] There are no ionic_radii or atomic_radii in ase.data
            atoms_info = "ASE"
            radii_dict = dict(vdw_radii=vdw_radii, covalent_radii=covalent_radii)
            if radii_type.strip().upper()[0] not in ["C", "V"]:
                raise IOError(f"{radii_type} not provided in ase.data")
        else:
            # atoms_info.strip().upper() == "ASE-JMOL"
            colors = jmol_colors
            #        [O] vdw_radii / covalent_radii read from ase.data
            #        [X] There are no ionic_radii or atomic_radii in ase.data
            atoms_info = "ASE"
            radii_dict = dict(vdw_radii=vdw_radii, covalent_radii=covalent_radii)
            if radii_type.strip().upper()[0] not in ["C", "V"]:
                raise IOError(f"{radii_type} not provided in ase.data")
    elif atoms_info.strip().upper().startswith("V"):
        if radii_type == "default":
            radii_type = "atomic"
        atoms_info = "VESTA"
        colors = vesta_colors
        radii_dict = dict(vdw_radii=vesta_vr, atomic_radii=vesta_ar, ionic_radii=vesta_ir)
    else:
        pass
        # TODO: When the data is provided manually.
    #         atoms_info_path = atoms_info
    #         atoms_info = "VESTA"
    #         colors, vdw_radii, covalent_radii, ionic_radii = read_vesta_ini()

    if radii_type.strip().lower().startswith("i"):
        radii = radii_dict.get("ionic_radii")
    elif radii_type.strip().lower().startswith("c"):
        radii = radii_dict.get("covalent_radii")
    elif radii_type.strip().lower().startswith("v"):
        radii = radii_dict.get("vdw_radii")
    elif radii_type.strip().lower().startswith("a"):
        radii = radii_dict.get("atomic_radii")
    else:
        if atoms_info.startswith("ASE"):
            raise IOError(f"{radii_type} not supported. Choose one between covalent radii or vdw radii for ASE type")
        else:
            raise IOError(
                f"{radii_type} not supported. Choose one among ionic radii, atomic radii, vdw radii for VESTA type")

    if radii is None:
        if atoms_info.startswith("ASE"):
            print("Fallback to ase covalent radii")
            radii = radii_dict.get("covalent_radii")
        elif atoms_info.startswith("V"):
            print("Fallback to VESTA atomic radii")
            radii = radii_dict.get("atomic_radii")

    assert xmax > xmin, f"xmax({xmax}) should be larger than xmin({xmin})."
    assert ymax > ymin, f"ymax({ymax}) should be larger than ymin({ymin})."
    #    poscar = AR.to_new_cell(poscar) # This should be checked!
    if diagonal_transform_for_hexagonal:
        if np.abs(np.diff(poscar.cell.cellpar()[:2])[0]) < numeric_precision:
            gamma = poscar.cell.cellpar()[-1]
            if np.abs(gamma - 60) < numeric_precision or np.abs(gamma - 120) < numeric_precision:
                poscar = AR.to_new_cell(make_supercell(poscar, np.array([[1, -1, 0], [1, 1, 0], [0, 0, 1]])))

    times_x, times_y, translation = _match_cell_to_image(poscar, xmin, xmax, ymin, ymax, numeric_precision)

    surf = AR.Surf(poscar, al_tol = al_tol)
    taglist = sorted(list(range(1, len(surf.dheights) + 1)))
    for i in range(len(taglist)-top_n_layers):
        surf.remove_layer(1)
    surf.sup_xy(times_x + 1, times_y + 1)
    surf.atoms.translate(translation)

    if ax is None:
        ax = plt.figure().gca()

    selected_taglist = sorted(np.unique(surf.atoms.get_tags()))[-top_n_layers:]
    for st in selected_taglist:
        model = surf.atoms[surf.atoms.get_tags() == st]
        for atom in model:
            plot_position = _cut_position_to_image_bounds(
                atom.position,
                xmin=xmin,
                xmax=xmax,
                ymin=ymin,
                ymax=ymax,
                numeric_precision=numeric_precision,
            )
            if plot_position is None:
                continue
            _c = colors[atom.number]
            _r = radii[atom.number]
            ax.plot(*plot_position, c=_c, ms=_r * radii_marker_scale,
                    marker="o")  # TODO: 나중에 원소 별로 한 번에 plot하면 좋긴 할텐데 ㅎ.. 귀찮네
    return ax


def plot_cell(poscar, xmax, ymax, xmin=0, ymin=0, ax=None, diagonal_transform_for_hexagonal=True,
              numeric_precision=1e-6, **kwargs):
    if ax is None:
        ax = plt.figure().gca()
    assert xmax > xmin, f"xmax({xmax}) should be larger than xmin({xmin})."
    assert ymax > ymin, f"ymax({ymax}) should be larger than ymin({ymin})."
    poscar = AR.to_new_cell(poscar)  # This should be checked!
    transformed = False
    if diagonal_transform_for_hexagonal:
        if np.abs(np.diff(poscar.cell.cellpar()[:2])[0]) < numeric_precision:
            gamma = poscar.cell.cellpar()[-1]
            if np.abs(gamma - 60) < numeric_precision or np.abs(gamma - 120) < numeric_precision:
                transformed = True
                poscar = AR.to_new_cell(make_supercell(poscar, np.array([[1, -1, 0], [1, 1, 0], [0, 0, 1]])))

    times_x, times_y, translation = _match_cell_to_image(poscar, xmin, xmax, ymin, ymax, numeric_precision)
    a_vec = poscar.cell[0][:2]
    b_vec = poscar.cell[1][:2]
    trans = translation[:2]

    if transformed:
        for i_nx in range(-times_y, times_x):
            # plot (1,1,0) direction
            origin = i_nx * a_vec + trans
            ax.plot(*np.array([origin, origin + (a_vec + b_vec) * times_y]).T, **kwargs)
        for i_ny in range(1, times_x + times_y):
            # plot (1,-1,0) direction
            origin = i_ny * a_vec + trans
            ax.plot(*np.array([origin, origin + (-a_vec + b_vec) * times_y]).T, **kwargs)

    else:
        for i_nx in range(times_x + 1):
            # plot b lattices
            origin = i_nx * a_vec + trans
            ax.plot(*np.array([origin, origin + b_vec * times_y]).T, **kwargs)
        for i_ny in range(times_y + 1):
            # plot a lattices
            origin = i_ny * b_vec + trans
            ax.plot(*np.array([origin, origin + a_vec * times_x]).T, **kwargs)
    return ax

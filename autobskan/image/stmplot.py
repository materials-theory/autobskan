##################################################################
## -------------------------About code------------------------- ##
## bSKAN automation ( IMAGE GENERATING )                        ##
## [ Giyeok Lee, Ethan / Inu Kim ] of MTG                       ##
## ------------------------------------------------------------ ##
##################################################################

from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import os
import copy
from ase.build import sort
from autobskan.image import AR
from ase.io.vasp import read_vasp, write_vasp
from ase.cell import Cell
from ase.build import make_supercell
from scipy.interpolate import interp1d, CubicSpline


def gy_norm(X, norm_factor=0, min_val=0, max_val=1):
    """
    X : np.array type
    norm_factor : -1 ~ 1, do not recommend outside of this range.
    """
    if norm_factor == 0:
        return (X - X.min()) / (X.max() - X.min()) * (max_val - min_val) + min_val
    else:
        tmp = gy_norm(X, 0, 0, 1)
        power = 10 ** norm_factor
        return gy_norm(tmp ** power, 0, min_val, max_val)


class Surf:
    # adopted from TMCN projects
    def __init__(self, model, al_dir = 2, al_tol = 0.5):
        if isinstance(model, str):
            model = read_vasp(model)

        model = model.copy()  # protect original model
        model.set_constraint(None)
        model.set_tags(None)

        coord_along_axis = model.get_positions()[:, al_dir]
        t = {}
        for z_sorted_index, atom_index in enumerate(np.argsort(coord_along_axis)):
            t[z_sorted_index] = atom_index

        al = 1
        for i in sorted(t):
            if i != 0 and coord_along_axis[t[i]] - coord_along_axis[t[i - 1]] > al_tol:
                al += 1
            model[t[i]].tag = al
        self.atoms  = model
        self.al_dir = al_dir
        self.al_tol = al_tol


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
            self.cell = Cell(cell)
            self.cellpar = Cell(cell).cellpar()
            self.natoms = natoms
            self.coord_type = coord_type
            self.coord = coord
            if selective:
                self.selective = selective
                self.selective_tags = selective_tags

            self.grids = grids
            self.cur = cur
            self.cur_3d = cur_3d
            # self.iso_max = np.min(cur_3d[0, :, :])
            # self.iso_min = np.max(cur_3d[nz - 1, :, :])
            self.iso_min = np.max(cur_3d[-1])
            self.iso_max = np.min(cur_3d[0])


def seek_z_surface(cur_3d, isosurface, real=False, constant_current=True, interpolate="exponential"):
    """
    [input]
    * cur_3d : np.array(shape=(nz, ny, nx)), which is whole current data
    * isosurface : wanted isosurface values
      |___ input this after selecting the possible range of insosurface values.
    * real : if you want the real relative heights from z1 in Angstrom. else, index value will be used
      |___ for the image generations, it doesn't matter.
    * constant_current : Constant current mode if True, else constant height mode will be selected

    [output]
    surface(2d data) of z values
    """
    if constant_current:
        # Return the matrix of height (Angstrom if real, else index)
        if isosurface >= np.min(cur_3d[0]):  # Minimum of maximum isosurfaces
            raise IOError(f"Isosurface should be less than {np.max(cur_3d)}")
        elif isosurface <= np.max(cur_3d[-1]):  # Maximum of minimum isosurfaces
            raise IOError(f"Isosurface should be larger than {np.min(cur_3d)}")
        else:
            pass

        i_low = np.argmin(cur_3d > isosurface, axis=0) - 1
        # In case of multiple occurences of the minimum values when using np.argmin,
        # corresponding to the first occurence are returned. That's why we can just use the argmin here.

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

        if not real:
            return z_index

        else:
            return z_index * 0.0529177 + 0.5  # distance from topmost z (Angstrom)

    else:
        # constant height mode: return the matrix of current
        if isosurface > len(cur_3d) * 0.0529177 + 0.5:
            raise IOError(f"The height should be less than {len(cur_3d) * 0.0529177 + 0.5}")

        elif isosurface < 0.5:
            raise IOError(f"The height should be higher than 0.5 Angstrom")

        else:
            pass

        if not real:
            z_index = isosurface
        else:
            z_index = (isosurface - 0.5) / 0.0529177

        i_low = np.floor(z_index)

        if interpolate == "exponential":
            lower = np.log(cur_3d[i_low])
            upper = np.log(cur_3d[i_low + 1])
            b = lower - upper
            return np.exp(b * z_index - (i_low + 1) + upper)

        elif interpolate == "cubic":
            cs = CubicSpline(range(len(cur_3d)), cur_3d, axis=0, extrapolate=False)
            return cs(z_index)

        elif interpolate == "linear":
            lin = interp1d([i_low, i_low + 1], (cur_3d[i_low], cur_3d[i_low + 1]), axis=0)
            return lin(z_index)

        else:
            raise NotImplementedError(f"{interpolate} method is not implemented. Available: exponential, linear, cubic")


def deprecated_seek_z(cur_z, incur, real=False):
    """
    * cur_z : np.array(shape=(nz,)), current values along Z axis
      |___ input this after selecting x, y coordinate.
           Don't input (nz, ny, nx) shape array

    * incur : wanted isosurface values
      |___ input this after selecting the possible range of insosurface values.

    * real : if you want the real relative height from z1 in Angstrom, we will calculate z_real = z_index * 0.0529177
      |___ for the image generations, it doesn't matter.
      |___ And increments of z is 0.0529177, which has been used in all BSKAN calculations.
      |___ But for now, absolute z value from topmost surface is unclear... (Need to be updated)
          ---> 어디가 시작점이고 어디가 종점인지는 잘 모르겠다. 0.5A부터 시작해서 5.7918 (0.5+5.2918)A 까지인가?
          ---> CURRENT의 z vector가 주어져 있으니, 원자 위치를 대조하여 z방향의 절댓값으로 활용하면 될 것 같다!!
    """
    if cur_z[cur_z == incur].size > 0:
        # if cur_z has incur value, we don't need interpolation
        z_index = int(np.where(cur_z == incur)[0])
    else:
        i_low = len(cur_z[(cur_z - incur) > 0]) - 1  # Since the current increases when z decreases
        i_high = i_low + 1
        b = (np.log(cur_z[i_low]) - np.log(cur_z[i_high]))
        z_index = i_high - (np.log(incur) - np.log(cur_z[i_high])) / b
    if not real:
        return z_index
    else:
        return z_index * 0.0529177


def deprecated_seek_z_surface(cur_3d, incur, real=False):
    """
    [input]
    * cur_3d : np.array(shape=(nz, ny, nx)), which is whole current data
    * incur : wanted isosurface values
      |___ input this after selecting the possible range of insosurface values.
    * real : if you want the real relative heights from z1 in Angstrom. else, index value will be used
      |___ for the image generations, it doesn't matter.
    
    [output]
    surface(2d data) of z values
    """
    ny, nx = cur_3d.shape[1:]
    surface = np.zeros((ny, nx))
    for i in range(ny):
        for j in range(nx):
            surface[i][j] = deprecated_seek_z(cur_3d[:, i, j], incur, real=real)
    return surface


# plot guide atoms
def vasp_to_bskan(str_vasp, current):
    """
    [input]
    * str_vasp = ase.Atoms, vasp structure (POSCAR, ASAMPLE)
    * current = Current instance
    
    [output]
    * str_bskan : ase.Atoms, rectangular structure written in CURRENT file, which is the output file of bSKAN
    * ctype : quadratic / hexagonal / oblique
    """
    if isinstance(str_vasp, str):
        str_vasp = read_vasp(str_vasp)
    poscar_cellpar = str_vasp.get_cell_lengths_and_angles()
    if np.round(poscar_cellpar[-1], 4) in [60., 120.]:
        # sometimes... 119.99999999991849 degree arises
        ctype = "hexagonal"
        A1, A2 = str_vasp.cell[0], str_vasp.cell[1]
        AR1 = AR.length(A1 + A2)
        AR2 = AR.length(A1 - A2)
        if np.round(AR1, 4) == np.round(current.cellpar[0]):
            X = np.array([[1, 1, 0], [1, -1, 0], [0, 0, 1]])
            # AR1 becomes new axis
        else:
            X = np.array([[1, -1, 0], [1, 1, 0], [0, 0, 1]])
            # AR2 becomes new axis        
        str_bskan = make_supercell(str_vasp, X)
    elif poscar_cellpar[-1] == 90:
        ctype = "quadratic"
        str_bskan = AR.to_new_cell(str_vasp.copy())
    else:
        ctype = "oblique"
        mono_l = AR.to_new_cell(str_vasp)
        str_bskan = mono_l.copy()
        x, y, z = mono_l.get_cell()[0, 0], mono_l.get_cell()[1, 1], mono_l.get_cell()[2, 2]
        rectangular_cell = np.diag((x, y, z))
        rectangular_coord = mono_l.get_positions().copy()
        for i in rectangular_coord:
            if i[0] < 0:  # for obtuse angle (angle>90)
                i[0] += x
            if i[0] > x:  # for acute angle (angle<90)
                i[0] -= x
        str_bskan.set_cell(rectangular_cell)
        str_bskan.set_positions(rectangular_coord)
    return str_bskan, ctype


# ----------------------------------------------------------------------------------------------------------------------------
def main(current, bskan_input, image_dir='.', save=True, ax_stm=None, plot_atoms=False):
    """
    * current : Current instance
    * bskan_input : Bskan_input instance (in input.py)
    * image_dir : where to store the images
    """
    if save:
        os.chdir(image_dir)
    else:
        data_to_return = []

    data = tqdm(bskan_input.iso) if isinstance(bskan_input.iso, list) else [bskan_input.iso]
    for iso in data:
        real_x, real_y = current.cellpar[:2]
        try:
            z_surface = seek_z_surface(current.cur_3d, iso)
        except IOError:
            print(f"{iso} is out of range")
            continue
        brightness = bskan_input.brightness
        resolution = bskan_input.contour_resolution
        norm_factor = bskan_input.contrast
        cmap = bskan_input.cmap
        x = np.linspace(0, real_x, current.grids[0])
        y = np.linspace(0, real_y, current.grids[1])
        X, Y = np.meshgrid(x, y)

        if not save:
            data_to_return.append([X, Y, z_surface])

        (brighter, darker) = (0, -brightness) if brightness < 0 else (brightness, 0)

        if ax_stm is None:
            ax_stm = plt.figure(figsize=(real_x, real_y)).gca()
        else:
            # TODO ax_stm의 figure size를 조절해야 할까?
            pass

        ax_stm.contourf(X, Y,
                        gy_norm(z_surface, norm_factor=norm_factor),
                        cmap=cmap,
                        levels=np.linspace(0 - brighter, 1 + darker, int(resolution * (1 + abs(brightness)))))
        ax_stm.axis("off")
        ax_stm.set_aspect("equal")
        if save:
            plt.savefig(f"{iso}.png", dpi=150, bbox_inches='tight', pad_inches=0)
            plt.close()

        if plot_atoms:
            if bskan_input.poscar is not None:
                radius_type = bskan_input.radius_type
                size_ratio = bskan_input.size_ratio
                size_and_color = dict()
                n_layers_for_plot = bskan_input.layers

                module_dir = os.path.dirname(os.path.abspath(__file__))
                with open(os.path.join(module_dir, 'elements_vesta.ini'), 'r') as vesta_data:
                    for line in vesta_data:
                        index, symbol, a_r, v_r, i_r, color_r, color_g, color_b = line.split()
                        if radius_type == "ATOMIC":
                            radius = float(a_r)
                        elif radius_type == "VDW":
                            radius = float(v_r)
                        else:
                            # ionic radius
                            radius = float(i_r)
                        color = tuple(map(lambda x: float(x), (color_r, color_g, color_b)))
                        size_and_color[symbol] = (radius, color)

                if bskan_input.atom_addinfo is not None:
                    # with open(bskan_input.atom_addinfo, 'r') as addinfo:
                    for line in bskan_input.atom_addinfo:
                        symbol, radius, c_r, c_g, c_b = line.split()
                        color_temp = (float(c_r) / 255, float(c_g) / 255, float(c_b) / 255)
                        size_and_color[symbol] = (float(radius), color)

                str_bskan = vasp_to_bskan(bskan_input.poscar, current)[0]
                if save:
                    write_vasp("bskan_structure.vasp", str_bskan, vasp5=True, direct=True, sort=True)
                # atom_species = np.array(AR.species(str_bskan.get_chemical_symbols(), overall=True))
                surf = Surf(str_bskan, al_tol=0.5)
                surf = AR.to_new_cell(sort(surf.atoms, tags=surf.atoms.get_tags()))
                plot_layers = AR.species(surf.get_tags(), overall=True)[-1:-n_layers_for_plot - 1:-1]

                if save:
                    if not os.path.isdir("plot_guide_atoms"):
                        os.mkdir("plot_guide_atoms")
                    os.chdir("plot_guide_atoms")
                    ax_plot_atoms = plt.figure(figsize=(real_x, real_y)).gca()
                    ax_plot_atoms.contourf(X, Y,
                                           gy_norm(z_surface, norm_factor=norm_factor),
                                           cmap=cmap,
                                           levels=np.linspace(0 - brighter, 1 + darker,
                                                              int(resolution * (1 + abs(brightness)))))

                for atom in surf:
                    if atom.tag in plot_layers:
                        size, color_temp = size_and_color[atom.symbol]
                        plot_coord = atom.position
                        if save:
                            ax_plot_atoms.plot(plot_coord[0], plot_coord[1],
                                               color=color_temp, marker=".", ms=size_ratio * size)
                            ax_plot_atoms.axis("off")
                            ax_plot_atoms.set_aspect("equal")
                        else:
                            ax_stm.plot(plot_coord[0], plot_coord[1],
                                        color=color_temp, marker=".", ms=size_ratio * size)

                if save:
                    plt.savefig(f"{iso}_guide.png", dpi=150, bbox_inches='tight', pad_inches=0)
                    os.chdir("../")
                    plt.close()

    if save:
        os.chdir("../")
    else:
        return data_to_return

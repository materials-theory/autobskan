##################################################################
## -------------------------About code------------------------- ##
## bSKAN automation ( IMAGE GENERATING )                        ##
## [ Giyeok Lee, Ethan / Inu Kim ] of MTG                       ##
## ------------------------------------------------------------ ##
##################################################################

import numpy as np
import matplotlib.pyplot as plt
import os
from ase.geometry.cell import cell_to_cellpar as ctc
from ase.build import sort
from  autobskan.image import AR
from ase.io.vasp import read_vasp, write_vasp
from ase.build import make_supercell


def gy_norm(X, norm_factor=0, min_val=0, max_val=1):
    '''
    X : np.array type
    norm_factor : -1 ~ 1, do not recommend outside of this range.
    '''
    if norm_factor==0:
        return (X - X.min()) / (X.max()-X.min()) * (max_val-min_val) + min_val
    else:
        tmp = gy_norm(X, 0, 0, 1)
        power = 10 ** norm_factor
        return gy_norm(tmp**power, 0, min_val, max_val)

class Surf:
    # adopted from TMCN projects
    def __init__(self, model, al_dir = 2, al_tol = 0.5):
        if isinstance(model, str):
            model = read_vasp(model)
        model = model.copy() #protect original model
        model.set_constraint(None)
        model.set_tags(None)

        coord_along_axis = model.get_positions()[:,al_dir]
        t = {}
        for z_sorted_index, atom_index in enumerate(np.argsort(coord_along_axis)):
            t[z_sorted_index] = atom_index
        for i in sorted(t):
            if i == 0:
                al = 1
            else:
                if coord_along_axis[t[i]] - coord_along_axis[t[i-1]] > al_tol:
                    al += 1
            model[t[i]].tag = al
        self.atoms = model
        self.al_dir = al_dir
        self.al_tol = al_tol

class Current:
    def __init__(self, filename = "CURRENT"):
        with open(filename, 'r') as f:
            ################################### Structural Information
            system_name = f.readline().strip()
            lattice = float(f.readline().strip())
            cell = []
            for i in range(3):
                cell.append(f.readline().split())
            cell = np.array(cell, dtype='d')
            natoms = np.array(f.readline().split(), dtype = int)
            tot_natoms = int(np.sum(natoms))
            
            # Selective Dynamics & Coord_Type
            what1=f.readline().strip()
            if what1.upper().startswith("S"):
                selective=True
                what2=f.readline().strip()
                if what2.upper().startswith("D"):
                    cartesian=False
                else:
                    cartesian=True
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
            if not(selective):
                for i in range(tot_natoms):
                    coord.append(f.readline().split())
                coord = np.array(coord, dtype = 'd')
            else:
                for i in range(tot_natoms):
                    line = f.readline().split()
                    coord.append(line[0:3])
                    selective_tags.append(line[3:])
                coord = np.array(coord, dtype = 'd')
                selective_tags = np.array(selective_tags, dtype = str)
                
            ################################### CURRENT Information
            if f.readline().strip() != '':
                raise IOError("something wrong in structural data. The number of positions seems wrong")
            grids = np.array(f.readline().split(), dtype = int)
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
            self.lattice     = lattice
            self.cell        = cell
            self.cellpar     = ctc(cell)
            self.natoms      = natoms
            self.coord_type  = coord_type
            self.coord       = coord
            if selective:
                self.selective      = selective
                self.selective_tags = selective_tags
                
            self.grids       = grids
            self.cur         = cur
            self.cur_3d      = cur_3d
            self.iso_max     = np.min(cur_3d[0, :, :])
            self.iso_min     = np.max(cur_3d[nz - 1, :, :])


def seek_z(cur_z, incur, real = False):
    '''
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
    '''
    if cur_z[cur_z==incur].size > 0 :
        #if cur_z has incur value, we don't need interpolation
        z_index = int(np.where(cur_z == incur)[0])
    else:
        i_low = len(cur_z[(cur_z-incur)>0]) - 1
        i_high = i_low + 1
        b = (np.log(cur_z[i_low])-np.log(cur_z[i_high]))
        z_index = i_high - (np.log(incur)-np.log(cur_z[i_high]))/b
    if not real:
        return z_index
    else:
        return z_index * 0.0529177


def seek_z_surface(current, incur, real = False):
    '''
    [input]
    * cur : np.array(shape=(nz, ny, nx)), which is whole current data
    * incur : wanted isosurface values
      |___ input this after selecting the possible range of insosurface values.
    * real : if you want the real relative heights from z1 in Angstrom. else, index value will be used
      |___ for the image generations, it doesn't matter.
    
    [output]
    surface(2d data) of z values
    '''
    ny, nx = current.shape[1:]
    surface = np.zeros((ny, nx))
    for i in range(ny):
        for j in range(nx):
            surface[i][j] = seek_z(current[:,i,j], incur, real = real)
    return surface

# plot guide atoms
def vasp_to_bskan(str_vasp, current):
    '''
    [input]
    * str_vasp = ase.Atoms, vasp structure (POSCAR, ASAMPLE)
    * current = Current instance
    
    [output]
    * str_bskan : ase.Atoms, rectangular structure written in CURRENT file, which is the output file of bSKAN
    * ctype : quadratic / hexagonal / oblique
    '''
    if isinstance(str_vasp, str):
        str_vasp = read_vasp(str_vasp)
    poscar_cellpar = str_vasp.get_cell_lengths_and_angles()
    if np.round(poscar_cellpar[-1], 4) in [60., 120.]:
        # sometimes... 119.99999999991849 degree arises
        ctype = "hexagonal"
        A1, A2 = str_vasp.cell[0], str_vasp.cell[1]
        AR1 = AR.length(A1+A2)
        AR2 = AR.length(A1-A2)
        if np.round(AR1, 4) == np.round(current.cellpar[0]):
            X = np.array([[1,1,0],[1,-1,0],[0,0,1]])
            # AR1 becomes new axis
        else:
            X = np.array([[1,-1,0],[1,1,0],[0,0,1]])
            # AR2 becomes new axis        
        str_bskan = make_supercell(str_vasp, X)
    elif poscar_cellpar[-1] == 90:
        ctype = "quadratic"
        str_bskan = AR.to_new_cell(str_vasp.copy())
    else:
        ctype = "oblique"
        mono_l = AR.to_new_cell(str_vasp)
        str_bskan = mono_l.copy()
        x, y, z = mono_l.get_cell()[0,0], mono_l.get_cell()[1,1], mono_l.get_cell()[2,2]
        rectangular_cell = np.diag((x,y,z))
        rectangular_coord = mono_l.get_positions().copy()
        for i in rectangular_coord:
            if i[0] < 0: # for obtuse angle (angle>90)
                i[0] += x
            if i[0] > x: # for acute angle (angle<90)
                i[0] -= x
        str_bskan.set_cell(rectangular_cell)
        str_bskan.set_positions(rectangular_coord)
    return str_bskan, ctype


#----------------------------------------------------------------------------------------------------------------------------
def main(current, bskan_input, image_dir='.'):
    '''
    * current : Current instance
    * bskan_input : Bskan_input instance (in input.py)
    * image_dir : where to store the images
    '''
    os.chdir(image_dir)
    if not os.path.exists("plot_guide_atoms"):
        os.makedirs("plot_guide_atoms")
    for iso in bskan_input.iso:
        print(iso)
        real_x, real_y = current.cellpar[:2]
        z_surface      = seek_z_surface(current.cur_3d, iso)
        brightness     = bskan_input.brightness
        resolution     = bskan_input.contour_resolution
        norm_factor    = bskan_input.contrast
        cmap           = bskan_input.cmap
        x = np.linspace(0, real_x, current.grids[0])
        y = np.linspace(0, real_y, current.grids[1])
        X, Y = np.meshgrid(x, y)

        (brighter, darker) = (0, -brightness) if brightness < 0 else (brightness, 0)
        plt.figure(figsize = (real_x, real_y))
        plt.contourf(X, Y, gy_norm(z_surface, norm_factor = norm_factor), cmap = cmap, levels = np.linspace(0-brighter, 1+darker, int(resolution*(1+abs(brightness)))))
        plt.axis("off")
        plt.savefig(f"{iso}.png", dpi=150, bbox_inches='tight', pad_inches=0) #linear의 경우 결과는 똑같다!
        plt.close()

        if bskan_input.poscar != None:
            radius_type    = bskan_input.radius_type
            size_ratio     = bskan_input.size_ratio
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
                    color = tuple(map(lambda x:float(x), (color_r, color_g, color_b)))
                    size_and_color[symbol] = (radius, color)

            if bskan_input.atom_addinfo != None:
                # with open(bskan_input.atom_addinfo, 'r') as addinfo:
                for line in bskan_input.atom_addinfo:
                    symbol, radius, c_r, c_g, c_b = line.split()
                    color_temp = (float(c_r)/255, float(c_g)/255, float(c_b)/255)
                    size_and_color[symbol] = (float(radius), color)

            str_bskan = vasp_to_bskan(bskan_input.poscar, current)[0]
            write_vasp("bskan_structure.vasp", str_bskan, vasp5=True, direct=True, sort=True)
            atom_species = np.array(AR.species(str_bskan.get_chemical_symbols(), overall=True))
            surf = Surf(str_bskan, al_tol = 0.5)
            surf = AR.to_new_cell(sort(surf.atoms, tags = surf.atoms.get_tags()))
            plot_layers = AR.species(surf.get_tags(), overall=True)[-1:-n_layers_for_plot-1:-1]

            os.chdir("plot_guide_atoms")
            plt.figure(figsize = (real_x, real_y))
            plt.contourf(X, Y, gy_norm(z_surface, norm_factor = norm_factor), cmap = cmap, levels = np.linspace(0-brighter, 1+darker, int(resolution*(1+abs(brightness)))))

            for atom in surf:
                if atom.tag in plot_layers:
                    size, color_temp = size_and_color[atom.symbol]
                    plot_coord = atom.position
                    plt.plot(plot_coord[0], plot_coord[1], color=color_temp, marker=".", ms=size_ratio*size)
            plt.axis("off")
            plt.savefig(f"{iso}_guide.png", dpi=150, bbox_inches='tight', pad_inches=0)
            plt.close()
            os.chdir("../")
    os.chdir("../")






from ase.io.vasp import read_vasp, write_vasp
from ase import Atoms
import numpy as np
# import os

class Chgcar:
    def __init__(self, ipf="LOCPOT", direction="Z"):
        if type(ipf) != str:
            if type(ipf) == tuple:
                if len(ipf) == 4:
                    if len(ipf[0]) == 5:
                        return ipf
            raise IOError("input file error! You must input file directory(str type) or already read data(tuple)")
        ip = open(ipf, 'r')
        # print ('* Name of System : ' + ip.readline())
        ip.readline()
        ip.readline()
        cell_vec = []
        cell_vec.append(ip.readline().split())
        cell_vec.append(ip.readline().split())
        cell_vec.append(ip.readline().split())
        cell_vec = np.array(cell_vec, dtype="d")
        symbols = ip.readline().split()
        natoms = list(map(lambda x: int(x), ip.readline().split()))
        tot_natoms = sum(natoms)
        coord_read = ip.readline()
        if coord_read.strip().startswith("D") or coord_read.startswith("d"):
            coord_type = "Direct"
        elif coord_read.strip().startswith("F") or coord_read.startswith("f"):
            coord_type = "Direct"
        elif coord_read.strip().startswith("C") or coord_read.startswith("c"):
            coord_type = "Cartesian"
        else:
            # print("------- Coordination Type Recognition Failed. Calculation ended. -------")
            raise IOError("Coordination type error")

        coord = []
        for i in range(tot_natoms):
            coord.append(list(map(lambda x: float(x), ip.readline().split())))
        ip.readline()
        grids = list(map(lambda x: int(x), ip.readline().split()))
        nx, ny, nz = grids
        # print("* Matrix : [ %d ] x [ %d ] x [ %d ]"%(nx, ny, nz))
        lines = ip.readlines()

        if ipf.endswith("LOCPOT"):
            magnetic = False if np.ceil(len(lines)) == np.ceil(nx * ny * nz / 5) else True
        else:
            # magnetic moment is only provided in LOCPOT case. other cases like CHGCAR, it should be studied
            magnetic = False

        if not magnetic:
            pot = []
            for x in lines:
                if "a" in x:
                    break
                else:
                    tmp = x.split()
                    for y in tmp:
                        pot.append(y)
        else:
            # For now, magnetic system is only supported in LOCPOT
            pot_up, pot_down, temp_ind = [], [], 0
            separate_index = int(np.ceil(nx * ny * nz / 5))
            for x1 in lines[:separate_index]:
                tmp = x1.split()
                for y in tmp:
                    pot_up.append(y)

            while lines[separate_index].split() != [str(nx), str(ny), str(nz)]:
                separate_index += 1

            for x2 in lines[separate_index + 1:]:
                tmp = x2.split()
                for y in tmp:
                    pot_down.append(y)

        if not magnetic:
            # ceil is used just in case nx*ny*nz is not divided into 5
            pot = np.reshape(np.array(pot, dtype="d"), (grids[::-1]))
        else:
            pot_up = np.reshape(np.array(pot_up, dtype="d"), (grids[::-1]))
            pot_down = np.reshape(np.array(pot_down, dtype="d"), (grids[::-1]))
        # if you want to undo reshape, use command, pot.flatten()

        # since chglist is nz ny nx sequence, so if we want nz direction, mean axis will be axis=(1,2)
        # if axis is in tuple type, they calculate average in both direction.
        if direction in "zZcC2":
            direction = "Z"
            ndir = 2
            mean_axis = (1, 2)
        elif direction in "yYbB1":
            direction = "Y"
            ndir = 1
            mean_axis = (2, 0)
        elif direction in "xXaA0":
            direction = "X"
            ndir = 0
            mean_axis = (0, 1)
        else:
            raise IOError("------- Wrong input of direction. Please set direction again -------")

        if not magnetic:
            potavg = np.mean(pot, axis=mean_axis)
        else:
            temp_dict = dict()
            temp_category = ["up", "down"]
            for pot, category in zip([pot_up, pot_down], temp_category):
                potavg = np.mean(pot, axis=mean_axis)
                temp_dict[category] = potavg
        ip.close()

        self.cell = cell_vec
        self.symbols = symbols
        self.natoms = natoms
        self.coord_type = coord_type
        self.coord = coord
        self.grids = grids
        self.magnetic = magnetic
        if not magnetic:
            self.pot = pot
            self.potavg = potavg
        else:
            self.pot_up = pot_up
            self.potavg_up = temp_dict[temp_category[0]]
            self.pot_down = pot_down
            self.potavg_down = temp_dict[temp_category[1]]
            if False not in (self.potavg_down == self.potavg_up):
                self.potavg = temp_dict[temp_category[0]]
            else:
                self.potavg = "WARNING (local potential of spin-up and spin-down is different)"
                # plz check whether potavg_up and potavg_down is same or not.
                # It was always same in my cases.

        try:
            self.atoms = read_vasp(ipf)
        except Exception:
            expanded_symbols = []
            for sym, n in zip(symbols, natoms):
                expanded_symbols.extend([sym] * int(n))

            if len(expanded_symbols) != tot_natoms:
                expanded_symbols = ["X"] * tot_natoms

            coord_array = np.array(coord, dtype=float)[:, :3]
            if coord_type.strip().lower().startswith("d"):
                self.atoms = Atoms(
                    symbols=expanded_symbols,
                    scaled_positions=coord_array,
                    cell=cell_vec,
                    pbc=True,
                )
            else:
                self.atoms = Atoms(
                    symbols=expanded_symbols,
                    positions=coord_array,
                    cell=cell_vec,
                    pbc=True,
                )
        abs_length = self.atoms.cell.cellpar()[ndir]
        if not magnetic:
            self.potavg_x = np.linspace(0, abs_length, len(self.potavg) + 1)[
                            :-1]  # Is this right? As far as I saw in CHGCAR format in VASP, they divide the grids in real-space from 0/max to (max-1)/max
        else:
            self.potavg_x = np.linspace(0, abs_length, len(self.potavg_up) + 1)[:-1]
            # since length of z axis is same, we can use one

    def vac_region(self, vac_width = 1.0, vac_Ediff = 2e-3, loc = "down", return_x = False, return_x_in_index = False):
        '''
        * loc : wanted region to investigate the vacuum level (up / down)
                especially import when vacuum levels of two surfaces are different (i.e. asymmetric slab)
        '''
        if return_x_in_index:
            return_x = True

        planar_data = self.potavg
        x_increments = np.diff(self.potavg_x)[0]
        num_x_index = int(vac_width / x_increments)

        # smoothing is just for determining the vacuum region.
        smoothing = np.zeros(len(planar_data))
        for i in range(-2,3):
            smoothing += np.roll(planar_data, i)
        smoothing /= 5

        vac_bool = abs(np.diff(smoothing)) < vac_Ediff
        vac_bool = np.hstack((vac_bool, np.array([np.abs(planar_data[0]-planar_data[-1])<vac_Ediff])))
        flats_y = {}
        flats_i = {}
        flats_ind = 0
        for i in range(len(vac_bool)):
            if vac_bool[i]:
                if not vac_bool[i-1]:
                    flats_ind += 1
                else:
                    # modified : vac_bool이 True로 시작하는 symmetric slab에서 문제가 있군요!! 0으로 시작하네요
                    if i == 0:
                        flats_ind += 1
                if flats_ind in flats_y:
                    flats_y[flats_ind].append(planar_data[i])
                    flats_i[flats_ind].append(i)
                else:
                    flats_y[flats_ind] = [planar_data[i]]
                    flats_i[flats_ind] = [i]

        long_enough_flats_y = {}
        long_enough_flats_i = {}
        len_ind = 0
        for key in range(1, len(flats_y)+1):
            if len(flats_y[key]) >= num_x_index:
                len_ind += 1
                long_enough_flats_y[len_ind] = flats_y[key]
                long_enough_flats_i[len_ind] = flats_i[key]
        if len(long_enough_flats_y)>2:
            raise ValueError("Too many flat regions are found. Check again")
        elif len(long_enough_flats_y) == 2:
            if loc.strip().upper().startswith("D"):
                vac_region = long_enough_flats_y[1]
                vac_region_i = long_enough_flats_i[1]
            else:
                vac_region = long_enough_flats_y[2]
                vac_region_i = long_enough_flats_i[2]
        elif len(long_enough_flats_y) == 1:
            vac_region = long_enough_flats_y[1]
            vac_region_i = long_enough_flats_i[1]
        else:
            # when no flat regions founded
            raise ValueError("Flat regions which meet the criteria not founded")

        if return_x:
            if return_x_in_index:
                return vac_region_i, vac_region
            else:
                return self.potavg_x[vac_region_i], vac_region
        else:
            return vac_region

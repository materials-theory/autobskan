'''
Package for arranging miscellaneous things
'''
__author__ = "Giyeok Lee"
__email__ = "lgy4230@yonsei.ac.kr"
__date__ = "Sep 14, 2020"
__maintainer__ = "Giyeok Lee"
__version__ = "2023.02.03"
__copyright__ = "Giyeok Lee"

from numpy import cos, sin, tan, sqrt
import numpy as np
import glob, os
from ase import Atoms
from ase.cell import Cell
from ase.io.vasp import read_vasp, write_vasp


def _as_float_triplet(line):
    return [float(x) for x in line.split()[:3]]


def _expand_symbols(symbols, natoms):
    expanded = []
    for sym, n in zip(symbols, natoms):
        expanded.extend([sym] * int(n))
    return expanded


def _atoms_from_parsed_vasp(symbols, natoms, coord_type, coord, cell_vec):
    expanded_symbols = _expand_symbols(symbols, natoms)
    coord_array = np.asarray(coord, dtype=float)[:, :3]
    if len(expanded_symbols) != len(coord_array):
        expanded_symbols = ["X"] * len(coord_array)
    if str(coord_type).strip().lower().startswith("d"):
        return Atoms(
            symbols=expanded_symbols,
            scaled_positions=coord_array,
            cell=np.asarray(cell_vec, dtype=float),
            pbc=True,
        )
    return Atoms(
        symbols=expanded_symbols,
        positions=coord_array,
        cell=np.asarray(cell_vec, dtype=float),
        pbc=True,
    )


def _atoms_from_vasp_header(path):
    with open(path, "r") as fileobj:
        lines = fileobj.readlines()
    if len(lines) < 8:
        raise ValueError("VASP structure file is too short.")

    scale_values = [float(x) for x in lines[1].split()]
    raw_cell = np.asarray([_as_float_triplet(line) for line in lines[2:5]], dtype=float)
    if len(scale_values) == 1:
        scale = scale_values[0]
        if scale < 0.0:
            volume = abs(scale)
            scale = (volume / abs(np.linalg.det(raw_cell))) ** (1.0 / 3.0)
        cell_vec = raw_cell * scale
        coord_scale = scale
    elif len(scale_values) == 3:
        cell_vec = raw_cell * np.asarray(scale_values, dtype=float)[:, None]
        coord_scale = 1.0
    else:
        raise ValueError("VASP scale line should have one or three values.")

    line_index = 5
    first = lines[line_index].split()
    try:
        natoms = [int(x) for x in first]
        symbols = lines[0].split()[: len(natoms)]
        if len(symbols) != len(natoms):
            symbols = ["X"] * len(natoms)
    except ValueError:
        symbols = first
        line_index += 1
        natoms = [int(x) for x in lines[line_index].split()]

    line_index += 1
    coord_type = lines[line_index].strip()
    if coord_type.lower().startswith("s"):
        line_index += 1
        coord_type = lines[line_index].strip()
    line_index += 1

    total = int(sum(natoms))
    coord = [_as_float_triplet(line) for line in lines[line_index: line_index + total]]
    atoms = _atoms_from_parsed_vasp(symbols, natoms, coord_type, coord, cell_vec)
    if not str(coord_type).strip().lower().startswith("d") and coord_scale != 1.0:
        atoms.set_positions(atoms.get_positions() * coord_scale)
    return atoms


def read_vasp_robust(path):
    try:
        return read_vasp(path)
    except Exception:
        return _atoms_from_vasp_header(path)


def _next_nonempty_line(fileobj):
    while True:
        line = fileobj.readline()
        if not line:
            return ""
        if line.strip():
            return line


def _read_vasp_volumetric_dataset(fileobj, value_count):
    """Read one VASP volumetric block without materializing Python tokens."""
    value_count = int(value_count)
    if value_count <= 0:
        raise ValueError("VASP volumetric grid should contain at least one value.")

    block_start = fileobj.tell()
    first_line = _next_nonempty_line(fileobj)
    first_values = np.fromstring(first_line, dtype=float, sep=" ")
    if first_values.size == 0:
        raise ValueError("VASP volumetric data block is missing or invalid.")
    if first_values.size > value_count:
        raise ValueError("VASP volumetric data contains too many values per row.")

    values_per_row = int(first_values.size)
    full_rows, remainder = divmod(value_count, values_per_row)
    fileobj.seek(block_start)

    try:
        if full_rows:
            values = np.loadtxt(
                fileobj,
                dtype=float,
                max_rows=full_rows,
                ndmin=2,
            ).reshape(-1)
        else:
            values = np.empty(0, dtype=float)

        expected_full = full_rows * values_per_row
        if values.size != expected_full:
            raise ValueError(
                "VASP volumetric rows do not have a consistent number of values."
            )

        if remainder:
            tail_line = _next_nonempty_line(fileobj)
            tail = np.fromstring(tail_line, dtype=float, sep=" ")
            if tail.size != remainder:
                raise ValueError(
                    "VASP volumetric final row has an unexpected number of values."
                )
            values = np.concatenate((values, tail))
    except (TypeError, ValueError):
        # Non-standard writers may wrap values at arbitrary columns. NumPy's
        # text reader is slower in this fallback path but remains low-memory.
        fileobj.seek(block_start)
        values = np.fromfile(
            fileobj,
            dtype=float,
            count=value_count,
            sep=" ",
        )

    if values.size != value_count:
        raise ValueError(
            "VASP volumetric data is incomplete: "
            f"expected {value_count} values, found {values.size}."
        )
    return values


def _find_vasp_grid_marker(fileobj, grids):
    expected = tuple(int(value) for value in grids)
    while True:
        line = fileobj.readline()
        if not line:
            return False
        fields = line.split()
        if len(fields) != 3:
            continue
        try:
            candidate = tuple(int(value) for value in fields)
        except ValueError:
            continue
        if candidate == expected:
            return True


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
    def __init__(self, model, al_dir=2, al_tol=0.5):
        '''
        al_dir : Choose among 0(x), 1(y), 2(z) axis
        al_tol : criteria for recognizing as a same atomic layer (default : 0.5AA)
        '''
        if isinstance(model, str):
            model = read_vasp(model)
        model = model.copy()
        model.set_constraint(None)
        model.set_tags(None)

        if isinstance(al_dir, str):
            if al_dir in "0xXaA":
                al_dir = 0
            elif al_dir in "1yYbB":
                al_dir = 1
            else:
                al_dir = 2
        else:
            al_dir = int(al_dir)

        if al_dir not in [0, 1, 2]:
            raise IOError("wrong input direction")

        self.atoms = model
        self.al_dir = al_dir
        self.al_tol = al_tol
        self.height_update()

    def remove_layer(self, layer=-1):
        '''
        remove layer. -1 is the topmost layer
        '''
        if layer == -1:
            layer = np.max(self.atoms.get_tags())

        indexes = []
        for atom in self.atoms:
            if atom.tag == layer:
                indexes.append(atom.index)
        del self.atoms[indexes]
        self.height_update()
        return self.atoms

    def height_update(self):
        model = self.atoms
        coord_along_axis = model.get_positions()[:, self.al_dir]
        interlayer_spacing = {}
        t = dict(enumerate(np.argsort(coord_along_axis)))
        for i in sorted(t):
            if i == 0:
                al = 1
            else:
                if coord_along_axis[t[i]] - coord_along_axis[t[i - 1]] > self.al_tol:
                    al += 1
            model[t[i]].tag = al
            if al not in interlayer_spacing:
                interlayer_spacing[al] = [coord_along_axis[t[i]]]
            else:
                interlayer_spacing[al].append(coord_along_axis[t[i]])
        dheights = []
        for als in interlayer_spacing:
            dheights.append(np.mean(interlayer_spacing[als]))

        self.al_heights = interlayer_spacing
        self.dheights = dheights

    def sup_xy(self, nx=1, ny=1, nz=1, diagonal=False):
        from ase.build import sort
        '''
        [Descriptions]
        : Simply make supercell only along x and y axis.

        [Args]
        nx, ny, nz : int, determines superecell
        '''
        from ase.build import make_supercell
        if not diagonal:
            X = np.diag((nx, ny, nz))
        else:
            X = np.array([[1, -1, 0], [1, 1, 0], [0, 0, 1]])
        self.atoms = sort(make_supercell(self.atoms, X))
        self.height_update()
        return self.atoms


def species(X, overall = True):
    '''
    For comparing just species and its' orders (not number of atoms), we use not np.array, but list.
    [io]
    X : list, tuple, np.array any iterable type.
    overall : don't care about sequences. True if you jlust want to know species
    [output]
    species : list type.
    '''
    species=[]
    for i in X:
        if len(species)==0:
            species.append(i)
        else:
            if overall:
                if i not in species:
                    species.append(i)
            if not overall:
                if i!=species[-1]:
                    species.append(i)
    return species

def npsort(X, axis=-1):
    '''
    Maintaining the coordinates of ions, and sort along axis.
    X = np.array type (n x 3 matrix)
    '''
    temp=[]
    for i in np.argsort(X[:,axis]):
        temp.append(X[i])
    return np.array(temp)
        
def tolsort(X, tol = 1e-5):
    '''
    When X have values with small difference,
    and if you want to remove that values, using this function

    for now, only 1D array is considered
    '''
    before_tol = np.sort(X)
    after_tol = before_tol[[True]+list(np.abs(np.diff(before_tol))>tol)]
    return after_tol        

def hor_translate(X, shift = 0):
    A = X[shift:]
    B = X[:shift]
    return np.hstack((X[shift:], X[:shift]))    

def dir2car(coord, cell_vec):
    '''Direct to Cartesian
    [io] : coord, cell_vec
    |-> coord : coordination. @np.array(dtype="d")
    |-> cell_vec : cell vector. @np.array(dtype='d')
    [output] : new coordination @np.array(dtype='d')
    '''
    return np.matmul(coord, cell_vec)

def car2dir(coord, cell_vec):
    '''Cartesian to Direct
    [io] : coord, cell_vec
    |-> coord : coordination. @np.array(dtype='d')
    |-> cell_vec : cell vector. @np.array(dtype='d')
    [output] : new coordination @np.array(dtype='d')
    '''
    inv_cell_vec = np.linalg.inv(cell_vec)
    return np.matmul(coord, inv_cell_vec)
  
def to_new_cell(model, decimals = 13):
    new_cell = np.zeros((3,3))
    a, b, c, alpha, beta, gamma = model.cell.cellpar()
    alpha, beta, gamma = np.array([alpha, beta, gamma]) * np.pi / 180
    new_cell[0,0] = a
    new_cell[1] = b * cos(gamma), b * sin(gamma), 0
    new_cell[2,0] = c * cos(beta)
    new_cell[2,1] = c * (cos(alpha) - cos(beta) * cos(gamma)) / (sin(gamma))
    new_cell[2,2] = sqrt(c ** 2 - new_cell[2,0] ** 2 - new_cell[2,1] ** 2)
    if decimals != None:
        new_cell = new_cell.round(decimals = decimals)
    tmp = model.copy()
    tmp.set_cell(new_cell)
    tmp.set_scaled_positions(model.get_scaled_positions())
    return tmp

def pop(model, index = [-1]):
    '''
    [Description]
    : Remove selected atoms from model(ase.Atoms instance)
      and return popped atoms with original unit cell.

    [args]
    * index : list of integer

    [Output]
    * return the structure after popping
    '''
    tmp = model.copy()
    removed = []
    for i in model:
        if i.tag in index:
            removed.append(i.index)
    del tmp[removed]
    return tmp

def dirlist():
    all = glob.glob("*")
    direc = [x for x in all if os.path.isdir(x)]
    return direc

def to_new_cell_onlycell(cell):
    if cell.shape != (3, 3):
        raise IOError("Cell input error")

    cell = Cell(cell)
    a, b, c, alpha, beta, gamma = cell.cellpar()
    alpha, beta, gamma = np.array([alpha, beta, gamma]) * np.pi / 180
    new_cell = np.zeros((3, 3))
    new_cell[0, 0] = a
    new_cell[1] = b * cos(gamma), b * sin(gamma), 0
    new_cell[2, 0] = c * cos(beta)
    new_cell[2, 1] = c * (cos(alpha) - cos(beta) * cos(gamma)) / (sin(gamma))
    new_cell[2, 2] = sqrt(c ** 2 - new_cell[2, 0] ** 2 - new_cell[2, 1] ** 2)
    new_cell *= np.abs(new_cell) > 1e-10  # Zero filter
    return Cell(new_cell)


def make_supercell_onlycell(cell, M: np.ndarray):
    if M.shape != (3, 3):
        raise IOError("Supercell matrix input error. The shape should be (3,3)")
    return Cell(M @ cell)


def bskan_cell_transform(cell, numeric_precision=1e-6, diagonal_transform_for_hexagonal=True):
    '''
    Cell transformation performed in bskan
    '''
    if cell.shape != (3, 3):
        raise IOError("Cell input error")
    cell = Cell(cell)
    a, b = cell[:2]
    assert np.abs(np.linalg.norm(a) - a[0]) < numeric_precision,\
        "Your 'a' vector is not in 'x' axis. Perform to_new_cell function before conducting this transform"
    assert np.abs(np.linalg.norm(b) - np.linalg.norm(b[:2])) < numeric_precision,\
        "Your 'b' vector is not in xy plane. Perform to_new_cell function before conducting this transform"
    la, lb = cell.cellpar()[:2]
    gamma  = cell.cellpar()[-1]
    if np.abs(gamma - 90) < numeric_precision:
        return Cell(cell)

    if diagonal_transform_for_hexagonal and np.abs(la - lb) < numeric_precision:
        if np.abs(gamma - 60) < numeric_precision or np.abs(gamma - 120) < numeric_precision:
            M_for_hex = np.array([[1, -1, 0], [1, 1, 0], [0, 0, 1]])
            return to_new_cell_onlycell(make_supercell_onlycell(cell, M_for_hex))
    else:
        return Cell(np.eye(3) * cell)


class Chgcar:
    # from ase.io.vasp import read_vasp, write_vasp
    # import numpy as np
    # import os
    def __init__(self, ipf = "LOCPOT", direction = "Z"):
        if type(ipf) != str:
            if type(ipf) == tuple:
                if len(ipf) == 4:
                    if len(ipf[0]) == 5:
                        return ipf
            raise IOError("input file error! You must input file directory(str type) or already read data(tuple)")
        # Binary buffering keeps file positions reliable after NumPy consumes
        # millions of rows and avoids an additional Unicode decoding layer.
        ip=open(ipf,'rb')
        # print ('* Name of System : ' + ip.readline())
        ip.readline()
        ip.readline()
        cell_vec=[]
        cell_vec.append(ip.readline().split())
        cell_vec.append(ip.readline().split())
        cell_vec.append(ip.readline().split())
        cell_vec=np.array(cell_vec,dtype="d")
        symbols=[value.decode("ascii") for value in ip.readline().split()]
        natoms=list(map(lambda x:int(x),ip.readline().split()))
        tot_natoms=sum(natoms)
        coord_read=ip.readline()
        coord_key = coord_read.strip().lower()
        if coord_key.startswith(b"d"):
            coord_type="Direct"
        elif coord_key.startswith(b"f"):
            coord_type="Direct"
        elif coord_key.startswith(b"c"):
            coord_type="Cartesian"
        else:
            # print("------- Coordination Type Recognition Failed. Calculation ended. -------")
            raise IOError("Coordination type error")

        coord=[]
        for i in range(tot_natoms):
            coord.append(list(map(lambda x:float(x),ip.readline().split())))
        ip.readline()
        grids = list(map(lambda x:int(x),ip.readline().split()))
        nx, ny, nz = grids
        # print("* Matrix : [ %d ] x [ %d ] x [ %d ]"%(nx, ny, nz))
        value_count = nx * ny * nz
        first_dataset = _read_vasp_volumetric_dataset(ip, value_count)
        magnetic = ipf.endswith("LOCPOT") and _find_vasp_grid_marker(ip, grids)

        if not magnetic:
            self._pot = first_dataset.reshape(grids[::-1])
        else:
            # For spin-resolved LOCPOT files, preserve the historical API of
            # exposing the first and second volumetric datasets separately.
            second_dataset = _read_vasp_volumetric_dataset(ip, value_count)
            self._pot_up = first_dataset.reshape(grids[::-1])
            self._pot_down = second_dataset.reshape(grids[::-1])
        # if you want to undo reshape, use command, pot.flatten()

        # since chglist is nz ny nx sequence, so if we want nz direction, mean axis will be axis=(1,2)
        # if axis is in tuple type, they calculate average in both direction.
        if direction in "zZcC2":
            direction="Z"
            ndir=2
            mean_axis=(1,2)
        elif direction in "yYbB1":
            direction="Y"
            ndir=1
            mean_axis=(2,0)
        elif direction in "xXaA0":
            direction="X"
            ndir=0
            mean_axis=(0,1)
        else:
            raise IOError("------- Wrong input of direction. Please set direction again -------")
        self.mean_axis = mean_axis

        ip.close()

        self.cell = cell_vec
        self.symbols = symbols
        self.natoms = natoms
        self.coord_type = coord_type
        self.coord = coord
        self.grids = grids
        self.magnetic = magnetic

        self._calculate_potavg()
        # if non-magnetic -> self.potavg
        # if magnetic -> self.potavg_up, self.potavg_down

        try:
            self.atoms = read_vasp(ipf)
        except Exception:
            self.atoms = _atoms_from_parsed_vasp(
                symbols,
                natoms,
                coord_type,
                coord,
                cell_vec,
            )

        abs_length = self.atoms.cell.cellpar()[ndir]
        if not magnetic:
            self.potavg_x = np.linspace(0, abs_length, len(self.potavg)+1)[:-1] # Is this right? As far as I saw in CHGCAR format in VASP, they divide the grids in real-space from 0/max to (max-1)/max
        else:
            self.potavg_x = np.linspace(0, abs_length, len(self.potavg_up)+1)[:-1]
            # since length of z axis is same, we can use one

    def _calculate_potavg(self):
        if not self.magnetic:
            self.potavg = np.mean(self._pot, axis = self.mean_axis)
        else:
            temp_dict = dict()
            temp_category = ["up", "down"]
            for pot, category  in zip([self._pot_up, self._pot_down], temp_category):
                potavg = np.mean(pot, axis=self.mean_axis)
                temp_dict[category] = potavg
            self.potavg_up = temp_dict[temp_category[0]]
            self.potavg_down = temp_dict[temp_category[1]]

            if np.all(self.potavg_down == self.potavg_up):
                self.potavg = temp_dict[temp_category[0]]
            else:
                self.potavg = None
                # plz check whether potavg_up and potavg_down is same or not.
                # For LOCPOT, it was always same in my cases.

    @property
    def pot(self):
        if self.magnetic:
            return [self._pot_up, self._pot_down]
        return self._pot

    @pot.setter
    def pot(self, new_pot):
        self._pot = np.array(new_pot, dtype=float)
        if self.magnetic:
            assert len(new_pot)==2, "new_pot should contain both pot_up and pot_down"
            self._pot_up = np.array(new_pot[0], dtype=float)
            self._pot_down = np.array(new_pot[1], dtype=float)
        self._calculate_potavg()

    @pot.deleter
    def pot(self):
        del self._pot
        del self.potavg
        if self.magnetic:
            del self._pot_up
            del self._pot_down
            del self.potavg_up
            del self.potavg_down

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

    def write_supercell(self, filename = "CHGCAR_x_y_z", x=1, y=1, z=1):
        '''
        Now testing... Plz be careful when using, and use this only for reference
        filename : the output CHGCAR name. default : f"CHGCAR_{x}_{y}_{z}"
        x, y, z  : determining how big is the supercell (https://www.c2x.org.uk/downloads/c2x_man.html 에서 CHGCAR 다루는 거 보고 따라하자)
        '''
        filename = f"CHGCAR_{str(x)}_{str(y)}_{str(z)}" if filename == "CHGCAR_x_y_z" else filename
        supercell_atoms = make_supercell(self.atoms, np.diag((x, y, z)))
        supercell_grids = np.array([x, y, z]) * self.grids

        chg_x = copy.deepcopy(self.pot)
        for iterate_x in range(x-1):
            chg_x = np.vstack((chg_x, self.pot))

        chg_y = copy.deepcopy(chg_x)
        for iterate_y in range(y-1):
            chg_y = np.hstack((chg_y, chg_x))

        chg_z = copy.deepcopy(chg_y)
        for iterate_z in range(z-1):
            chg_z = np.dstack((chg_z, chg_y))

        write_vasp(f"POSCAR_{x}_{y}_{z}", supercell_atoms, vasp5=True, direct=True)
        self.write(filename, atoms=supercell_atoms, pot=chg_z, grids=supercell_grids)


    def write(self, filename = "CHGCAR_gen", atoms = None, pot = None, grids = None):
        if pot is None:
            pot = self.pot
        pot = pot.flatten()
        if atoms is None:
            atoms = self.atoms
        if grids is None:
            grids = self.grids

        write_vasp(filename, atoms, vasp5=True, direct=True)
        with open(filename, "a") as chgcar:
            print("", file=chgcar)
            print(f"{grids[0]} {grids[1]} {grids[2]}", file=chgcar)
            for i, value in enumerate(pot):
                # Let't think how to optimize this part for reducing time.
                endtag = " " if (i+1)%5 != 0 else "\n"
                print(f"{value:17.11E}", end = endtag, file = chgcar)
        # How to write 'augmentation occupancies' which is in CHGCAR?
        # Above one, and other additional things are needed for actual VASP non-scf calculation...

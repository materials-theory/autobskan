# First version written Giyeok Lee, 2018_08_21
# Modified by Giyeok Lee, 2018_10_07 (Using argparse), 2019_01_10 (Adding Monoclinic_iteration)
# Modified by Giyeok Lee, 2023_02_03 (Revising array_iter -> to general gamma value)

import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from PIL import Image
from ase.cell import Cell
from autobskan.image import AR

from scipy.interpolate import griddata
from scipy.ndimage import zoom


def bskan_cell_transform(cell, numeric_precision=1e-6, diagonal_transform_for_hexagonal=True):
    '''
    Cell transformation performed in bskan
    (From VASP to BSKAN)
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
        return cell

    if diagonal_transform_for_hexagonal and np.abs(la - lb) < numeric_precision:
        if np.abs(gamma - 60) < numeric_precision or np.abs(gamma - 120) < numeric_precision:
            M_for_hex = np.array([[1, -1, 0], [1, 1, 0], [0, 0, 1]])
            return AR.to_new_cell_onlycell(AR.make_supercell_onlycell(cell, M_for_hex))
        else:
            # TODO: a==b 이면서, 60도나 120도가 아닌 경우
            raise RuntimeError("아직 안됐어용")
    else:
        return Cell(np.eye(3) * cell)

def array_stack(Z_orig, nx = 2, ny = 2):
    # _l = Z_orig.copy()
    # __l = np.hstack([_l]*nx)
    # ___l = np.vstack([__l]*ny)
    # return ___l
    return np.tile(Z_orig, (ny, nx))

def array_iter(Z_orig, nx=2, ny=2, cell=None,
               gamma=90, l_a: float = None, l_b: float = None,
               diagonal_transform_for_hexagonal=True, numeric_precision=1e-6,
               input_filetype: str = "BSKAN",  # 1) BSKAN CURRENT format: "BSKAN" / 2) VASP CHGCAR format: "VASP"
               accurate_cut_by_interpolation: bool = True,
               filled_array_for_contourf=True, verbose=True
               ):
    '''
    Iterate the array

    :param Z_orig: Array to repeat. (2D matrix)
    :param nx: Iterate along cartesian x axis (lattice vector, a should be aligned with x axis)
    :param ny: Iterate along cartesian y axis (y component of lattice vector, b)

    * Choose one between (a) Put cell itself "Cell" / (b) manual input: "Gamma, l_a, l_b"
    (a) Input cell
        (Here, when input_filetyle == "BSKAN", you need to input structure file's cell info, not the cell info from Current)
    :param cell: np.ndarray with shape of (3,3) or ase.cell.Cell obj

    (b) Manual input of cell params.
        (Note that l_a, l_b should be orthogonalized lattice vector from BSKAN CURRENT file)
    :param gamma: [From POSCAR file] Angle of a and b lattice. Default = 90 Degree (Rectangular in this case)
    :param l_a:   [From BSKAN CURRENT file] Length of a lattice vector, not needed when input_filetype == "VASP"
    :param l_b:   [From BSKAN CURRENT file] Length of b lattice vector, not needed when input_filetype == "VASP"


    :param diagonal_transform_for_hexagonal: If (gamma in [60, 120]) and (l_a==l_b), perform diagonal transform, and set iterate without cutting
    :param numeric_precision: Numerical precision when comparing l_a==l_b, and angle
    :param input_filetype: ["BSKAN" / "VASP"], case insensitive 1) BSKAN CURRENT format / 2) VASP CHGCAR format
    :param accurate_cut_by_interpolation: Since the array data is always discrete (grids), we lose lots of data when cut the images if gamma!=90.
                                          So we perform cubic interpolation, and then cut the image, and iterate it.
    :param accurate_cut_by_interpolation:
    :return:
    '''
    Z_orig = np.array(Z_orig, dtype=float)
    # ------------------------------------------------------- input_filetype == "VASP"
    if input_filetype.upper().startswith("V"):
        if diagonal_transform_for_hexagonal:
            # Just want ot notice that diagonal_transform is meaningless for VASP CHGCAR format
            if verbose:
                print("[INFO] VASP CHGCAR format doesn't perform any diagonal transform."\
                      "Automatically fallback to diagonal_transform_for_hexagonal=False")
            diagonal_transform_for_hexagonal = False
        if accurate_cut_by_interpolation:
            # Just want ot notice that accurate_cut is meaningless for VASP CHGCAR format
            if verbose:
                print(
                    "[INFO] VASP CHGCAR format doesn't require accurate_cut_by_interpolation."\
                    "Automatically fallback to accurate_cut_by_interpolation=False")
            accurate_cut_by_interpolation = False
        if cell is None:
            if verbose:
                print("[Warning] When input_filetype==VASP, This iter requires cell info."\
                      "The output will be stacked the arrays!")
            if (l_a is not None) | (l_b is not None):
                if verbose:
                    print("[INFO] l_a and l_b is not needed when input_filetype==VASP. l_a and l_b will be ignored.")
            if filled_array_for_contourf:
                Z_new = array_stack(Z_orig, nx=nx, ny=ny)
                Z_new = np.hstack([Z_new, Z_new[:, 0].reshape(-1, 1)])
                Z_new = np.vstack([Z_new, Z_new[0].reshape(1, -1)])
                return Z_new
            else:
                return array_stack(Z_orig, nx=nx, ny=ny)
        else:
            assert cell.shape == (3, 3), "[Error] Cell should be either np.ndarray or ase.cell.Cell obj," \
                                         "which has the shape of (3,3)"

        a = cell[0].copy()
        b = cell[1].copy()
        assert np.abs(np.linalg.norm(a) - a[0]) < numeric_precision
        assert np.abs(np.linalg.norm(b) - np.linalg.norm(b[:2])) < numeric_precision

        #         x = np.linspace(0, a[0]*nx, Z_orig.shape[1]*nx+1)[:-1]
        #         y = np.linspace(0, b[1]*ny, Z_orig.shape[0]*ny+1)[:-1]
        #         b_x_shift = np.linspace(0, b[0]*ny, Z_orig.shape[0]*ny+1)[:-1]
        #         X, Y = np.meshgrid(x, y)
        #         X += b_x_shift.reshape(-1, 1)
        #         return X, Y, array_stack(Z_orig, nx = nx, ny = ny)

        if filled_array_for_contourf:
            x = np.linspace(0, a[0] * nx, Z_orig.shape[1] * nx + 1)
            y = np.linspace(0, b[1] * ny, Z_orig.shape[0] * ny + 1)
            b_x_shift = np.linspace(0, b[0] * ny, Z_orig.shape[0] * ny + 1)
        else:
            x = np.linspace(0, a[0] * nx, Z_orig.shape[1] * nx + 1)[:-1]
            y = np.linspace(0, b[1] * ny, Z_orig.shape[0] * ny + 1)[:-1]
            b_x_shift = np.linspace(0, b[0] * ny, Z_orig.shape[0] * ny + 1)[:-1]

        X, Y = np.meshgrid(x, y)
        X += b_x_shift.reshape(-1, 1)
        Z_new = array_stack(Z_orig, nx=nx, ny=ny)
        if not filled_array_for_contourf:
            return X, Y, Z_new

        # When filled_array_for_contourf
        Z_new = np.hstack([Z_new, Z_new[:, 0].reshape(-1, 1)])
        Z_new = np.vstack([Z_new, Z_new[0].reshape(1, -1)])
        return X, Y, Z_new

    else:
        assert input_filetype.upper().startswith(
            "B"), f"[Error] Incompatible filetype: {input_filetype}. Choose between BSKAN or VASP"
    # ------------------------------------------------------- input_filetype == "BSKAN"
    if cell is None:
        if np.abs(gamma - 90) > numeric_precision:
            if l_a is None or l_b is None:
                l_b, l_a = Z_orig.shape
                if verbose:
                    print("[Warning] When gamma != 90, I recommend to put either [cell] or [l_a & l_b]."\
                          "Here, 'l_b, l_a = Z_orig.shape' will be used")
                if (np.abs(gamma - 60) < numeric_precision) or (np.abs(gamma - 120) < numeric_precision):
                    # Without cell, without l_a, l_b... But gamma in [60, 120]...
                    if diagonal_transform_for_hexagonal:
                        raise IOError("[Error] When (gamma in [60, 120]) and diagonal_transform,"\
                                      "l_a and l_b should be given")
                    else:
                        # Go to "cutting part"
                        pass
        else:
            # when gamma == 90, same as VASP format in here :D
            return array_iter(Z_orig[:-1, :-1], nx=nx, ny=ny, cell=None, input_filetype="VASP",
                              filled_array_for_contourf=filled_array_for_contourf, verbose=False)
    #             Z_new = array_stack(Z, nx=nx, ny=ny)
    #             if not filled_array_for_contourf:
    #                 return Z_new
    #             else:
    #                 Z_new = np.hstack([Z_new, Z_new[:,0].reshape(-1,1)])
    #                 Z_new = np.vstack([Z_new, Z_new[0].reshape(1,-1)])
    #                 return Z_new
    else:
        if np.abs(gamma - 90) > numeric_precision or l_a is not None or l_b is not None:
            print("[Warning] When cell is given, gamma, l_a, l_b will be ignored.")

        gamma = Cell(cell).cellpar()[-1]
        #         l_a, l_b = np.linalg.norm(cell, axis=1)[:2] # Not this!! The cell should be orthogonalized to match BSKAN CURRENT
        if np.abs(gamma - 90) < numeric_precision:
            # cell is given, gamma == 90.
            return array_iter(Z_orig[:-1, :-1], nx=nx, ny=ny, cell=cell,
                              input_filetype="VASP", filled_array_for_contourf=filled_array_for_contourf,
                              numeric_precision=numeric_precision, verbose=False)

        bskan_cell = bskan_cell_transform(cell, diagonal_transform_for_hexagonal=diagonal_transform_for_hexagonal, numeric_precision = numeric_precision)
        l_a, l_b = bskan_cell.cellpar()[:2]
        l_a_poscar, l_b_poscar = cell.cellpar()[:2]
        if diagonal_transform_for_hexagonal and (np.abs(l_a_poscar - l_b_poscar) < numeric_precision):
            if (np.abs(gamma - 60) < numeric_precision) or (np.abs(gamma - 120) < numeric_precision):
                return array_iter(Z_orig, nx=nx, ny=ny, cell=bskan_cell,
                                  filled_array_for_contourf=filled_array_for_contourf,
                                  numeric_precision=numeric_precision,
                                  accurate_cut_by_interpolation=accurate_cut_by_interpolation,
                                  input_filetype="BSKAN", verbose=verbose)

    ## Cutting part
    # input_filetype == "BSKAN" / gamma!=90 / No diagonal_transform
    if accurate_cut_by_interpolation & ((cell is not None) | ((l_a is not None) & (l_b is not None))):
        # When l_a, l_b == Z_orig.shape, this is not needed
        Z_orig = zoom(Z_orig, 3)
    Z = Z_orig[:-1, :-1]

    size_y, size_x = Z.shape

    if (cell is not None) | (l_a is not None and l_b is not None):
        _dgrid = l_a / size_x
        x_cut = int(np.round(cell[1][0] / _dgrid, 0))
    else:
        # This is not that accurate. Reason: Each grids' height and width is not 100% identical.
        if verbose:
            print(
                "[Warning] We recommend to put cellinfo inside. The height & Width cannot be 100% identical in every cases.")
        _G = np.radians(gamma)
        _s_g, _c_g = np.sin(_G), np.cos(_G)
        x_cut = int(np.round(size_y * _c_g / _s_g, 0))

    Z_result = np.hstack([Z] * nx)
    Z_hor = Z_result.copy()
    for i in range(1, ny):
        if x_cut != 0:
            attach = np.roll(Z_hor, axis=1, shift=x_cut * i)
        else:
            attach = Z_hor.copy()
        Z_result = np.vstack([Z_result, attach])

    if not filled_array_for_contourf:
        X = np.linspace(0, l_a * nx, size_x * nx + 1)[:-1]
        Y = np.linspace(0, l_b * ny, size_y * ny + 1)[:-1]
        X, Y = np.meshgrid(X, Y)
        return X, Y, Z_result

    Z_result = np.hstack([Z_result, Z_result[:, 0].reshape(-1, 1)])
    top = list(np.roll(Z_hor, axis=1, shift=x_cut * ny)[0])
    top.append(top[0])
    top = np.array(top, dtype=float).reshape(1, -1)
    Z_result = np.vstack([Z_result, top])

    if l_a is not None:
        assert l_b is not None
        X = np.linspace(0, l_a * nx, size_x * nx + 1)
        Y = np.linspace(0, l_b * ny, size_y * ny + 1)
        X, Y = np.meshgrid(X, Y)
        return X, Y, Z_result
    else:
        return Z_result

# -------------------------------------------------------------------------------- From imagefile. not from array itself
# [image iterating] -> generate png and save it to generated_iter/ directory
def image_iter(file, x, y, gamma, name, savedir="generated_iter"):
    im = Image.open(file)
    im = im.resize(np.array(np.round(np.array(im.size)/max(x, y)), dtype=int)) # somewhat numerical error due to integer pixels..
    size_x, size_y = im.size
    if np.round(gamma, 4) in [60., 120.]:
        gamma = 90
    G=np.radians(gamma)
    s_g, c_g = np.sin(G), np.cos(G)

    # when gamma=90, x_cut=0
    x_cut = int(round(size_y*c_g/s_g, 0))
    new_image = Image.new("RGBA", (size_x * x, size_y * y), (255, 255, 255, 0))

    if x_cut==0:
        for ix in range(x):
            for iy in range(y):
                new_image.paste(im, (size_x * ix, size_y * iy, size_x * (ix + 1), size_y * (iy + 1)))
    else:
        for iy in range(y):
            ix_cut = (iy*x_cut) % size_x
            im_left  = im.crop(((ix_cut, 0, size_x, size_y)))
            im_right = im.crop(((0, 0, ix_cut, size_y)))
            new_image.paste(im_left, (0, size_y * iy, size_x-ix_cut, size_y * (iy+1)))
            new_image.paste(im_right, (size_x*x-ix_cut, size_y*iy, size_x*x, size_y * (iy+1)))
            for ix in range(x-1):
                new_image.paste(im, ((size_x-ix_cut)+size_x*ix, size_y*iy, (size_x-ix_cut)+size_x*(ix+1), size_y*(iy+1)))
    if not os.path.exists(savedir):
        os.makedirs(savedir)
    new_image.save(f"{savedir}/{name}_{x}x{y}.png", "PNG")
    return f"{savedir}/{name}_{x}x{y}.png"

# [blurring] -> generate png and save it to generated_blur/ directory
def image_blur(file, blur_sigma, name, savedir = "generated_blur"):
    '''blur_sigma = blur strength, only blur task -> onlyblur='t' '''
    imgb=plt.imread(file)
    imgb2 = ndimage.gaussian_filter(imgb, sigma=(blur_sigma, blur_sigma, 0), order=0)
    plt.imshow(imgb2, interpolation='nearest')
    if not os.path.exists(savedir):
        os.makedirs(savedir)
    plt.imsave(f"{savedir}/{name}_{blur_sigma}.png", imgb2)

# --------------------------------------------------------------------------------
def main(raw_image_dir, bskan_input):
    '''
    [for next step of image generations]
    * raw_image_dir = directory path
    * bskan_input = Bskan_input instance of bskan automation code.
    '''
    if hasattr(bskan_input, "option_dict"):
        options = bskan_input.option_dict
        ext = str(options["EXT"]).lower()
        iteration = options["ITERATION"]
        gamma = options["GAMMA"]
        blur_sigma = options["BLUR_SIGMA"]
    else:
        ext = bskan_input.ext
        iteration = bskan_input.iteration
        gamma = bskan_input.gamma
        blur_sigma = bskan_input.blur_sigma

    for i in ["generated_iter", "generated_blur"]:
        if not os.path.exists(i):
            os.makedirs(i)
    raw_files = glob.glob(f"{raw_image_dir}/*.{ext}")
    for raw_file in raw_files:
        name = raw_file.split("/")[-1].replace("." + ext, "")
        iterated = image_iter(raw_file, iteration[0], iteration[1], gamma, name, savedir = "generated_iter")
        image_blur(iterated, blur_sigma, iterated.split("/")[-1].replace(f".{ext}", ""), savedir = "generated_blur")

# # -------------------------------------------------------------------------------- [operate] : run functions for all files in current directory
# if __name__ == "__main__":
#     import argparse

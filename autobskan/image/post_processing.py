# Giyeok Lee, 2018_08_21
# Modified by Giyeok Lee, 2018_10_07 (Using argparse), 2019_01_10 (Adding Monoclinic_iteration)
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sin, cos
import os
import glob
from PIL import Image

# -------------------------------------------------------------------------------- [image iterating] -> generate png and save it to generated_iter/ directory
def image_iter(file, x, y, gamma, name, savedir="generated_iter"):
    im = Image.open(file)
    im = im.resize(np.array(np.array(im.size)/max(x, y), dtype=int)) # somewhat numerical error due to integer pixels..
    size_x, size_y = im.size
    if np.round(gamma, 4) in [60., 120.]:
        gamma = 90
    G=gamma*pi/180
    s_g, c_g = sin(G), cos(G)

    # when gamma=90, x_cut=0
    x_cut = int(round(size_y*c_g/s_g,0))
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


# --------------------------------------------------------------------------------  [blurring] -> generate png and save it to generated_blur/ directory
def image_blur(file, blur_sigma, name, savedir = "generated_blur"):
    '''blur_sigma = blur strength, only blur task -> onlyblur='t' '''
    imgb=plt.imread(file)
    imgb2 = ndimage.gaussian_filter(imgb, sigma=(blur_sigma, blur_sigma, 0), order=0)
    plt.imshow(imgb2, interpolation='nearest')
    if not os.path.exists(savedir):
        os.makedirs(savedir)
    plt.imsave(f"{savedir}/{name}_{blur_sigma}.png", imgb2)

def main(raw_image_dir, bskan_input):
    '''
    [for next step of image generations]
    * raw_image_dir = directory path
    * bskan_input = Bskan_input instance of bskan automation code.
    '''
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
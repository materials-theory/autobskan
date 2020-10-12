# remove black edge of tiff image(simulated stm image from openDX) and make images(png : no black edge) iterated and blurring
# Giyeok Lee, 2018_08_21
# Modified by Giyeok Lee, 2018_10_07 (Using argparse), 2019_01_10 (Adding Monoclinic_iteration)
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt
from math import pi, sin, cos
import os
import glob
from PIL import Image
import argparse
import numpy as np
from ase.io.vasp import read_vasp

def rm_blackbox(file,name="cut_image"):
    im=Image.open(file)
    size_x,size_y=im.size
    directory = os.path.dirname("rm_blackbox/")
    if not os.path.exists(directory):
        os.makedirs(directory)
    
    pix=im.load()
    
    def blackbox_check(a,b):
        while True:
            if pix[a,b]==(0,0,0):
                break
            else:
                if a==size_x-1:
                    a=0; b+=1
                else:
                    a+=1

        #find the x size of black edge box
        for i in range(1,int(size_x/2)):
            h_a=size_x-i
            if pix[h_a,b]==(0,0,0):
                break
        #find the y size of black edge box
        for p in range(1,int(size_y/2)):
            v_b=size_y-p
            if pix[a,v_b]==(0,0,0):
                break

        # box test
        if pix[h_a,v_b]==(0,0,0):
            return (a,b,h_a,v_b,"ok")
        else:
            return (a,b,h_a,v_b,"fail")

    (a,b,h_a,v_b,c)=blackbox_check(0,0)

    if c=="ok":
        #im.crop((left,upper,right,lower))
        img2=im.crop(((a+1,b+1,h_a,v_b)))
        img2.save("rm_blackbox/"+name+".png","PNG")
    else:
        print("rm_blackbox failed in \"%s\". Can't find black edge."%(file))
    
def image_iter(file, x, y, gamma, name, savedir="generated_iter"):
    im = Image.open(file)
    im = im.resize(np.array(np.round(np.array(im.size)/max(x, y)), dtype=int)) # somewhat numerical error due to integer pixels..
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

def image_blur(file, blur_sigma, name, savedir = "generated_blur"):
    blur_sigma = np.float(blur_sigma)
    imgb=plt.imread(file)
    imgb2 = ndimage.gaussian_filter(imgb, sigma=(blur_sigma, blur_sigma, 0), order=0)
    plt.imshow(imgb2, interpolation='nearest')
    if not os.path.exists(savedir):
        os.makedirs(savedir)
    plt.imsave(f"{savedir}/{name}_{blur_sigma}.png", imgb2)

def main():
    pars = argparse.ArgumentParser()
    pars.add_argument('-n', type = int,nargs=2,default=[4, 2], help='number of iteration. nx ny')
    pars.add_argument('-gamma', type=str, default=90.0, help="gamma value or the name of structure")
    pars.add_argument('-m', type = int, help='Which task do you want. 1=Do all tasks, 2=iteration+blur, 3=remove black edge+iteration, 4=only remove black edge, 5=only iterating, 6=only blur', default=2)
    pars.add_argument('-bm',type=int,help='blur mode. (how to select blur sigma) 1:single blur sigma value(bsv), 2:range of blur sigma value(bsl)',choices=[1,2], default=1)
    pars.add_argument('-bsv',type=float, default='5',help='single blur sigma value (when you choose blur mode=1)')
    pars.add_argument('-bsl',type=str,nargs=2,default=['20', '25'],help="range of blur sigma value. (when you choose blur mode=2)")
    pars.add_argument('-ext',type=str,default='png',help="format of input images (ex. png, tiff, jpg)")
    args = pars.parse_args()

    mode, blur_mode, blur_sigma, blur_sigma_r, ext, gamma = args.m, args.bm, args.bsv, args.bsl, args.ext, str(args.gamma)
    if blur_mode==2:
        blur_sigma=" ".join(blur_sigma_r)
    x,y=int(args.n[0]),int(args.n[1])

    if os.path.exists(gamma):
        gamma = read_vasp(gamma).cell.cellpar()[-1]
    else:
        try:
            gamma = np.float(gamma)
        except:
            raise IOError("wrong input of gamma value. put filename of structure or exact value of gamma")

    files=glob.glob("*."+ext)
    print("Run for %s file(s)."%(len(files)))    
    for file in files:
        name=file.replace(f".{ext}", "")
        pre_path = ""
        # mode 1 = remove black edge(box) + Image iteration + Image blurring
        if mode in [1, 3, 4]:
            rm_blackbox(file, name)
            pre_path = "rm_blackbox/"
        if mode in [1, 2, 3, 5]:
            try:
                image_iter(f"{pre_path}{name}.png", x, y, gamma, name)
            except FileNotFoundError as e:
                print(e)
                image_iter(file, x, y, gamma, name)
            pre_path = "generated_iter/"
        if mode in [1, 2, 6]:
            try:
                image_blur(f"{pre_path}{name}_{x}x{y}.png", blur_sigma, name)
            except FileNotFoundError as e:
                print(e)
                image_blur(file,blur_sigma, name)
        # exist just because english ordinal
        if files.index(file)+1 in [1,2,3]:
            ind=['st','nd','rd'][files.index(file)]
        else:
            ind='th'
        print("%2s%s file done among %s files" %(files.index(file)+1,ind,len(files)))
    print("job's done")

if __name__ == "__main__":
    main()
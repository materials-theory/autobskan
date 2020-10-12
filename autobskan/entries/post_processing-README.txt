* It uses python3

[Functions]
There are 3 functions in this script.
(1) rm_blackbox : extract images inside black line box.
(2) image_iter : iterate extracted images, iteration number is from -n argument.
(3) image_blur : blurring images. 

[main arguments]
-n : number of iteration. If you want to get 4x1 images, add "-n 4 1"
-m : type of mode. 1=Run all functions(defulat), 2=iteration+blur, 3=only extraction, 4=only iteration, 5=only blur
-ext : input file's extension. Default is 'tiff'.

[Operating]
Runs on all files that are in current path(directory) and have that extension (by -ext argument.)


[Example]
$python3 gyb_1007.py
= $python3 gyb_1007.py -n 4 2 -m 1 -ext tiff


$python3 gyb_1007.py -m 3
 > only extraction from tiff files
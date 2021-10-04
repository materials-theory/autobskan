

=========================
(1) Installation
=========================

Using pip command to download autobskan ::

  $ pip install autobskan


Or download via github, and unzip it. ::

  $ pip install -r requirements.txt

  $ python setup.py install






==================================================
(2) Explanation of bskan.in input file.
==================================================

``MODE`` = CALCULATION / **IMAGE** / POST_PROCESSING

  [tip] you can also type ca/im/po




------------------------------------------------
a. Calculation Parts (Not supported yet)
------------------------------------------------

``VASP`` = command of executing VASP binary

``BSKAN`` = command of excuting BSKAN binary

  (ex. mpirun -np 8 ~/programs/VASP/norm_std.x > stdout.log)

``METHOD`` = TH (Tersoff Hamann) / CHEN / BARDEEN (Numerical)

  [tip] you can also write th(or te) / ch / ba(or nu)

``BIAS`` = start_end_steps & value

  1) If you put separator "_", it is regarded as MIN_MAX_INCREMENTS to make list of input bias values.
       For example, -0.01_0.03_0.02 will be regarded as [-0.01, 0.01, 0.03].
       
  2) If you put separator "&", it is regarded as ingredients of lists.
       For example, 0.02 & 0.05 will be regarded as [0.02, 0.05].
       
  3) Or you can write the combination of 1) and 2).
       For example, "-0.05_0.01_0.02 & 0.00" will give you array([-0.05, -0.03, -0.01, 0., 0.01])
       
  +) Space will be ignored. Do not worry about this.
 

------------------------
b. Image Parts
------------------------

``CURRENT`` = filename of CURRENT. [Default = CURRENT]

  * For filenames, multiple choices are also possible. (Must be divided by "&")
  * Or, following Regular expressions. (ex. '\*current' indicates all files which have 'current' words in the last.)


``ISO_AUTO`` = True / False / **LOGSCALE** [Default = LOGSCALE]

  [tip] you can also type t/f/l
  
``ISO`` = value(s) / numbers of wanted isosurface values.

  1) When ISO_AUTO = True, ISO will be the number of isosurface values. [Default = 5]
                     * ex) generate 5 images when ISO = 5
                     
  2) When ISO_AUTO = False, ISO will be the exact value of isosurface. (Same format of BIAS input)
                     * ex1) 2500_1e4_2500 will be regarded as [2500., 5000., 7500., 10000.]
                     * ex2) 1e3 & 1e4 & 1e5 will be regarded as [1000., 10000., 100000.]
                     
  3) When ISO_AUTO = LOGSCALE, ISO will be possible 10^x value
                    * here, x will be set automatically from minimum to maximum. (2020.08.29 updated)
                    * +) input ISO value will be ignored

``CMAP`` = name of colormap [Default = afmhot]

  Colomaps are following matplotlib.pyplot cmaps

``CONTRAST`` = 0

  * normalization factor (float). 0 is default, which is linear normalization.

  * value from -1 to 1 is recommended, and usually absolte value within 0.2 is enough.
  
``BRIGHTNESS``  = -1 ~ 0[Default] ~ 1

``CONTOUR_RESOLUTION`` = 200

``EXT`` = Wanted extension type or raw_images. [Default = png]

``POSCAR`` = filename of structure file. vasp5 POSCAR format is supported.

``ATOM_ADDINFO`` = If you want manual setting of Atomic size and colors, you can put the filename with atomic informations.

  * Default setting is identical to default setting of VESTA program
  
  * For example, if you want to change the size and color of hydrogen atom, to 1.5 Angstrom and white color,
  
  * You can set ATOM_ADDINFO = element.txt in bskan.in which includes information of H by following commands,
  
  * $echo 'H 1.5 255 255 255' > elements.txt
  
``LAYERS`` = number of layers to plot from surface. [Default = 1]

``RADIUS_TYPE`` = **atomic** / van der Waals / ionic radius [Default = ATOMIC]

  [tip] you can also type a / v / i
  
``SIZE_RATIO`` = marker size ratio of atoms [Default = 30]


------------------------
c. Postprocessing Parts
------------------------

``POST_PROCESSING`` = Precede to the post_processing process or not. [Default = FALSE]

``ITERATION`` = nx, ny [Default = 4, 4] # iterations along x / y axis

``BLUR_METHOD`` = GAUSSIAN [Default] # For now, there is only one choice (Gaussian).

``BLUR_SIGMA`` = Postive number [Default = 10]

``GAMMA`` = Manual input gamma value of lattice parameter. Using in iteration process. [Default = 90]

  Only when there is no POSCAR file. If POSCAR file exists, it will automatically calculate this value.




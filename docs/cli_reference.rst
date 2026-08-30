CLI reference
=============

Commands
--------

``autobskan --input PATH``
  Read an AutoBSKAN input file (default ``bskan.in``) and run its ``TASK``.

``autobskan pre-check POSCAR``
  Classify bSKAN symmetry translations before SCF. Detailed stages are visible
  by default; ``--quiet`` prints only the decision.

``autobskan-gui``
  Start the interactive local application. Available flags are ``--host``,
  ``--port``, ``--no-browser``, ``--debug``, and ``--version``.

Input syntax
------------

``bskan.in`` uses one ``KEY = VALUE`` pair per line. Text after ``#`` is a
comment. Numeric lists use commas. A three-part underscore expression denotes
``start_end_step``; for example ``-0.5_0.5_0.5`` yields ``-0.5, 0.0, 0.5``.

Task and calculation options
----------------------------

``TASK``
  ``CALCULATION``, ``IMAGE``, or ``POST_PROCESSING_ONLY``. Default ``IMAGE``.

``VASP_COMMAND`` / ``BSKAN_COMMAND``
  Shell commands for the external executables.

``SCF_PATH``
  Existing VASP SCF directory. AutoBSKAN validates its POSCAR and CHGCAR before
  preparing NONSCF/bSKAN files.

``TIP_PATH``
  Tip directory. Chen expects ``ATIP``, ``WAVTIP``, and ``PROCARtip``.

``METHOD``
  ``TH``, ``CHEN``, or ``BARDEEN``.

``BIAS``
  Bias value, comma-separated values, or ``start_end_step`` expression.

For ``METHOD=CHEN``, AutoBSKAN constructs each ``CURRENT`` file directly from
the calculation's ``ASAMPLE``, ``CURSAVE``, and ``INSCAN`` files. The generated
header follows the bSKAN grid and spacing conventions, and the body retains the
established scientific numeric format. A separate Tersoff-Hamann calculation is
not required solely to supply a CURRENT header.

Image options
-------------

``SIMULATION``
  ``STM``, ``PHI_APP``, or ``LWF``. The ``LWF`` input keyword denotes the
  :math:`\Phi_{\rm loc}` LOCPOT workflow.

``INPUT_SOURCE``
  ``BSKAN`` or ``VASP``.

``CURRENT`` / ``VOLUME``
  ``CURRENT`` path or glob for bSKAN; ``PARCHG``, ``LOCPOT``, or ``CHGCAR`` path
  for VASP. ``VOLUME`` falls back to ``CURRENT`` when omitted.

``FERMI_LEVEL``
  Fermi level in eV for ``SIMULATION = LWF``. Empty, ``AUTO``, or ``NONE`` reads
  a sibling ``OUTCAR``.

``IMAGE_MODE``
  ``CONSTANT CURRENT`` or ``CONSTANT HEIGHT``.

``ISO_AUTO`` / ``ISO``
  ``ISO_AUTO=TRUE`` is the default logarithmic mode and is equivalent to
  ``LOGSCALE``; it uses powers of ten strictly inside each file's valid
  isosurface range. ``ISO_AUTO=LINEAR`` uses ``ISO`` evenly spaced current
  values inside that range. With ``FALSE``, ``ISO`` contains explicit values.
  If no complete decade is available, use ``LINEAR`` or provide an explicit
  ``ISO``.

``FIT_RADIUS``
  :math:`\Phi_{\rm app}` fit half-width in A. Default ``0.5``.

``CMAP`` / ``CONTRAST`` / ``BRIGHTNESS``
  Matplotlib colormap and display-normalisation settings.

``CONTOUR_RESOLUTION`` / ``EXT``
  Contour level count and output extension for non-interactive image output.

Post-processing and overlays
----------------------------

``POSCAR``
  Optional structure path used for iteration and atom/cell overlays.

``ITERATION``
  Repetition counts along the first and second lattice vectors. Default ``1,1``.

``GAMMA``
  Lattice gamma in degrees when no structure file is available.

``BLUR_METHOD`` / ``BLUR_SIGMA``
  Gaussian filtering and sigma in grid pixels. Sigma zero disables filtering.

``DISPLAY_ATOMS`` / ``LAYERS``
  Show atoms and select the number of outermost layers.

``ATOMS_INFO`` / ``ATOM_ADDINFO`` / ``RADIUS_TYPE`` / ``SIZE_RATIO``
  Atom style, optional element data, radius basis, and marker scale.

``DISPLAY_CELL`` / ``DISPLAY_CBAR`` / ``DISPLAY_CBAR_ON_SEPARATE_PLOT``
  Cell and colour-scale output switches.

``DIAGONAL_TRANSFORM_FOR_HEXAGONAL``
  Enable the established bSKAN-compatible diagonal transform for equal-length
  60/120-degree in-plane cells.

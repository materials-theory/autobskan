Scientific scope and assumptions
================================

STM surface maps
----------------

For a bSKAN ``CURRENT`` file, constant-current mode locates the uppermost
crossing of the selected current at each XY point. Constant-height mode
interpolates between adjacent z planes. VASP ``PARCHG`` and ``CHGCAR`` use the
same surface-extraction interface; positive density values are interpolated in
the logarithmic domain where exponential interpolation is selected.

The stored VASP grid follows the conventional ``(nz, ny, nx)`` array order.
Physical coordinates span grid points ``0`` through ``(n-1)/n`` of each lattice
period. Repeated boundaries are half-open so atoms at ``xmax`` or ``ymax`` are
not drawn a second time.

Apparent barrier height
-----------------------

AutoBSKAN reports

.. math::

   \Phi_{\rm app} = 0.952495\left(\frac{\partial \ln Q}{\partial z}\right)^2

in eV under the implemented Angstrom-based decay convention. ``Q`` is bSKAN
``CURRENT`` or VASP ``PARCHG``. The derivative is obtained by least-squares
fitting ``ln(Q)`` in a window centred at the selected height; non-positive
samples are excluded and at least three valid z points are required.

:math:`\Phi_{\rm app}` is an apparent tunnelling barrier. It depends on the
decay observable, fit window, tip/sample state, and numerical setup and is not
interchangeable with an electrostatic local work function.

Finite-height local work function
---------------------------------

For VASP ``LOCPOT`` AutoBSKAN evaluates

.. math::

   \Phi_{\rm loc}(x,y;z_0) = V(x,y,z_0)-E_{\rm F}

using linear z interpolation. ``z0`` is measured above the topmost atom in the
active embedded or explicitly confirmed structure. ``E_F`` is user supplied or
read from the last matching line in a sibling ``OUTCAR``.

This is a finite-height local potential reference, not automatically the
asymptotic macroscopic work function. Its interpretation depends on vacuum
thickness, dipole corrections, slab asymmetry, electrostatic gauge, and height.
For an electrostatic map, generate ``LOCPOT`` with ``LVHAR = .TRUE.`` and report
the chosen height and Fermi-level source.

Geometry domain
---------------

VASP surface extraction and the bSKAN symmetry pre-check currently require the
``c`` lattice vector to be parallel to the Cartesian surface normal. The
pre-check explicitly returns ``UNSUPPORTED_GEOMETRY`` when ``c`` is not
perpendicular to the ``ab`` plane.

The symmetry pre-check reconstructs bSKAN-compatible lattice and basis
operations, including fractional translations. A proposed origin shift is
applied with ``ase.Atoms.translate()``, wrapped, and checked again. It cannot
repair intrinsic glide or screw translations; these return ``UNSAFE``.

Filtering and display resampling
--------------------------------

Gaussian filtering is applied on the scalar array with periodic ``wrap``
boundaries. Coordinates and unit-cell bounds are unchanged. The filtered field
is used consistently for map colours, line profiles, FFT, and export.

GUI display quality can lower only the raster sent to Plotly. The cached
analysis field retains the full post-interpolation grid, so changing display
quality must not alter physical XY resolution or exported numerical profiles.

Reproducibility record
----------------------

A publication should record at least the input source and filename, calculation
method, bias, simulation quantity, acquisition mode, isosurface or height,
fit half-width for :math:`\Phi_{\rm app}`, Fermi level and LOCPOT settings for
:math:`\Phi_{\rm loc}`, repetition, interpolation factor, filter sigma, and
software version.

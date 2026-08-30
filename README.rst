AutoBSKAN
=========

.. image:: https://github.com/materials-theory/autobskan/actions/workflows/ci.yml/badge.svg
   :target: https://github.com/materials-theory/autobskan/actions/workflows/ci.yml
   :alt: Continuous integration status

.. image:: https://img.shields.io/badge/license-MIT-green.svg
   :target: LICENSE
   :alt: MIT license

AutoBSKAN is a Python workflow and interactive analysis application for
simulated scanning tunnelling microscopy (STM). It connects VASP and bSKAN
calculation output to reproducible surface maps, apparent-barrier analysis,
finite-height local work-function analysis, overlays, line profiles, Fourier
analysis, and publication-oriented export.

Statement of need
-----------------

STM simulation workflows often combine electronic-structure output, a separate
tunnelling-current code, ad hoc conversion scripts, and plotting notebooks.
That fragmentation makes it difficult to reproduce grid conventions, height
references, symmetry handling, and figure settings. AutoBSKAN provides one
tested interface for preparing bSKAN calculations and analysing native bSKAN
``CURRENT`` files or VASP volumetric files without changing their physical XY
sampling during interactive display optimisation.

The package does not distribute VASP or bSKAN. Calculation automation requires
separately installed and licensed executables; analysis of existing files does
not.

Capabilities
------------

* bSKAN ``CURRENT`` and VASP ``PARCHG``, ``LOCPOT``, and ``CHGCAR`` input.
* Constant-current and constant-height STM maps.
* Apparent barrier height, :math:`\Phi_{\rm app}`, from the fitted logarithmic
  decay of ``CURRENT`` or ``PARCHG``.
* Finite-height local work function, :math:`\Phi_{\rm loc}(x,y;z_0)`, from
  ``LOCPOT`` referenced to the Fermi level.
* Symmetry pre-check and origin-shift proposal before an SCF calculation.
* Interactive atom and unit-cell overlays, line profiles, FFT, Gaussian
  filtering, manual colour limits, and independent high-resolution export.
* Chunked browser upload for large volumetric files and cached full-resolution
  analysis after the initial parse.

Scientific quantities
---------------------

AutoBSKAN keeps two quantities that are sometimes both called a local work
function separate:

.. math::

   \Phi_{\rm app} = 0.952495\left(\partial_z \ln Q\right)^2,

where ``Q`` is the bSKAN current or VASP partial charge density under the
implemented decay convention, and

.. math::

   \Phi_{\rm loc}(x,y;z_0) = V_{\rm LOCPOT}(x,y,z_0) - E_{\rm F}.

``PARCHG`` decay is therefore reported as :math:`\Phi_{\rm app}`, never as
:math:`\Phi_{\rm loc}`. The latter requires ``LOCPOT`` and an explicit Fermi
level or a sibling ``OUTCAR``. See `Scientific scope`_ for assumptions and
limitations.

Installation
------------

AutoBSKAN requires Python 3.10 or newer. Install the published package with::

  python -m pip install autobskan

For an editable source installation::

  git clone https://github.com/materials-theory/autobskan.git
  cd autobskan
  python -m pip install -e ".[test]"

Static image export uses ``plotly>=6.1,<6.8`` and ``kaleido>=1.0,<2``. This
compatibility range avoids the Plotly 6.8 header API, which requires Kaleido
1.3 or newer. Kaleido 1 also needs a local Chrome or Chromium installation;
``plotly_get_chrome`` installs a compatible browser when required.

Interactive GUI
---------------

Start the local application with::

  autobskan-gui

The version badge identifies the running package. Select the data source,
simulation, and volumetric file in the **Data** and **Model** tabs. VASP files
provide their embedded structure automatically; an explicitly selected
POSCAR/CONTCAR can replace it after confirmation.

Slider changes render live, while adjacent numeric fields allow exact values.
The fast and balanced display modes reduce only the number of browser display
pixels. Line analysis and export continue to use the full rendered scalar
field. Manual ``vmin``/``vmax`` applies to the map, colour scale, and line
profile; returning to automatic range re-enables brightness.

Useful launch options are::

  autobskan-gui --host 127.0.0.1 --port 8050
  autobskan-gui --no-browser
  autobskan-gui --debug
  autobskan-gui --version

The server is intended for local use. Closing the final GUI window stops the
server. Do not expose it on an untrusted network: the local-path field is
designed to read files available to the launching user.

Pre-check before SCF
--------------------

Run the visible four-stage bSKAN symmetry check before starting SCF::

  autobskan pre-check POSCAR

The command reports one of four outcomes:

``READY``
  No origin shift is needed. ``ASAMPLE`` is written in VASP 4 style, with
  Selective dynamics and Direct coordinates.

``SHIFT_REQUIRED``
  ``POSCAR-autobskan_shifted.vasp`` is written using
  ``ase.Atoms.translate()`` and validated again. Start SCF from this structure;
  do not reuse the previous ``WAVECAR``, ``CHGCAR``, or ``WAVSAMPLE``.

``UNSAFE``
  No single global origin shift removes all detected translations. Use a
  symmetry-free/full-Brillouin-zone workflow or correct the bSKAN symmetry
  expansion. No output structure is written.

``UNSUPPORTED_GEOMETRY``
  The ``c`` lattice vector is not perpendicular to the ``ab`` plane, so the
  current pre-check cannot guarantee symmetry safety. No output is written.

Use ``--quiet`` for the final decision only. Calculation automation repeats the
check against an existing SCF and refuses post-SCF shifting.

Calculation automation
----------------------

After a ``READY`` pre-check, configure the external executables and calculation
directories in ``bskan.in``::

  TASK = CALCULATION
  VASP_COMMAND = mpirun -np 8 /path/to/vasp_std
  BSKAN_COMMAND = mpirun -np 8 /path/to/bskan
  SCF_PATH = /absolute/path/to/scf
  TIP_PATH = /absolute/path/to/tip
  METHOD = CHEN
  BIAS = -0.5_0.5_0.5

``METHOD`` accepts ``TH``, ``CHEN``, or ``BARDEEN``. ``BIAS`` accepts one
value, comma-separated values, or ``start_end_step`` notation. AutoBSKAN
prepares the non-SCF and bSKAN directories, tracks completed stages, and writes
resulting ``CURRENT`` files under ``3_result``.

For Chen calculations, the ``CURRENT`` header and body are generated directly
from ``ASAMPLE``, ``CURSAVE``, and ``INSCAN`` using bSKAN-compatible grid,
spacing, and numeric formatting. A preliminary Tersoff-Hamann run is therefore
not required solely to provide the header. Use absolute ``SCF_PATH`` and
``TIP_PATH`` values so the workflow remains independent of the launch directory.

Image CLI
---------

The non-interactive image path reads ``bskan.in``::

  TASK = IMAGE
  SIMULATION = STM
  INPUT_SOURCE = BSKAN
  CURRENT = CURRENT
  IMAGE_MODE = CONSTANT CURRENT
  ISO_AUTO = FALSE
  ISO = 1e-5

Run it with::

  autobskan --input bskan.in

For batch STM rendering, ``ISO_AUTO = TRUE`` is equivalent to ``LOGSCALE`` and
uses valid powers of ten. ``ISO_AUTO = LINEAR`` instead uses ``ISO`` evenly
spaced current values per CURRENT file. GUI exports write source-specific
``CURRENT`` or ``VOLUME`` entries and accessible structure paths as absolute
paths.

For :math:`\Phi_{\rm app}`, set ``SIMULATION = PHI_APP`` and optionally set
``FIT_RADIUS`` (default ``0.5`` A). For :math:`\Phi_{\rm loc}`, use
``SIMULATION = LWF``, ``INPUT_SOURCE = VASP``, ``VOLUME = LOCPOT``, and either
``FERMI_LEVEL`` or a sibling ``OUTCAR``. The complete option list is in the
`CLI reference`_.

Documentation
-------------

* `Documentation index`_
* `GUI guide`_
* `CLI reference`_
* `Scientific scope`_
* `Python API`_

.. _Documentation index: docs/index.rst
.. _GUI guide: docs/gui.rst
.. _CLI reference: docs/cli_reference.rst
.. _Scientific scope: docs/scientific_scope.rst
.. _Python API: docs/api.rst

Testing and reproducibility
---------------------------

Run the complete test suite with::

  python -m pytest

The tests cover VASP volumetric streaming, bSKAN ``CURRENT`` parsing and
conversion, symmetry pre-check outcomes, monoclinic plotting, GUI callback
contracts, cache and browser lifecycle, line-profile export, and image-export
configuration. Scientific changes should include a small redistributable
fixture and a numerical regression against the original code or an independent
reference.

Contributing, support, and citation
-----------------------------------

Development and validation requirements are in `CONTRIBUTING.rst`_. Report
reproducible problems through the `issue tracker`_. Security-sensitive reports
should follow `SECURITY.md`_. Cite the version used according to
``CITATION.cff``.

.. _CONTRIBUTING.rst: CONTRIBUTING.rst
.. _issue tracker: https://github.com/materials-theory/autobskan/issues
.. _SECURITY.md: SECURITY.md

License
-------

AutoBSKAN is distributed under the MIT License. See ``LICENSE``.

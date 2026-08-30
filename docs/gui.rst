GUI guide
=========

Launch and lifecycle
--------------------

Run ``autobskan-gui`` from the directory containing the files to analyse. The
default address is ``http://127.0.0.1:8050``. A second **New window** session has
independent browser state. Closing the final window stops the local server.

Use ``--no-browser`` in a headless environment. ``--debug`` exposes detailed
callback and parse status and should not be used for routine analysis.

Data
----

Choose the source explicitly:

* **bSKAN** reads ``CURRENT``.
* **VASP** reads ``PARCHG``, ``LOCPOT``, or ``CHGCAR``.

An accessible local path avoids copying the file. Browser selection transfers a
large volumetric file in 8 MB chunks and displays transfer progress. Uploaded
copies live in a process-owned temporary directory and disappear when the GUI
process exits.

A VASP volumetric file supplies its embedded atomic structure automatically.
Loading a separate POSCAR or CONTCAR opens a confirmation dialog before that
structure replaces the embedded one for overlays and cell-dependent processing.

Model
-----

**STM** supports constant-current and constant-height maps. For bSKAN
constant-current data, the isosurface control is ``log10(current)``. A VASP
``LOCPOT`` selected as STM is restricted to constant-height potential slices.

:math:`\Phi_{\rm app}` uses a constant-height fit centre and a configurable fit
half-width. The complete interval must remain within the available vacuum grid.

:math:`\Phi_{\rm loc}` accepts only VASP ``LOCPOT``. Enter ``E_F`` in eV or
leave it empty to read the last Fermi-energy record from an ``OUTCAR`` beside
the LOCPOT.

Display and analysis
--------------------

Sliders update while dragged and every precision slider has a numeric field.
Controls that cannot affect the active image are disabled and do not trigger a
render.

**Fast** and **Balanced** resample only the browser display raster. They do not
change the physical XY bounds, the full cached scalar field, line-profile
sampling, FFT input, or export data. **Full** sends the full raster to the
browser and can be slower for large repeated grids.

Gaussian filtering uses periodic boundary conditions and leaves the coordinate
axes unchanged. When enabled, the filtered scalar field is used consistently by
the displayed map, colour limits, line profile, FFT, and export.

Select P1 by clicking once on the map; a second click adds P2 and the line.
**Snap to extrema** moves each endpoint to the nearest local extremum. The CSV
export contains ``x`` (distance in A) and ``y`` (the plotted value).

Colour range and export
-----------------------

Automatic range follows the active scalar field and permits brightness
adjustment. Manual ``vmin`` and ``vmax`` disable brightness and apply the exact
same y limits to the line profile. **Auto range** restores automatic values.

Image export is independent of the viewport. Choose PNG, JPEG, WebP, SVG, or
PDF and set an output dimension up to 12,000 px, subject to a 60-megapixel
limit. Width and height remain locked to the physical XY aspect ratio. Raster
surface exports have a transparent background where the format supports it.
The separate colourbar export is an opaque, margin-free strip without ticks or
labels.

Failures are shown in a compact error alert in normal mode and logged to the
terminal. Correct the selected source, file, or parameter and render again; a
successful operation clears the alert.

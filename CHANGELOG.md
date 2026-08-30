# Changelog

## 1.3.5 - 2026-08-31

- Replaced the legacy Tkinter window with a responsive Dash workspace for
  bSKAN and VASP volumetric analysis.
- Added interactive STM, apparent-barrier-height, and finite-height local
  work-function workflows with live controls, exact numeric inputs, overlays,
  line profiles, FFT analysis, and independent figure export.
- Generated Chen CURRENT files directly from ASAMPLE, CURSAVE, and INSCAN while
  preserving bSKAN-compatible header and body formatting, so Chen calculations
  no longer require a preliminary Tersoff-Hamann run.
- Added large-file upload progress, bounded surface caching, multi-window
  lifecycle handling, scientific documentation, package-data validation, and
  continuous integration.

## 1.3.4 - 2026-08-25

- Converted constant-current surface indices to physical heights using the
  topmost atomic position, the CURRENT cell height, and `GRID_Z`.
- Selected the outermost tip-side solution of `I(x, y, z) = I_iso` when a
  non-monotonic current profile produces multiple crossings, preventing holes
  in reconstructed STM maps.
- Generalized periodic image reconstruction for rectangular, hexagonal,
  oblique, and monoclinic in-plane cells, including unequal `a` and `b` lengths.
- Matched bSKAN's `1e-6` lattice-classification precision in image repetition,
  atom overlays, and unit-cell overlays.
- Applied contrast through the colormap while retaining the calculated height
  array and color scale in angstrom.
- Added bSKAN and VASP volumetric processing, automatic isosurface generation,
  apparent-barrier analysis, and finite-height local work-function analysis to
  the command-line backend.
- Added structure pre-checking for nonzero symmetry translations and VASP 4
  ASAMPLE generation before bSKAN calculations.

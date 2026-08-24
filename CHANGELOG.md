# Changelog

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

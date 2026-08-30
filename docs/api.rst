Python API
==========

The command-line interfaces are the primary supported interfaces. The functions
below are also maintained as importable Python APIs. Other underscore-prefixed
helpers and GUI callback internals are implementation details.

Symmetry pre-check
------------------

``autobskan.calculation.symmetry_precheck.precheck_bskan_structure(structure, *, type_ids=None, tolerance=..., lattice_symprec=..., surface_normal_tolerance_degrees=...)``
  Return a ``PrecheckResult`` with status, detected operations, translations,
  and a proposed ASE displacement when one exists.

``apply_bskan_origin_shift(atoms, result)``
  Return a translated and periodically wrapped copy for a ``SHIFT_REQUIRED``
  result.

``write_prechecked_poscar(source_path, output_path, *, force=False, ...)``
  Apply, write, and post-validate a proposed shift atomically.

``write_bskan_asample(source_path, output_path=None, *, force=False)``
  Write a VASP 4 count-only ASAMPLE in Direct coordinates while preserving
  existing Selective-dynamics flags or adding ``T T T``.

``validate_existing_scf_for_bskan(scf_directory)``
  Validate an existing SCF structure and raise
  ``ExistingScfRequiresRestartError`` when restarting from pre-check is needed.

Chen CURRENT conversion
-----------------------

``autobskan.calculation.current_header.asample_to_current(asample_path='ASAMPLE', *, cursave_path=None, inscan_path=None, current_path=None, cell_repeat=None, z_vacuum_angstrom=None, delz_bohr=...)``
  Convert Chen ``CURSAVE`` output into an atomically written ``CURRENT`` file
  using bSKAN-compatible header and ``E14.4`` body formatting.

``build_current_head_from_asample(asample_path, grid_shape, *, inscan_path=None, cell_repeat=None, z_vacuum_angstrom=None, scan_origin_bohr=(0, 0), delz_bohr=...)``
  Return only the source-compatible CURRENT header.

Surface quantities
------------------

``autobskan.image.lwfplot.local_work_function_slice(volumetric, height, fermi_level, topmost=0.0)``
  Return the interpolated LOCPOT slice minus ``E_F``.

``fermi_energy_from_outcar(path)`` / ``resolve_fermi_level(volume_path, fermi_level=None)``
  Read the last reported Fermi level or resolve an explicit value.

``apparent_barrier_from_vasp_density(volumetric, height, topmost=0.0, fit_radius=0.5)``
  Fit :math:`\Phi_{\rm app}` on a VASP density grid.

``apparent_barrier_from_current(current, height=None, fit_radius=0.5, isosurface=None, reference='CONSTANT HEIGHT')``
  Fit :math:`\Phi_{\rm app}` from bSKAN CURRENT at a constant-height or
  constant-current reference.

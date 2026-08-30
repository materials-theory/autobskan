Contributing to AutoBSKAN
========================

AutoBSKAN welcomes reproducible bug reports and focused pull requests. Use the
GitHub issue templates for problems and feature proposals. General usage
questions may also be opened as issues when the documentation does not answer
them.

Development setup
-----------------

Use Python 3.10 or newer in an isolated environment::

  python -m pip install -e ".[dev,docs]"
  python -m ruff check autobskan/gui autobskan/main.py tests
  python -m pytest
  python -m sphinx -W -b html docs docs/_build/html

Before opening a pull request, also verify the package metadata::

  python -m build
  python -m twine check dist/*

Contribution scope
------------------

Keep changes local to the affected workflow. Do not combine a numerical change,
GUI redesign, and unrelated refactor in one pull request. Update user
documentation and ``CHANGELOG.md`` when public behaviour changes.

Scientific changes require all of the following:

* a written equation, upstream source location, or publication establishing the
  intended convention;
* a small redistributable input fixture;
* a numerical regression against the original bSKAN/VASP-derived result or an
  independent implementation;
* explicit tolerances and units; and
* tests for invalid domains rather than silent fallback to another method.

Do not commit licensed VASP binaries, POTCAR files, private structures, or large
raw output. Reduce a failure to a synthetic or redistributable fixture.

GUI changes
-----------

GUI work should preserve the physical scalar field independently of browser
display optimisation. Check at least one desktop viewport and one window below
1260 px. Interactive controls must remain keyboard reachable, exact numeric
input must agree with sliders, and long labels must not overlap. Slow callbacks
need a visible busy state; failures need a visible user-facing message in normal
mode and a terminal traceback for diagnosis.

Versioning and review
---------------------

AutoBSKAN uses three-part versions. Increment the patch component for each
released correction (for example ``1.2.26`` to ``1.2.27``), and keep
``autobskan.__version__``, package metadata, ``CITATION.cff``, and the changelog
consistent.

Pull requests are reviewed for behavioural scope, test evidence, scientific
provenance, backward compatibility, documentation, and generated-file churn.

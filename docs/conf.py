from __future__ import annotations

from importlib.metadata import version


project = "AutoBSKAN"
author = "AutoBSKAN contributors"
copyright = "2026, AutoBSKAN contributors"
release = version("autobskan")
version = ".".join(release.split(".")[:2])

extensions: list[str] = []
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
html_theme = "alabaster"
html_title = f"AutoBSKAN {release}"

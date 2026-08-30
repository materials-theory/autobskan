from pathlib import Path

from setuptools import find_namespace_packages, setup


ROOT = Path(__file__).resolve().parent


setup(
    name="autobskan",
    version="1.3.5",
    description="STM, apparent-barrier, and local work-function analysis",
    long_description=(ROOT / "README.rst").read_text(encoding="utf-8"),
    long_description_content_type="text/x-rst",
    url="https://github.com/materials-theory/autobskan",
    author="Giyeok Lee",
    author_email="lgy4230@yonsei.ac.kr",
    license="MIT",
    packages=find_namespace_packages(include=["autobskan", "autobskan.*"]),
    package_data={
        "": ["elements_vesta.ini", "bskan.in", "bskan.in_example"],
        "autobskan.gui": ["assets/*.css", "assets/*.js"],
    },
    include_package_data=True,
    zip_safe=False,
    python_requires=">=3.10",
    keywords="STM bSKAN DFT VASP ASE chen bardeen",
    install_requires=[
        "numpy",
        "scipy",
        "ase",
        "spglib>=2.0",
        "matplotlib",
        "Pillow",
        "tqdm",
        "dash",
        "dash-mantine-components",
        "plotly>=6.1,<6.8",
        "kaleido>=1.0,<2",
    ],
    extras_require={
        "test": ["pytest>=8"],
        "dev": ["build>=1", "pytest>=8", "ruff>=0.9", "twine>=5"],
        "docs": ["sphinx>=8"],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Operating System :: MacOS",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    project_urls={
        "Documentation": "https://github.com/materials-theory/autobskan#readme",
        "Issues": "https://github.com/materials-theory/autobskan/issues",
        "Source": "https://github.com/materials-theory/autobskan",
    },
    entry_points={
        "console_scripts": [
            "autobskan = autobskan.main:main",
            "autobskan-gui = autobskan.gui.frontend:main",
            "autobskan-post = autobskan.cli.post_processing_cli:main",
            "autobskan-collect = autobskan.cli.file_collector:main",
        ]
    },
)

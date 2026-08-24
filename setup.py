import io

from setuptools import find_packages, setup

from autobskan import __version__


# Read in the README for the long description on PyPI
def long_description():
    with io.open('README.rst', 'r', encoding='utf-8') as f:
        readme = f.read()
    return readme

setup(
    name='autobskan',
    version=__version__,
    description='STM and local electronic-property analysis for bSKAN and VASP',
    long_description=long_description(),
    long_description_content_type='text/x-rst',
    url='https://github.com/materials-theory/autobskan',
    author='Giyeok Lee',
    author_email='lgy4230@yonsei.ac.kr',
    license='MIT',
    packages=find_packages(),
    package_data={
        'autobskan.image': ['elements_vesta.ini'],
        'autobskan.io': ['*.in', 'bskan.in_example'],
    },
    include_package_data=True,
    zip_safe=False,
    python_requires='>=3.10',
    keywords='STM bSKAN DFT VASP ASE Chen Bardeen',
    install_requires=[
        'numpy',
        'scipy',
        'ase',
        'spglib>=2.0',
        'matplotlib',
        'Pillow',
        'tqdm',
    ],
    extras_require={'test': ['pytest>=8']},
    classifiers=[
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',
    ],
    entry_points={'console_scripts': [
        'autobskan = autobskan.main:main',
        'autobskan-gui = autobskan.gui.frontend:main',
        'autobskan-post = autobskan.cli.post_processing_cli:main',
        'autobskan-collect = autobskan.cli.file_collector:main',
    ]},
)

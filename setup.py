import io
from setuptools import find_packages, setup, Extension
from autobskan.main import __version__


# Read in the README for the long description on PyPI
def long_description():
    with io.open('README.rst', 'r', encoding='utf-8') as f:
        readme = f.read()
    return readme

setup(name='autobskan',
      version=__version__,
      description='Image generation and post processing code of bSKAN',
      long_description=long_description(),
      long_description_content_type='text/x-rst',
      url='https://github.com/materials-theory/autobskan',
      author='Giyeok Lee',
      author_email='lgy4230@yonsei.ac.kr',
      license='MIT',
      packages=find_packages(), # this cannot find packages without __init__ (Important!!)
      package_data={'':["elements_vesta.ini"]},
      include_package_data = True,
      zip_safe = False,
      python_requires='>=3.10',
      keywords='STM bSKAN DFT vasp ase chen bardeen',
      classifiers=[
          'Programming Language :: Python :: 3 :: Only',
          'Programming Language :: Python :: 3.10',
          'Programming Language :: Python :: 3.11',
          'Programming Language :: Python :: 3.12',
          'Programming Language :: Python :: 3.13'
          ],
      install_requires=['numpy', 'scipy', 'ase', 'spglib>=2.0', 'matplotlib',
                        'Pillow', 'tqdm'],
      entry_points = {'console_scripts' : [
      'autobskan = autobskan.main:main',
      'autobskan-gui = autobskan.gui.frontend:main',
      'autobskan-post = autobskan.cli.post_processing_cli:main',
      'autobskan-collect = autobskan.cli.file_collector:main']}
      )

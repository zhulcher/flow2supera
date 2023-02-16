from skbuild import setup  # This line replaces 'from setuptools import setup'
import argparse

import io,os,sys
this_directory = os.path.abspath(os.path.dirname(__file__))
with io.open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="larnd2supera",
    version="0.0.1",
    #cmake_source_dir='src/',
    include_package_data=True,
    #cmake_args=[
    #    #'-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON',
    #    '-DCMAKE_OSX_DEPLOYMENT_TARGET:STRING=10.9',
    #],
    author=['Kazuhiro Terao, Zach Hulcher, Andrew Morgan'],
    author_email='kterao@slac.stanford.edu, zhulchero@slac.stanford.edu, andrew.mogan@colostate.edu',
    description='Supera interface for larnd-sim data files',
    license='MIT',
    keywords='supera larnd-sim larcv larcv3 neutrinos deep learning lartpc_mlreco3d',
    project_urls={
        'Source Code': 'https://github.com/DeepLearnPhysics/larnd2supera'
    },
    url='https://github.com/DeepLearnPhysics/larnd2supera',
    scripts=['bin/run_larnd2supera.py'],
    packages=['larnd2supera','larnd2supera.pdg_data','larnd2supera.config_data'],
    package_dir={'': 'python'},
    package_data={'larnd2supera': ['pdg_data/pdg.npz']},
    install_requires=[
        'numpy',
        'scikit-build',
        'supera',
        'edep2supera',
    ],
    long_description=long_description,
    long_description_content_type='text/markdown',
)

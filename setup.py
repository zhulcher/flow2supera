from setuptools import setup  # This line replaces 'from setuptools import setup'
import argparse

import io,os,sys
this_directory = os.path.abspath(os.path.dirname(__file__))
with io.open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="flow2supera",
    version="0.0.1",
    include_package_data=True,
    author=['Kazuhiro Terao, Zach Hulcher, Andrew Mogan, Sindhujha Kumaran'],
    author_email='kterao@slac.stanford.edu, zhulcher@slac.stanford.edu, andrew.mogan@colostate.edu, s.kumaran@uci.edu',
    description='Supera interface for ndlar_flow data files',
    license='MIT',
    keywords='supera larnd-sim larcv larcv3 neutrinos deep learning lartpc_mlreco3d ndlar_flow h5flow',
    project_urls={
        'Source Code': 'https://github.com/DeepLearnPhysics/flow2supera'
    },
    url='https://github.com/DeepLearnPhysics/flow2supera',
    scripts=['bin/run_flow2supera.py'],
    packages=['flow2supera','flow2supera.config_data'],
    package_dir={'': 'src'},
    package_data={'flow2supera': ['config_data/*.yaml']},
    install_requires=[
        'numpy',
        #'scikit-build',
        'supera',
        'edep2supera',
        'h5flow',
    ],
    long_description=long_description,
    long_description_content_type='text/markdown',
)

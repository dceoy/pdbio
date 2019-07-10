#!/usr/bin/env python

from setuptools import find_packages, setup

from pdbio import __version__

setup(
    name='pdbio',
    version=__version__,
    description='Pandas-based Data Handler for VCF and BED Files',
    packages=find_packages(),
    author='Daichi Narushima',
    author_email='dnarsil+github@gmail.com',
    url='https://github.com/dceoy/pdbio',
    include_package_data=True,
    install_requires=['docopt', 'pandas'],
    entry_points={
        'console_scripts': ['pdbio=pdbio.cli:main'],
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Plugins',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    long_description="""\
# pdbio

Pandas-based Data Handler for VCF and BED Files
"""
)

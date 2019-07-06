#!/usr/bin/env python

from setuptools import find_packages, setup

from pdvcf import __version__

setup(
    name='pdvcf',
    version=__version__,
    description='Pandas-based Data Handlers for DNA-sequencing',
    packages=find_packages(),
    author='Daichi Narushima',
    author_email='dnarsil+github@gmail.com',
    url='https://github.com/dceoy/pdvcf',
    include_package_data=True,
    install_requires=['docopt', 'pandas'],
    entry_points={
        'console_scripts': ['pdvcf=pdvcf.cli:main'],
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
# pdvcf

Pandas data handler for VCF files
"""
)

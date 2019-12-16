#!/usr/bin/env python

from setuptools import find_packages, setup

from pdbio import __version__

with open('README.md', 'r') as f:
    long_description = f.read()

setup(
    name='pdbio',
    version=__version__,
    author='Daichi Narushima',
    author_email='dnarsil+github@gmail.com',
    description='Pandas-based Data Handler for VCF, BED, and SAM Files',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/dceoy/pdbio',
    packages=find_packages(),
    include_package_data=True,
    install_requires=['docopt', 'pandas'],
    entry_points={'console_scripts': ['pdbio=pdbio.cli:main']},
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Plugins',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    python_requires='>=3.6'
)

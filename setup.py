#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup

# Define package info
name = "pycoMeth"
description = "DNA methylation analysis downstream to Nanopolish for Oxford Nanopore DNA sequencing datasets"
with open("README.md", "r") as fh:
    long_description = fh.read()

# Collect info in a dictionnary for setup.py
setup(
    name=name,
    description="Differential methylation calling suite for Nanopore methylation calls PycoMeth",
    version="2.1.1",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/snajder-r/pycoMeth,",
    author="Rene Snajder",
    author_email="r.snajder@dkfz-heidelberg.de",
    license="GPL",
    python_requires=">=3.7",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python :: 3"
    ],
    install_requires=[
        "numpy>=1.19.2",
        "scipy==1.4.1",
        "statsmodels>=0.11.1",
        "pandas>=1.1.3",
        "Jinja2>=2.11.1",
        "plotly>=4.7.1",
        "pyfaidx>=0.5.8",
        "tqdm>=4.60.0",
        "colorlog>=4.1.0",
        "nbformat>=4.2.0",
        "meth5>=1.0.1"
    ],
    packages=["pycoMeth", "pycoMeth.meth_seg"],
    package_dir={"pycoMeth": "pycoMeth", "pycoMeth.meth_seg":"pycoMeth/meth_seg"},
    package_data={name: ["templates/*"]},
    entry_points={"console_scripts": ["pycometh=pycoMeth.__main__:main", "pycoMeth=pycoMeth.__main__:main"]},
)

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
    version="2.2.1",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/snajder-r/pycoMeth,",
    author="Rene Snajder",
    author_email="r.snajder@dkfz-heidelberg.de",
    license="GPL",
    python_requires="==3.7",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3"
    ],
    install_requires=[
        "numpy==1.22.2",
        "scipy==1.4.1",
        "statsmodels==0.13.2",
        "pandas==1.4.1",
        "Jinja2==3.0.3",
        "plotly==5.6.0",
        "pyfaidx==0.6.4",
        "tqdm==4.62.3",
        "colorlog==6.6.0",
        "nbformat==5.1.3",
        "meth5==1.1.1"
    ],
    packages=["pycoMeth", "pycoMeth.meth_seg"],
    package_dir={"pycoMeth": "pycoMeth", "pycoMeth.meth_seg":"pycoMeth/meth_seg"},
    package_data={name: ["templates/*"]},
    entry_points={"console_scripts": ["pycometh=pycoMeth.__main__:main", "pycoMeth=pycoMeth.__main__:main"]},
)

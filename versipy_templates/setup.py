#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup

# Define package info
description = "DNA methylation analysis downstream to Nanopolish for Oxford Nanopore DNA sequencing datasets"
with open("README.md", "r") as fh:
    long_description = fh.read()

# Collect info in a dictionnary for setup.py
setup(
    name="__package_name__",
    description="__package_description__",
    version="__package_version__",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="__package_url__",
    author="__author_name__",
    author_email="__author_email__",
    license="__package_licence__",
    python_requires=">=__minimal_python__",
    classifiers=[
        __@{,
        ::"classifiers"}__
    ],
    install_requires=[
        __@{,
        ::"dependencies"}__
    ],
    packages=["__package_name__"],
    package_dir={"__package_name__": "__package_name__"},
    package_data={"__package_name__": ["templates/*"]},
    entry_points={"console_scripts": ["pycoMeth=pycoMeth.__main__:main"]},
)

#!/bin/bash
# -*- coding: utf-8 -*-

repostring=""
if [[ $# -eq 1 ]]; then
    if [ "$1" = "test" ]; then
        repostring="--repository testpypi"
    fi
fi

echo "compile package from setup.py"
python setup.py sdist

echo "Uploading to pypi..."
twine upload $repostring dist/*

if [[ ! -z $repostring ]]; then
	echo "Testing complete - not attempting to upload to conda"
	exit
fi

echo "Build noarch package for conda..."
conda-build meta.yaml --python 3.7 --output-folder conda_build -c bioconda -c conda-forge -c plotly -c snajder-r --no-include-recipe --no-test

echo "Logging in to conda..."
anaconda login

echo "Deploying to conda..."
anaconda upload conda_build/**/*.tar.bz2

echo "Cleaning up"

rm -Rf dist
rm -Rf conda_build
rm -Rf build
rm -Rf *.egg-info
rm -Rf .eggs

exit 0

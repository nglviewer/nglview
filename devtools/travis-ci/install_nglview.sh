#!/bin/sh

if [ "$CONDA" = "yes" ]; then
    echo "test conda build"
    conda config --add channels conda-forge
    conda build devtools/travis-ci/conda-recipe --py=3.5
    conda install /home/travis/miniconda/conda-bld/linux-64/nglview-*.bz2
else
    echo "test pip build"
    python setup.py sdist
    pip install dist/*gz
fi

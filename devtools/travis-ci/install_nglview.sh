#!/bin/sh

if [ "$CONDA" = "yes" ]; then
    echo "test conda build"
    conda build devtools/travis-ci/conda-recipe
    # conda install 
else
    echo "test pip build"
    python setup.py sdist
    pip install dist/*gz
fi

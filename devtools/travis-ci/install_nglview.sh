#!/bin/sh

if [ "$CONDA" = "yes" ]; then
    echo "test conda build"
    pip uninstall ipywidgets -y
    conda remove ipywidgets --yes
    conda config --add channels conda-forge
    conda build devtools/travis-ci/conda-recipe --py=3.5
    conda install /home/travis/miniconda/conda-bld/linux-64/nglview-*.bz2
    # force to install new ipywidgets
    # I really don't know why travis still picks up ipywidgets 4.1.1
    # although we did specify ipywidgets > 5.1.1
    conda install ipywidgets -c conda-forge --yes
else
    echo "test pip build"
    python setup.py sdist
    pip install dist/*gz
fi

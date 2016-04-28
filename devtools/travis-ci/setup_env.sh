#!/bin/sh

if [ "$TRAVIS_OS_NAME" = "osx" ]; then
    wget http://repo.continuum.io/miniconda/Miniconda-3.7.0-MacOSX-x86_64.sh -O miniconda.sh;
else
    wget http://repo.continuum.io/miniconda/Miniconda-3.7.0-Linux-x86_64.sh -O miniconda.sh;
fi

bash miniconda.sh -b

export PATH=$HOME/miniconda/bin:$PATH
# install stable version
pip install conda

conda install --yes conda-build jinja2 anaconda-client pip

# create myenv
conda create -y -n myenv python=$PYTHON_VERSION jupyter notebook nose numpy mock

source activate myenv
pip install pytraj
conda install parmed -c ambermd
if [ "$PYTHON_VERSION" = "2.7" ]; then
    conda install mdanalysis -c MDAnalysis
fi
# conda install mdtraj -c omnia --yes

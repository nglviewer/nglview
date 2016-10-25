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

conda update -n root conda-build --yes
conda install --yes conda-build jinja2 anaconda-client pip

# create myenv
conda create -y -n myenv python=$PYTHON_VERSION jupyter notebook nose numpy mock coverage cython netcdf4

source activate myenv

# pytraj
pip install pytraj

# mdanalysis
conda config --add channels MDAnalysis
conda install mdanalysis -c kain88-de --yes

# mdtraj
conda install mdtraj -c omnia --yes

# ParmEd
pip install https://github.com/ParmEd/ParmEd/archive/2.5.1.tar.gz

# ase
pip install https://github.com/hainm/ase/archive/3.11.0.tar.gz

# simpletraj
pip install git+https://github.com/arose/simpletraj

# pytest
pip install pytest
pip install pytest-cov

# coveralls
pip install coveralls

# rdkit
conda install rdkit -c rdkit --yes

if [ "$TEST_NOTEBOOK" = "yes" ]; then
    npm install -g nightwatch
fi

# download ngl data
git clone https://github.com/arose/ngl
mv ngl ../

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
conda create -y -n myenv python=$PYTHON_VERSION jupyter notebook nose numpy mock coverage cython netcdf4

source activate myenv
pip install pytraj
conda install parmed -c ambermd --yes
if [ "$PYTHON_VERSION" = "2.7" ]; then
    conda install mdanalysis -c MDAnalysis --yes
else
   git clone https://github.com/MDAnalysis/mdanalysis
   cd mdanalysis/package
   python setup.py install
   cd ../../
fi
conda install mdtraj -c omnia --yes

# simpletraj
pip install git+https://github.com/arose/simpletraj

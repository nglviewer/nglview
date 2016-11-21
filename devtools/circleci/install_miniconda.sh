#!/bin/sh

miniconda=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
wget $miniconda -O miniconda.sh
bash miniconda.sh -b
export PATH=/root/miniconda3/bin:$PATH
conda update --yes --all
conda install --yes conda-build anaconda-client numpy matplotlib
conda info
conda install --yes jupyter notebook nose numpy mock coverage cython netcdf4
# pytraj
pip install pytraj
# mdanalysis
conda config --add channels MDAnalysis
conda install mdanalysis -c kain88-de --yes
# mdtraj
conda install mdtraj -c omnia --yes
# ParmEd
conda install -c ambermd parmed --yes
# ase
pip install https://github.com/hainm/ase/archive/3.11.0.tar.gz
# simpletraj
pip install git+https://github.com/arose/simpletraj
# rdkit
# conda install rdkit -c rdkit --yes

# update a bunch of stuff
conda install --yes -c conda-forge ipywidgets notebook pyzmq

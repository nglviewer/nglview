#!/bin/sh

# FIXME: Update to requirements.txt once those packages have binary
# distribution on pypi.
conda install mdanalysis=0.16.2 -c conda-forge --yes
conda install mdtraj=1.8.0 -c conda-forge --yes
conda install rdkit -c rdkit --yes
conda install simpletraj -c omnia --yes

# > 100 MB
# conda install ambertools=17 -c http://ambermd.org/downloads/ambertools/conda/ --yes

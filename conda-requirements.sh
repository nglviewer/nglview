#!/bin/sh

# FIXME: Update to requirements.txt once those packages have binary
# distribution on pypi.
conda install mdanalysis -c conda-forge --yes
conda install mdtraj=1.8.0 -c conda-forge --yes
conda install rdkit -c rdkit --yes

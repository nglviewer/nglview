#!/bin/sh

# FIXME: Update to requirements.txt once those packages have binary
# distribution on pypi.
conda install mdanalysis -c kain88-de --yes
conda install mdtraj=1.8.0 -c omnia --yes
conda install rdkit=2017.03.3 -c rdkit --yes

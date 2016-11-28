#!/bin/sh

# require: npm
# python setup.py install --npm
pip install . -v
jupyter labextension install nglview --py --sys-prefix
jupyter labextension enable nglview --py --sys-prefix

#!/bin/sh
python setup.py install --npm
jupyter labextension install nglview --py --sys-prefix
jupyter labextension enable nglview --py --sys-prefix

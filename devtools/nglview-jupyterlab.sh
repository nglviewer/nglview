#!/bin/sh

nglview_src=.
env=jupyterlab_0252

conda create -n $env python=3.6 numpy -y
source activate $env

pip install ipywidgets==7.0.0b0
jupyter nbextension enable --py --sys-prefix widgetsnbextension

conda install jupyterlab==0.25.2 -c conda-forge -y
jupyter labextension install @jupyter-widgets/jupyterlab-manager

cd $nglview_src

# nglview-js-widgets
cd js
npm install --save
cd ../

# nglview
python setup.py install

# nglview-jupyterlab
jupyter-labextension install jslab

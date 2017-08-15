#!/bin/sh

nglview_src=.
env=jupyterlab_0252
jupyterlab_version=0.25.2
python_version=3.6
notebook_version=5.0.0
ipywidgets_version=7.0.0b0
lab_manager_version=0.24.7

conda create -n $env python=${python_version} numpy nomkl -y
source activate $env

conda install notebook==${notebook_version} -c conda-forge
pip install ipywidgets==${ipywidgets_version}
pip install widgetsnbextension==3.0.0b2 --upgrade
jupyter nbextension install --py --sys-prefix widgetsnbextension
jupyter nbextension enable --py --sys-prefix widgetsnbextension

conda install jupyterlab==${jupyterlab_version} -c conda-forge -y
jupyter labextension install @jupyter-widgets/jupyterlab-manager@${lab_manager_version}

cd $nglview_src

# nglview-js-widgets
cd js
npm install --save
cd ../

# nglview
python setup.py install

# nglview-jupyterlab
jupyter-labextension install jslab

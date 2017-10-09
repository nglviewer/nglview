#!/bin/sh

nglview_src=`pwd`
env=jupyterlab_0265
jupyterlab_version=0.26.5
python_version=3.6
notebook_version=5.0.0
lab_manager_version=0.27.0 # need to match to ./jslab/package.json
ipywidgets_version=7.0.0
widgetsnbextension_version=3.0.0

echo "nglview root folder: $nglview_src"
echo "Creating env $env"

conda create -n $env python=${python_version} numpy nomkl -y
source activate $env
conda install setuptools -c conda-forge --force -y
conda install notebook==${notebook_version} -c conda-forge -y
pip install ipywidgets==${ipywidgets_version}
pip install widgetsnbextension==${widgetsnbextension_version} --upgrade
pip install pytraj # comment if using windows
jupyter nbextension install --py --sys-prefix widgetsnbextension
jupyter nbextension enable --py --sys-prefix widgetsnbextension

conda install jupyterlab==${jupyterlab_version} -c conda-forge -y
jupyter labextension install @jupyter-widgets/jupyterlab-manager@${lab_manager_version}

cd $nglview_src

## Developer should do this
## nglview-js-widgets
# cd js
# npm install && npm publish
# cd ../

# nglview
python setup.py install

# nglview-jupyterlab
jupyter-labextension install ./jslab

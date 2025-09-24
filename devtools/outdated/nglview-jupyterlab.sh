#!/bin/sh

nglview_src=`pwd`
env=jupyterlab_0_31
jupyterlab_version=0.31
lab_manager_version=${jupyterlab_version} # need to match to ./jslab/package.json too too.
python_version=3.6
notebook_version=5.0.0
ipywidgets_version=7.1.2
# widgetsnbextension_version=3.0.0

echo "nglview root folder: $nglview_src"
echo "Creating env $env"

conda create -n $env python=${python_version} numpy<2.3 nomkl -y
source activate $env
conda install setuptools -c conda-forge --force -y
conda install notebook==${notebook_version} -c conda-forge -y
pip install ipywidgets==${ipywidgets_version}
# pip install widgetsnbextension==${widgetsnbextension_version} --upgrade
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

conda create -n lab python=3.8 -y 
source  activate lab
conda install jupyterlab=3 -c conda-forge -y
pip install jupyterlab_widgets
pip install nglview # Python part
pip install https://github.com/nglviewer/nglview/releases/download/v2.7.7/nglview-js-widgets-2.7.7.tar.gz # JS part for Lab

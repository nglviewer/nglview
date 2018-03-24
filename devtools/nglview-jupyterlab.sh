conda create -n lab python=3.6 -y
source activate lab
conda install jupyterlab=0.31.12 -y -c conda-force
pip install nglview==1.1.2
nglview enable # might need this to enable nglview-js-widgets extension for notebook
conda install ipywidgets=7.1.2 -c conda-forge -y
jupyter-labextension install @jupyter-widgets/jupyterlab-manager@0.33.2
jupyter-labextension install nglview-js-widgets@1.1.2

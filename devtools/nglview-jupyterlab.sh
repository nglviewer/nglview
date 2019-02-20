pythonversion=3.7
nglviewversion=1.1.8
labversion=0.34.12
labmanagerversion=0.37.4
ipywidgetsversion=7.4.2

conda create -n lab python=$pythonversion -y
source activate lab
conda install ipywidgets=$ipywidgetsversion -c conda-forge -y
pip install nglview==$nglviewversion
nglview enable # might need this to enable nglview-js-widgets extension for notebook
conda install jupyterlab=$labversion  -y -c conda-forge
jupyter-labextension install @jupyter-widgets/jupyterlab-manager@$labmanagerversion
jupyter-labextension install nglview-js-widgets@$nglviewversion

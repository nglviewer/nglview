# docker build . -f Dockerfile-jupyterlab -t hainm/juypterlab:0.25.2
#
FROM  continuumio/miniconda3

RUN conda install notebook -c conda-forge -y
RUN pip install ipywidgets==7.0.0b0
RUN jupyter nbextension enable --py --sys-prefix widgetsnbextension

RUN conda install nodejs -c conda-forge -y
RUN conda install jupyterlab==0.25.2 -c conda-forge -y
RUN jupyter labextension install @jupyter-widgets/jupyterlab-manager # 0.24.7

# Expose Jupyter port & cmd
EXPOSE 8888
RUN mkdir -p /opt/app/data
CMD jupyter lab --ip=* --port=8888 --no-browser --notebook-dir=/opt/app/data --allow-root --NotebookApp.iopub_data_rate_limit=100000000

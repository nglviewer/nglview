# docker build . -t hainm/nglview:1.0.a0
FROM  continuumio/miniconda

ADD . /opt/app/nglview

RUN conda install notebook -c conda-forge -y
RUN pip install ipywidgets==7.0.0b0
RUN jupyter nbextension enable --py --sys-prefix widgetsnbextension

RUN conda install nodejs -c conda-forge -y
RUN conda install jupyterlab==0.25.2 -c conda-forge -y
RUN jupyter labextension install @jupyter-widgets/jupyterlab-manager # 0.24.7

# nglview-js-widgets
RUN cd /opt/app/nglview/js
RUN npm install --save
RUN cd /opt/app/nglview

# nglview
RUN python setup.py install

# nglview-jupyterlab
RUN jupyter-labextension install jslab

# Expose Jupyter port & cmd
EXPOSE 8888
RUN mkdir -p /opt/app/data
CMD jupyter lab --ip=* --port=8888 --no-browser --notebook-dir=/opt/app/data

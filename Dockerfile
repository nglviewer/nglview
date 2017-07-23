# (optional) docker build . -f Dockerfile-jupyterlab -t hainm/juypterlab:0.25.2
# docker build . -t hainm/nglview:1.0.a0

# How to run?
# docker -it --rm -p 8888:8888 hainm/nglview:1.0.a0

FROM  hainm/jupyterlab:0.25.2
RUN pip install numpy
RUN pip install pytraj

ADD . /opt/app/nglview

# nglview-js-widgets
RUN cd /opt/app/nglview/js && npm install . --save
RUN cd /opt/app/nglview && python setup.py install

# nglview-jupyterlab
RUN cd /opt/app/nglview && jupyter-labextension install jslab
CMD jupyter lab --ip=* --port=8888 --no-browser --notebook-dir=/opt/app/data --allow-root --NotebookApp.iopub_data_rate_limit=100000000

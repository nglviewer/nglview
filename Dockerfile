# docker build . -t hainm/nglview:1.0.a0
FROM  hainm/jupyterlab-dev

ADD . /opt/app/nglview
RUN pip install numpy

# nglview-js-widgets
RUN cd /opt/app/nglview/js && npm install . --save
RUN cd /opt/app/nglview && python setup.py install

# nglview-jupyterlab
RUN cd /opt/app/nglview && jupyter-labextension install jslab

# (optional) docker build . -f Dockerfile-jupyterlab -t hainm/juypterlab:0.25.2
# docker build . -t hainm/nglview:1.0.a0

# How to run?
# docker run -it --rm -p 8888:8888 hainm/nglview:1.0.a0

FROM hainm/jupyterlab:0.25.2
ADD ./jslab /opt/jslab

RUN pip install numpy==1.13.1
RUN pip install pytraj==2.0.2
RUN pip install matplotlib==2.0.2
RUN pip install moviepy==0.2.2.11
RUN pip install imageio==1.6
RUN pip install pillow==4.2.1
RUN pip install nglview==1.0.a0
RUN nglview install
RUN nglview enable
RUN jupyter-labextension install /opt/jslab
RUN mkdir -p /root/.imageio/ffmpeg/
RUN curl https://github.com/imageio/imageio-binaries/raw/master/ffmpeg/ffmpeg.linux64 --output /root/.imageio/ffmpeg/ffmpeg.linux64

CMD jupyter lab --ip=* --port=8888 --no-browser --notebook-dir=/opt/app/data --allow-root --NotebookApp.iopub_data_rate_limit=100000000

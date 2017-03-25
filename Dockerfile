FROM  jupyter/notebook

# 2.7
RUN python -m pip install pip --upgrade
RUN python -m pip install numpy
RUN python -m pip install nglview==0.6.2.3
RUN python -m pip install ipywidgets==5.2.2
RUN python -m pip install pytraj==1.0.9
RUN python -m jupyter nbextension enable --py --sys-prefix widgetsnbextension
RUN python -m jupyter nbextension install --py --sys-prefix nglview
RUN python -m jupyter nbextension enable --py --sys-prefix nglview

# 3.4
RUN python3 -m pip install pip --upgrade
RUN python3 -m pip install numpy
RUN python3 -m pip install nglview==0.6.2.3
RUN python3 -m pip install ipywidgets==5.2.2
RUN python3 -m pip install pytraj==1.0.9
RUN python3 -m jupyter nbextension enable --py --sys-prefix widgetsnbextension
RUN python3 -m jupyter nbextension install --py --sys-prefix nglview
RUN python3 -m jupyter nbextension enable --py --sys-prefix nglview

CMD ["jupyter", "notebook", "--no-browse"]

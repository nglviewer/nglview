FROM  jupyter/notebook
RUN python -m pip install pip --upgrade
RUN python -m pip install numpy
RUN python -m pip install nglview
RUN python -m pip install ipywidgets
RUN python -m pip install pytraj
RUN python -m jupyter nbextension enable --py --sys-prefix widgetsnbextension
RUN python -m jupyter nbextension install --py --sys-prefix nglview
RUN python -m jupyter nbextension enable --py --sys-prefix nglview

CMD ["jupyter", "notebook", "--no-browse"]

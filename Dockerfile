FROM  condaforge/linux-anvil
ADD . /root/nglview
RUN     cd /root/nglview \
    &&  sh devtools/circleci/install_miniconda.sh \
    &&  /opt/conda/bin/python setup.py install
# CMD ["nglview", "demo"]

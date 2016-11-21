# FROM  ambermd/manylinux-extra
FROM  condaforge/linux-anvil

# To get the AmberTools16.tar.bz file, fill out the form
# at the site below and click Download.
ADD . /root/nglview

RUN     cd /root/nglview \
    &&  sh devtools/circleci/install_miniconda.sh \
    &&  /opt/conda/bin/python setup.py install
# CMD ["nglview", "demo"]

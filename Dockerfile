FROM  ambermd/manylinux-extra

# To get the AmberTools16.tar.bz file, fill out the form
# at the site below and click Download.
ADD . /root/nglview

RUN     cd /root/nglview \
    &&  sh devtools/circleci/install_miniconda.sh \
    &&  /root/miniconda3/bin/python setup.py install
# CMD ["nglview", "demo"]

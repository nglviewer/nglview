FROM  ambermd/manylinux-extra

# To get the AmberTools16.tar.bz file, fill out the form
# at the site below and click Download.
ADD ../nglview /root/

RUN     cd /root/nglview \
    &&  sh devtools/circleci/install_miniconda.sh \
    &&  python setup.py install

General
=======

-  Please follow `pep8 <https://www.python.org/dev/peps/pep-0008/>`__ to
   format Python code.

Install all required packages for testing
=========================================

.. code:: bash

    pip install -r pip-requirements-test.txt
    sh conda-requirements-test.sh

Making a Release
================

We use `Versioneer <https://github.com/warner/python-versioneer>`__ to
automatically update the version string (of a release but also in
development). This means for a release a new git tag should be created.
The tag should be of the form vX.Y or vX.Y.Z and generally follow
`pep440 <https://www.python.org/dev/peps/pep-0440/>`__ with a prefixed
"v".

.. code:: bash

    git tag  -a vX.Y -m "version X.Y"
    git push
    git push origin --tags

    python setup.py sdist
    # Upload with twine:
    twine upload dist/nglview*gz

Install
=======

Install developer mode so you don't need to re-install after chaning
your code.

.. code:: python

    pip install -e .

Test
====

.. code:: bash

    pytest -vs nglview/tests/

Sync versions
=============

There is not yet an elegant way to sync version in nglview/widget.py;
js/package.json and jslab/package.json So make sure to update those

.. code:: python

    nglview/wiget.py: __frontend_version__ = '0.5.4-dev.13'
    js/package.json:              "version": "0.5.4-dev.13",
    jslab/package.json: "@jupyter-widgets/jupyterlab-manager": "^0.25.11", # make sure this compat with ipywidgets
    jslab/package.json: "nglview-js-widgets": "0.5.4-dev.13"

Update Javascript build
=======================

.. code:: bash

    # install nodejs
    # conda install -c javascript nodejs

    cd js
    npm install

    # known working versions:
    # node: v6.6.0
    # npm: 4.4.4

    # then python setup.py install or pip install -e . # development, you can edit the source code without re-installing

    # Tips
    # - After changing js code:
    #     - Kernel --> Restart and Clear Output
    #     - Refresh webpage (F5)

    # run quick test in terminal
    nglview demo

My (Hai) workflow
=================

.. code:: bash

    cd /nglview/root/folder

    # install development version for Python code
    # so we can update Python source code without reinstall
    pip install -e .
    # Now you can start changing Python source code.

    # make symlink the js code (nglview/static/*js)
    # to $PREFIX/share/jupyter/nbextensions
    # Example of $PREFIX: $HOME/miniconda3/
    # Double-check
    # $ ll $HOME/miniconda3/share/jupyter/nbextensions/
    # Will see something like 
    # nglview-js-widgets@ -> $HOME/3d/nglview/nglview/static

    nglview install --symlink
    nglview enable

    # Now, you can update JS code in js/src folder
    cd js
    npm install

    # If your notebook is openning and you did above step, you need to clear the web cache via two steps
    - Restart your notebook (Kernel -> Restart and Clear Output)
    - Refresh browser (F5)

Using ``NGL`` locally
=====================

1. Change ``var NGL = require('ngl');`` to
   ``var NGL = require('./ngl');``
   https://github.com/arose/nglview/blob/master/js/src/widget\_ngl.js#L2

2. Then, `build
   NGL <https://github.com/arose/ngl/blob/master/DEVELOPMENT.md#building>`__,
   then copy ``ngl.js`` (or ``ngl.dev.js``) to ``nglview/js/src/``

3. Rebuild js code

   ::

       cd js
       npm install
       nglview install # install updated js code
       nglview enable # enable again, (not sure if needed)

You need to install ``nodejs`` (which includes ``npm``). Tips:
``conda install nodejs -c conda-forge`` (and so on)

Test notebook
=============

-  [edit to add more notebooks] and update notebook files

   .. code:: bash

       python ./devtools/make_test_js.py --api

-  install chromedriver:
   https://chromedriver.storage.googleapis.com/index.html?path=2.30/

   Download, unzip and copy chromedriver to /use/local/bin or anywhere
   in your PATH (tested on MacOS 10.12.5)

-  install nightwatch

::

    npm install -g nightwatch

-  install notebook runner

::

    source devtools/travis-ci/clone_nbtest.sh # only once

-  (may be):

To avoid entering notebook token or password, you might want to update

::

    c.NotebookApp.token = '' in $HOME/.jupyter/jupyter_notebook_config.py

-  Run notebook server

::

    jupyter notebook --port=8889 &

-  Run all tests

   .. code:: bash

       nightwatch

-  Run a single test

   .. code:: bash

       # nightwatch /path/to/your/file.js
       nightwatch nglview/tests/js/render_image.js

More stuff
==========

`wiki <https://github.com/arose/nglview/wiki>`__

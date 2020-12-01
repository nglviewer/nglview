`Installation <#installation>`__ \| `Example <#example>`__ \|
`Usage <#usage>`__ \| `Command line <#command-line>`__ \| `API
doc <#api-doc>`__ \| `Interface classes <interface_classes.html>`__ \|
`Website <#website>`__ \| `GUI <#show-gui>`__ \|
`Acknowledgment <#acknowledgment>`__

|Binder| |DOI| |Build Status| |bioconda-badge| |Coverage Status|

An `IPython/Jupyter <http://jupyter.org/>`__ widget to interactively
view molecular structures and trajectories. Utilizes the embeddable `NGL
Viewer <https://github.com/arose/ngl>`__ for rendering. Support for
showing data from the file-system, `RCSB PDB <http:www.rcsb.org>`__,
`simpletraj <https://github.com/arose/simpletraj>`__ and from objects of
analysis libraries `mdtraj <http://mdtraj.org/>`__,
`pytraj <http://amber-md.github.io/pytraj/latest/index.html>`__,
`mdanalysis <http://www.mdanalysis.org/>`__,
`ParmEd <http://parmed.github.io/ParmEd/>`__,
`rdkit <https://github.com/rdkit/rdkit>`__,
`ase <https://wiki.fysik.dtu.dk/ase/>`__,
`HTMD <https://www.htmd.org>`__,
`biopython <https://github.com/biopython/biopython.github.io/>`__,
`cctbx <https://cci.lbl.gov/cctbx_docs/iotbx/>`__,
`pyrosetta <http://pyrosetta.org>`__, `schrodinger's
Structure <http://content.schrodinger.com/Docs/r2015-4/python_api/api/schrodinger.structure.Structure-class.html>`__

Should work with Python 2 and 3. If you experience problems, please file
an `issue <https://github.com/arose/nglview/issues>`__.

Ask question about usage? Please post
`here <https://github.com/arose/nglview/issues/589>`__

.. figure:: https://github.com/arose/nglview/blob/master/examples/images/membrane.gif?raw=true
   :alt: membrane

   membrane

Table of contents
=================

-  `Installation <#installation>`__
-  `Example <#example>`__
-  `Showcase from users <#showcase-from-users>`__
-  `Usage <#usage>`__
-  `Contributing <#contributing>`__
-  `Command line <#command-line>`__
-  `API doc <#api-doc>`__
-  `Interface classes <interface_classes.html>`__
-  `Changelog <changelog.html>`__
-  `FAQ <#faq>`__
-  `Website <#website>`__
-  `Acknowledgment <#acknowledgment>`__
-  `License <#license>`__

Installation
============

Released version
----------------

-  Available on ``bioconda`` channel

   .. code:: bash

       conda install nglview -c bioconda
       # might need: jupyter-nbextension enable nglview --py --sys-prefix

       # if you already installed nglview, you can `upgrade`
       conda upgrade nglview --force
       # might need: jupyter-nbextension enable nglview --py --sys-prefix

-  Available on `PyPI <https://pypi.python.org/pypi/nglview/>`__

.. code:: bash

       pip install ipywidgets==5.2.2 # if you don't have ipywidgets yet
       pip install nglview
       # might need: jupyter-nbextension enable nglview --py --sys-prefix

Version Compatibility
---------------------

+-----------+--------------+
| nglview   | ipywidgets   |
+===========+==============+
| < 1.0     | 5.2.2        |
+-----------+--------------+
| 1.0       | 7.0          |
+-----------+--------------+

Notes
-----

If you are using ``notebook`` v5.0, you need to increase the
``iopub_data_rate_limit`` to `visualize big structure (e.g: solvated
system) <https://github.com/arose/nglview/issues/633>`__

::

    jupyter notebook --NotebookApp.iopub_data_rate_limit=10000000

Development version
-------------------

Requirement: ``ipywidgets >= 7.0.0``, ``notebook >= 4.2``

The development version can be installed directly from github:

notebook user
~~~~~~~~~~~~~

.. code:: bash

        git clone https://github.com/arose/nglview
        cd nglview
        python setup.py install
        
        # if you edit files in ./js folder, make sure to rebuild the code
        cd js
        npm install

        # probably need to activate widgetsnbextension
        # python -m ipykernel install --sys-prefix
        # jupyter nbextension enable --py --sys-prefix widgetsnbextension
        # jupyter nbextension enable --py --sys-prefix nglview
        
        # tested with ipywidgets 5.2.2, notebook 4.2.1

jupyterlab user
~~~~~~~~~~~~~~~

HIGHLY EXPERIMENTAL. Things might not work as you are expecting.

::

    # Please read its content before running
    source devtools/nglview-jupyterlab.sh

Example
=======

-  Notebooks: please see our `Jupyter notebook
   examples <https://github.com/arose/nglview/blob/master/examples/README.md>`__
-  Simple demo for trajectory (take time to load):
   `biomembrane <http://amber-md.github.io/pytraj/latest/ngl_player.html>`__

Showcase from users
===================

Please check `user examples <examples/user_examples.md>`__. Feel free to
contribute.

Usage
=====

Open a notebook

::

    jupyter notebook

and issue

.. code:: python

    import nglview
    view = nglview.show_pdbid("3pqr")  # load "3pqr" from RCSB PDB and display viewer widget
    view

A number of convenience functions are available to quickly display data
from the file-system, `RCSB PDB <http:www.rcsb.org>`__,
`simpletraj <https://github.com/arose/simpletraj>`__ and from objects of
analysis libraries `mdtraj <http://mdtraj.org/>`__,
`pytraj <http://amber-md.github.io/pytraj/latest/index.html>`__,
`mdanalysis <http://www.mdanalysis.org/>`__,
`ParmEd <http://parmed.github.io/ParmEd/>`__,
`rdkit <https://github.com/rdkit/rdkit>`__,
`HTMD <https://github.com/Acellera/htmd>`__,
`biopython <https://github.com/biopython/biopython.github.io/>`__.

+---------------------------------+------------------------------------------+
| Function                        | Description                              |
+=================================+==========================================+
| ``show_file(path)``             | Shows any NGL supported file formats     |
|                                 | (pdb, gro, mol2, sdf, dx, ..) in         |
|                                 | ``path``                                 |
+---------------------------------+------------------------------------------+
| ``show_pdbid(pdbid)``           | Shows ``pdbid`` fetched from RCSB PDB    |
+---------------------------------+------------------------------------------+
| ``show_simpletraj(struc_path, t | Shows structure & trajectory loaded with |
| raj_path)``                     | ``simpletraj``                           |
+---------------------------------+------------------------------------------+
| ``show_mdtraj(traj)``           | Shows ``MDTraj`` trajectory ``traj``     |
+---------------------------------+------------------------------------------+
| ``show_pytraj(traj)``           | Shows ``PyTraj`` trajectory ``traj``     |
+---------------------------------+------------------------------------------+
| ``show_parmed(structure)``      | Shows ``ParmEd`` structure               |
+---------------------------------+------------------------------------------+
| ``show_mdanalysis(univ)``       | Shows ``MDAnalysis`` Universe or         |
|                                 | AtomGroup ``univ``                       |
+---------------------------------+------------------------------------------+
| ``show_rdkit(mol)``             | Shows ``rdkit`` rdkit.Chem.rdchem.Mol    |
+---------------------------------+------------------------------------------+
| ``show_ase(atoms)``             | Shows ``ase`` Atoms                      |
+---------------------------------+------------------------------------------+
| ``show_asetraj(traj)``          | Shows ``ase`` trajectory ``traj``        |
+---------------------------------+------------------------------------------+
| ``show_htmd(mol)``              | Shows ``HTMD`` Molecules                 |
+---------------------------------+------------------------------------------+
| ``show_biopython(mol)``         | Shows ``Biopython`` structural entities  |
+---------------------------------+------------------------------------------+
| ``show_iotbx(mol)``             | Shows ``cctbx's iotbx`` structure        |
+---------------------------------+------------------------------------------+
| ``show_rosetta(pose)``          | Shows ``pyrosetta's Pose``               |
+---------------------------------+------------------------------------------+

API
===

Representations
---------------

.. code:: python

    view.add_representation(repr_type='cartoon', selection='protein')

    # or shorter
    view.add_cartoon(selection="protein")
    view.add_surface(selection="protein", opacity=0.3)

    # specify color
    view.add_cartoon(selection="protein", color='blue')

    # specify residue
    view.add_licorice('ALA, GLU')

    # clear representations
    view.clear_representations()

    # update parameters for ALL cartoons of component 0 (default)
    view.update_cartoon(opacity=0.4, component=0)

    # remove ALL cartoons of component 0 (default)
    view.remove_cartoon(opacity=0.4, component=0)

And many more, please check `NGL
website <http://nglviewer.org/ngl/api/index.html>`__

Representations can also be changed by overwriting the
``representations`` property of the widget instance ``view``. The
available ``type`` and ``params`` are described in the NGL Viewer
`documentation <http://nglviewer.org/ngl/api/index.html>`__.

.. code:: python

    view.representations = [
        {"type": "cartoon", "params": {
            "sele": "protein", "color": "residueindex"
        }},
        {"type": "ball+stick", "params": {
            "sele": "hetero"
        }}
    ]

The widget constructor also accepts a ``representation`` argument:

.. code:: python

    initial_repr = [
        {"type": "cartoon", "params": {
            "sele": "protein", "color": "sstruc"
        }}
    ]

    view = nglview.NGLWidget(struc, representation=initial_repr)
    view

Properties
----------

.. code:: python

    # set the frame number
    view.frame = 100

.. code:: python

    # parameters for the NGL stage object
    view.parameters = {
        # "percentages, "dist" is distance too camera in Angstrom
        "clipNear": 0, "clipFar": 100, "clipDist": 10,
        # percentages, start of fog and where on full effect
        "fogNear": 0, "fogFar": 100,
        # background color
        "backgroundColor": "black",
    }

    # note: NGLView accepts both origin camel NGL keywords (e.g. "clipNear")
    # and snake keywords (e.g "clip_near")

.. code:: python

    # parameters to control the `delay` between snapshots
    # change `step` to play forward (positive value) or backward (negative value)
    # note: experimental code
    view.player.parameters = dict(delay=0.04, step=-1)

.. code:: python

    # update camera type
    view.camera = 'orthographic'

.. code:: python

    # change background color
    view.background = 'black'

Trajectory
----------

.. code:: python

    # adding new trajectory
    view.add_trajectory(traj)
    # traj could be a `pytraj.Trajectory`, `mdtraj.Trajectory`, `MDAnalysis.Universe`, 
    # `parmed.Structure`, `htmd.Molecule` or derived class of `nglview.Trajectory`

    # change representation
    view.trajectory_0.add_cartoon(...) # equal to view.add_cartoon(component=0)
    view.trajectory_1.add_licorice(...) # equal to view.add_licorice(component=1)

Add extra component
-------------------

.. code:: python

    # Density volumes (MRC/MAP/CCP4, DX/DXBIN, CUBE)
    # Or adding derived class of `nglview.Structure`
    view.add_component('my.ccp4')

    # add component from url
    view.add_component('rcsb://1tsu.pdb')
    # NOTE: Trajectory is a special case of component.

Mouse
-----

.. code:: python

    # coot mouse style (https://en.wikipedia.org/wiki/Coot_(software))
    view.stage.set_parameters(mouse_preset='coot')

Interaction controls
--------------------

-  `Mouse <https://github.com/arose/ngl/blob/master/doc/usage/interaction-controls.md#mouse>`__
-  `Keyboard <https://github.com/arose/ngl/blob/master/doc/usage/interaction-controls.md#keyboard>`__

Display more than two widgets
-----------------------------

.. code:: python

    # 1st cell
    import ipywidgets
    vbox = ipywidgets.VBox([view1, view2])
    vbox # display

    # 2nd cell
    view1.sync_view()
    view2.sync_view()

Show GUI
--------

Notes: Unstable feature. `See
also <https://github.com/arose/nglview/blob/master/examples/README.md#unstable-features>`__

.. figure:: https://github.com/arose/nglview/blob/master/examples/images/nglview_gui.png?raw=true
   :alt: 

Movie making
------------

Notes: Unstable feature.

.. code:: python

    from nglview.contrib.movie import MovieMaker
    movie = MovieMaker(view, output='my.gif')
    movie.make()

API doc
=======

-  `Latest version <http://nglviewer.org/nglview/latest/api.html>`__
-  `All releases
   versions <http://nglviewer.org/nglview/release/index.html>`__
-  `Development version <http://nglviewer.org/nglview/dev/api.html>`__

Command line
============

.. code:: bash


    # open a notebook and import nglview
    nglview 

    # Require installing pytraj (PR for other backends is welcome)
    # open notebook, load `my.pdb` to pytraj's trajectory then display `view`
    nglview my.pdb

    # load density data
    nglview my.ccp4

    # open notebook, create trajectory with given topology `my.parm7` and trajecotry file `traj.nc`,
    # then display `view`
    nglview my.parm7 -c traj.nc

    # load all trajectories with filename ending with 'nc'
    # make sure to use quote " "
    nglview my.parm7 -c "*.nc"

    # open notebook, copy content from `myscript.py` then execute it
    nglview myscript.py

    # open notebook and execute 1st cell
    nglview mynotebook.ipynb

    # create a remote notebook
    # just follow its instruction
    nglview my.pdb --remote
    nglview my.parm7 -c traj.nc --remote
    nglview mynotebook.ipynb --remote

    # demo (don't need pytraj)
    nglview demo

    # disable autorun the 1st cell of the notebook
    nglview my.pdb --disable-autorun

    # specify web browser
    nglview my.pdb --browser=google-chrome

FAQ
===

`Q&A <https://github.com/arose/nglview/wiki/Q&A>`__

Website
=======

-  http://nglviewer.org/nglview/latest
-  http://nglviewer.org/nglview/dev

Talks
=====

`Talks about NGL and nglview <./talks.html>`__

Contributing
============

`Join us here <./contributing.html>`__

Projects using NGLView
======================

(Feel free to make a PR to add/remove your project here)

-  `AMBER <http://ambermd.org/>`__ - A package of programs for molecular
   dynamics simulations of proteins and nucleic acids
-  `mbuild <https://github.com/iModels/mbuild>`__ - A hierarchical,
   component based molecule builder
-  `deepchem <https://github.com/deepchem/deepchem>`__ - Deep-learning
   models for Drug Discovery and Quantum Chemistry
-  `htmd <https://github.com/Acellera/htmd>`__ - High throughput
   molecular dynamics simulations
-  `Moleidoscope <https://github.com/kbsezginel/Moleidoscope>`__ -
   Molecular kaleidoscope
-  `ssbio <https://github.com/nmih/ssbio>`__ - Tools for enabling
   structural systems biology
-  `hublib <https://github.com/martin-hunt/hublib>`__ - hublib is a
   Python library for the `HUBzero <https://hubzero.org/>`__ science
   gateway platform.
-  `molPX <https://github.com/markovmodel/molPX>`__: ipython API to
   visualize MD-trajectories along projected trajectories
-  `nanoribbon <https://github.com/oschuett/nanoribbon>`__
-  `ase <https://github.com/rosswhitfield/ase>`__: Atomic Simulation
   Environment
-  `pida <https://github.com/jharman25/pida>`__: Software for analyzing
   multiple protein-protein interaction docking solutions,

Acknowledgment
==============

-  Funding: Hai Nguyen is supported by NIH Grant GM103297, "The Center
   for HIV RNA Studies" (2015 to 02-2017).
-  Many thanks to ``nglview``
   `contributors <https://github.com/arose/nglview/graphs/contributors>`__
-  `dunovank/jupyter-themes <https://github.com/dunovank/jupyter-themes>`__:
   for ``oceans16`` theme
-  `base64-arraybuffer <https://github.com/niklasvh/base64-arraybuffer>`__
-  `ipywidgets <https://github.com/jupyter-widgets/ipywidgets>`__

License
=======

Generally MIT, see the LICENSE file for details.

.. |Binder| image:: http://mybinder.org/assets/images/logo.svg
   :target: http://mybinder.org/repo/hainm/nglview-notebooks
.. |DOI| image:: https://zenodo.org/badge/11846/arose/nglview.svg
   :target: https://zenodo.org/badge/latestdoi/11846/arose/nglview
.. |Build Status| image:: https://travis-ci.org/arose/nglview.svg?branch=master
   :target: https://travis-ci.org/arose/nglview
.. |bioconda-badge| image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square
   :target: http://bioconda.github.io
.. |Coverage Status| image:: https://coveralls.io/repos/github/arose/nglview/badge.png?branch=master
   :target: https://coveralls.io/github/arose/nglview

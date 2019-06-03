[Citation](#cite) | [Installation](#installation) | [Example](#example) | [Usage](#usage) | [Command line](#command-line) | [API doc](#api-doc) | [Interface classes](docs/interface_classes.md) | [Website](#website) | [Acknowledgment](#acknowledgment)

[![DOI](https://zenodo.org/badge/11846/arose/nglview.svg)](https://zenodo.org/badge/latestdoi/11846/arose/nglview)
[![Build Status](https://travis-ci.org/arose/nglview.svg?branch=master)](https://travis-ci.org/arose/nglview)
[![Coverage Status](https://coveralls.io/repos/github/arose/nglview/badge.png?branch=master)](https://coveralls.io/github/arose/nglview)

## Try nglview online: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/arose/nglview/master?urlpath=lab)


An [IPython/Jupyter](http://jupyter.org/) widget to interactively view molecular structures and trajectories. Utilizes the embeddable [NGL Viewer](https://github.com/arose/ngl) for rendering. Support for showing data from the file-system, [RCSB PDB](http:www.rcsb.org), [simpletraj](https://github.com/arose/simpletraj) and from objects of analysis libraries [mdtraj](http://mdtraj.org/), [pytraj](http://amber-md.github.io/pytraj/latest/index.html), [mdanalysis](http://www.mdanalysis.org/), [ParmEd](http://parmed.github.io/ParmEd/), [rdkit](https://github.com/rdkit/rdkit), [ase](https://wiki.fysik.dtu.dk/ase/), [HTMD](https://www.htmd.org), [biopython](https://github.com/biopython/biopython.github.io/), [cctbx](https://cci.lbl.gov/cctbx_docs/iotbx/), [pyrosetta](http://pyrosetta.org), [schrodinger's Structure](http://content.schrodinger.com/Docs/r2015-4/python_api/api/schrodinger.structure.Structure-class.html)

Should work with Python 3. If you experience problems, please file an [issue](https://github.com/arose/nglview/issues).

Ask question about usage? Please post [here](https://github.com/arose/nglview/issues/785)

![membrane](examples/images/membrane.gif)

Table of contents
=================

* [Installation](#installation)
* [Example](#example)
* [Showcase from users](#showcase-from-users)
* [Usage](#usage)
* [Contributing](#contributing)
* [Command line](#command-line)
* [API doc](#api-doc)
* [Interface classes](docs/interface_classes.md)
* [Changelog](CHANGELOG.md)
* [FAQ](#faq)
* [Website](#website)
* [Acknowledgment](#acknowledgment)
* [Cite](#cite)
* [License](#license)


Installation
============

Released version
----------------

- Available on `conda-forge` channel

    ```bash
    conda install nglview -c conda-forge
    # might need: jupyter-nbextension enable nglview --py --sys-prefix

    # if you already installed nglview, you can `upgrade`
    conda upgrade nglview --force
    # might need: jupyter-nbextension enable nglview --py --sys-prefix
    ```

- Available on [PyPI](https://pypi.python.org/pypi/nglview/)

```bash
   pip install nglview
   jupyter-nbextension enable nglview --py --sys-prefix
```

## Version Compatibility

| nglview | ipywidgets | ipykernel |
| --------|:----------:|----------:|
| < 1.0   | 5.2.2      | N/A
| 1.0     | 7.0.0      | 4.6

## Notes

If you are using `notebook` v5.0, you need to increase the `iopub_data_rate_limit`
to [visualize big structure (e.g: solvated system)](https://github.com/arose/nglview/issues/633)

```
jupyter notebook --NotebookApp.iopub_data_rate_limit=10000000
```

Development version
-------------------

Requirement: `ipywidgets >= 7.0`, `notebook >= 4.2`

The development version can be installed directly from github:

### notebook user

```bash
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
```
### jupyterlab user

Not supported yet.

Example
=======

- Notebooks: please see our [Jupyter notebook examples](examples/README.md)
- Simple demo for trajectory (take time to load): [biomembrane](http://amber-md.github.io/pytraj/latest/ngl_player.html)

Showcase from users
===================
Please check [user examples](examples/user_examples.md). Feel free to contribute.

Usage
=====

Open a notebook

    jupyter notebook

and issue

```Python
import nglview
view = nglview.show_pdbid("3pqr")  # load "3pqr" from RCSB PDB and display viewer widget
view
```

A number of convenience functions are available to quickly display data from
the file-system, [RCSB PDB](http:www.rcsb.org), [simpletraj](https://github.com/arose/simpletraj) and from objects of analysis libraries [mdtraj](http://mdtraj.org/), [pytraj](http://amber-md.github.io/pytraj/latest/index.html), [mdanalysis](http://www.mdanalysis.org/), [ParmEd](http://parmed.github.io/ParmEd/), [rdkit](https://github.com/rdkit/rdkit), [HTMD](https://github.com/Acellera/htmd), [biopython](https://github.com/biopython/biopython.github.io/).

| Function                                 | Description                                           |
|------------------------------------------|-------------------------------------------------------|
| `show_file(path)`                        | Shows any NGL supported file formats (pdb, gro, mol2, sdf, dx, ..) in `path`       |
| `show_pdbid(pdbid)`                      | Shows `pdbid` fetched from RCSB PDB                   |
| `show_simpletraj(struc_path, traj_path)` | Shows structure & trajectory loaded with `simpletraj` |
| `show_mdtraj(traj)`                      | Shows `MDTraj` trajectory `traj`                      |
| `show_pytraj(traj)`                      | Shows `PyTraj` trajectory `traj`                      |
| `show_parmed(structure)`                 | Shows `ParmEd` structure                              |
| `show_mdanalysis(univ)`                  | Shows `MDAnalysis` Universe or AtomGroup `univ`       |
| `show_rdkit(mol)`                        | Shows `rdkit` rdkit.Chem.rdchem.Mol                   |
| `show_ase(atoms)`                        | Shows `ase` Atoms                                     |
| `show_asetraj(traj)`                     | Shows `ase` trajectory `traj`                         |
| `show_pymatgen(struct)`                  | Shows `pymatgen` Structure                            |
| `show_htmd(mol)`                         | Shows `HTMD` Molecules                                |
| `show_biopython(mol)`                    | Shows `Biopython` structural entities                 |
| `show_iotbx(mol)`                        | Shows `cctbx's iotbx` structure                       |
| `show_rosetta(pose)`                     | Shows `pyrosetta's Pose`                              |
| `show_iodata(obj)`                       | Shows `iodata's IOData`                               |
| `show_psi4(obj)`                         | Shows `psi4's Molecule`                               |
| `show_qcelementary`                      | Shows `QCelementary's Molecule`                       |
| `show_openbabel`                         | Shows `openbabel's OMol`                              |
| `show_prody`                             | Shows `prody's Ensemble or AtomGroup`                 |


API
===

Representations
---------------

```python
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
```

And many more, please check [NGL website](http://arose.github.io/ngl/api/index.html)

Representations can also be changed by overwriting the `representations` property
of the widget instance `view`. The available `type` and `params` are described
in the NGL Viewer [documentation](http://arose.github.io/ngl/api/index.html).

```Python
view.representations = [
    {"type": "cartoon", "params": {
        "sele": "protein", "color": "residueindex"
    }},
    {"type": "ball+stick", "params": {
        "sele": "hetero"
    }}
]
```

The widget constructor also accepts a `representation` argument:

```Python
initial_repr = [
    {"type": "cartoon", "params": {
        "sele": "protein", "color": "sstruc"
    }}
]

view = nglview.NGLWidget(struc, representation=initial_repr)
view
```


Properties
----------

```Python
# set the frame number
view.frame = 100
```

```Python
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
```

```python
# parameters to control the `delay` between snapshots
# change `step` to play forward (positive value) or backward (negative value)
# note: experimental code
view.player.parameters = dict(delay=0.04, step=-1)
```

```python
# update camera type
view.camera = 'orthographic'
```

```python
# change background color
view.background = 'black'
```
Trajectory
----------

```python
# adding new trajectory
view.add_trajectory(traj)
# traj could be a `pytraj.Trajectory`, `mdtraj.Trajectory`, `MDAnalysis.Universe`, 
# `parmed.Structure`, `htmd.Molecule` or derived class of `nglview.Trajectory`

# change representation
view[0].add_cartoon(...) # equal to view.add_cartoon(component=0)
view[1].add_licorice(...) # equal to view.add_licorice(component=1)
```

Add extra component
-------------------

```python
# Density volumes (MRC/MAP/CCP4, DX/DXBIN, CUBE)
# Or adding derived class of `nglview.Structure`
view.add_component('my.ccp4')

# add component from url
view.add_component('rcsb://1tsu.pdb')
# NOTE: Trajectory is a special case of component.
```

Mouse
-----
```python
# coot mouse style (https://en.wikipedia.org/wiki/Coot_(software))
view.stage.set_parameters(mouse_preset='coot')
```

Interaction controls
--------------------
- [Mouse](https://github.com/arose/ngl/blob/master/docs/usage/interaction-controls.md#mouse)
- [Keyboard](https://github.com/arose/ngl/blob/master/docs/usage/interaction-controls.md#keyboard)


Movie making
------------

Notes: Unstable feature.

```python
from nglview.contrib.movie import MovieMaker
movie = MovieMaker(view, output='my.gif')
movie.make()
```

Embed widget
------------
[embed](docs/embed.md)

API doc
=======
- [Latest version](http://arose.github.io/nglview/latest/api.html)
- [All releases versions](http://arose.github.io/nglview/release/index.html)
- [Development version](http://arose.github.io/nglview/dev/api.html)

Command line
============

```bash

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
```

FAQ
===

[Q&A](docs/FAQ.md)

Website
=======

- http://arose.github.io/nglview/latest
- http://arose.github.io/nglview/dev

Talks
=====
[Talks about NGL and nglview](./talks.md)

Contributing
============

[Join us here](./CONTRIBUTING.md)

Projects integrating NGLView
============================

(Feel free to make a PR to add/remove your project here)

- [AMBER](http://ambermd.org/) -  A package of programs for molecular dynamics simulations of proteins and nucleic acids
- [mbuild](https://github.com/iModels/mbuild) - A hierarchical, component based molecule builder
- [deepchem](https://github.com/deepchem/deepchem) - Deep-learning models for Drug Discovery and Quantum Chemistry
- [htmd](https://github.com/Acellera/htmd) - High throughput molecular dynamics simulations
- [Moleidoscope](https://github.com/kbsezginel/Moleidoscope) - Molecular kaleidoscope
- [ssbio](https://github.com/nmih/ssbio) - Tools for enabling structural systems biology
- [hublib](https://github.com/martin-hunt/hublib) - hublib is a Python library for the [HUBzero](https://hubzero.org/) science gateway platform.
- [molPX](https://github.com/markovmodel/molPX): ipython API to visualize MD-trajectories along projected trajectories
- [nanoribbon](https://github.com/oschuett/nanoribbon)
- [ase](https://github.com/rosswhitfield/ase): Atomic Simulation Environment
- [pida](https://github.com/jharman25/pida): Software for analyzing multiple protein-protein interaction docking solutions,
- [pytim](https://github.com/Marcello-Sega/pytim) 
- [MobleyLab/drug-computing](https://github.com/MobleyLab/drug-computing) Educational materials for, and related to, UC Irvine's Drug Discovery Computing Techniques course.
- [pyiron](https://pyiron.github.io/): an integrated development environment for implementing, testing, and running simulations in computational materials science.
- [BioSimSpace](https://biosimspace.org/): An interoperable framework for biomolecular simulation
- [pyrod](https://github.com/schallerdavid/pyrod): PyRod - Tracing water molecules in molecular dynamics simulations

Acknowledgment
==============
- Funding: Hai Nguyen is supported by NIH Grant GM103297, "The Center for HIV RNA Studies" (2015 to 02-2017).
- Many thanks to `nglview` [contributors](https://github.com/arose/nglview/graphs/contributors)
- [dunovank/jupyter-themes](https://github.com/dunovank/jupyter-themes): for `oceans16` theme
- [base64-arraybuffer](https://github.com/niklasvh/base64-arraybuffer)
- [ipywidgets](https://github.com/jupyter-widgets/ipywidgets)

Cite
====
If you would like to acknowledge our work, feel free to cite:

Hai Nguyen, David A Case, Alexander S Rose; NGLview - Interactive molecular graphics for Jupyter notebooks, Bioinformatics, , btx789, https://doi.org/10.1093/bioinformatics/btx789

License
=======

Generally MIT, see the LICENSE file for details.

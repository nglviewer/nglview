
[![Binder](http://mybinder.org/images/logo.svg)](http://mybinder.org/repo/hainm/nglview-notebooks)
[![DOI](https://zenodo.org/badge/11846/arose/nglview.svg)](https://zenodo.org/badge/latestdoi/11846/arose/nglview)
[![Build Status](https://travis-ci.org/arose/nglview.svg?branch=master)](https://travis-ci.org/arose/nglview)

![nglview](nglview.png)


An [IPython/Jupyter](http://jupyter.org/) widget to interactively view molecular structures and trajectories. Utilizes the embeddable [NGL Viewer](https://github.com/arose/ngl) for rendering. Support for showing data from the file-system, [RCSB PDB](http:www.rcsb.org), [simpletraj](https://github.com/arose/simpletraj) and from objects of analysis libraries [mdtraj](http://mdtraj.org/), [pytraj](http://amber-md.github.io/pytraj/latest/index.html), [mdanalysis](http://www.mdanalysis.org/).

Should work with Python 2 and 3. If you experience problems, please file an [issue](https://github.com/arose/nglview/issues).

Table of contents
=================

* [Installation](#installation)
* [Usage](#usage)
* [Interface classes](#interface classes)
* [Changelog](CHANGELOG.md)
* [License](#license)


Installation
============

From PyPI:

    pip install nglview
Note: The above will try to install `jupyter`, `traitlets` and `ipywidgets`  as dependencies. If that fails install it manually `pip install jupyter`.

From Conda

    conda install -c omnia nglview


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
the file-system, [RCSB PDB](http:www.rcsb.org), [simpletraj](https://github.com/arose/simpletraj) and from objects of analysis libraries [mdtraj](http://mdtraj.org/), [pytraj](http://amber-md.github.io/pytraj/latest/index.html), [mdanalysis](http://www.mdanalysis.org/), [ParmEd](http://parmed.github.io/ParmEd/).

| Function                                 | Description                                           |
|------------------------------------------|-------------------------------------------------------|
| `show_structure_file(path)`              | Shows structure (pdb, gro, mol2, sdf) in `path`       |
| `show_pdbid(pdbid)`                      | Shows `pdbid` fetched from RCSB PDB                   |
| `show_simpletraj(struc_path, traj_path)` | Shows structure & trajectory loaded with `simpletraj` |
| `show_mdtraj(traj)`                      | Shows `MDTraj` trajectory `traj`                      |
| `show_pytraj(traj)`                      | Shows `PyTraj` trajectory `traj`                      |
| `show_parmed(structure)`                 | Shows `ParmEd` structure
| `show_mdanalysis(univ)`                  | Shows `MDAnalysis` Universe or AtomGroup `univ`       |


Representations
---------------

```python
view.add_cartoon("protein", color="residueindex")
view.add_surface("protein", opacity=0.3)
```

Representations can also be changed by overwriting the `representations` property
of the widget instance `view`. The available `type` and `params` are described
in the NGL Viewer [documentation](http://arose.github.io/ngl/doc).

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
nglview.NGLWidget(struc, representation=initial_repr)
```


Structures
----------

The above convenience functions first create an `adaptor` that implements an [interface](#Interface classes) for communication with the IPython/Jupyter widget.

```Python
import nglview
struc = nglview.PdbIdStructure("3pqr")  # load file from RCSB PDB
view = nglview.NGLWidget(struc)            # create widget
view                                       # display widget
```


Trajectories
------------

To enable trajectory access pass a second `Trajectory` argument to the widget
constructor or supply a combined `Structure`/`Trajectory` object as the first
argument.

Seperate `Structure` and `Trajectory` objects using `FileStructure` and
`SimpletrajStructure` (requires the [`simpletraj`](https://github.com/arose/simpletraj)
package):

```Python
import nglview
struc = nglview.FileStructure(nglview.datafiles.GRO)
traj = nglview.SimpletrajStructure(nglview.datafiles.XTC)
nglview.NGLWidget(struc, traj)
```

Combined `Structure`/`Trajectory` object utilizing `MDTrajTrajectory` which
wraps a trajectory loaded with [MDTraj](http://mdtraj.org/):

```Python
import nglview
import mdtraj
traj = mdtraj.load(nglview.datafiles.XTC, top=nglview.datafiles.GRO)
strucTraj = nglview.MDTrajTrajectory(traj)
nglview.NGLWidget(strucTraj)
```

The displayed frame can be changed by setting the `frame` property of the
widget instance `w`:

```Python
view.frame = 100  # set to frame no 100
```


Adaptors
--------

A number of adaptor classes are available to make structures and trajectories available to the widget.
They can support either the `Structure` (S) or the `Trajectory` (T) interface as well as both combined.

| Class                        | Description                                       | Interface |
|------------------------------|---------------------------------------------------|-----------|
| `FileStructure(path)`        | Loads `path` from filesystem                      | S         |
| `PdbIdStructure(pdbid)`      | Fetches `pdbid` from RCSB PDB                     | S         |
| `SimpletrajTrajectory(path)` | Uses `simpletraj` to access trajectory at `path`  | T         |
| `MDTrajTrajectory(traj)`     | Wraps `MDTraj` trajectory `traj`                  | S and T   |
| `PyTrajTrajectory(traj)`     | Wraps `PyTraj` trajectory `traj`                  | S and T   |
| `MDAnalysisTrajectory(univ)` | Wraps `MDAnalysis` Universe or AtomGroup `univ`   | S and T   |


Multiple widgets
----------------

You can have multiple widgets per notebook cell:

```Python
from ipywidgets.widgets import Box
w1 = NGLWidget(...)
w2 = NGLWidget(...)
Box(children=(w1,w2))
```


API
===

NGLWidget
---------

### Constructor

```Python
view = NGLWidget(structure, trajectory=None, representations=None)
```


### Properties

```Python
# set the frame number
view.frame = 100
```

```Python
# list of representations
view.representations = [{"type": "cartoon"}]
```

```Python
# parameters for the NGL stage object
view.parameters = {
    # "percentages, "dist" is distance too camera in Angstrom
    "clipNear": 0, "clipFar": 100, "clipDist": 10,
    # percentages, start of fog and where on full effect
    "fogNear": 0, "fogFar": 100
}
```

Interface classes
=================

You can create your own adaptors simply by following the interfaces for `Structure` and `Trajectory`, which can also be combined into a single class.


Structure
---------

```Python
class MyStructure(nglview.Structure):
    ext = "pdb"  # or gro, cif, mol2, sdf
    params = {}  # loading options passed to NGL
    def get_structure_string(self):
        return "structure in the self.ext format"
```


Trajectory
----------

```Python
class MyTrajectory(nglview.Trajectory):
    def get_coordinates(self, index):
        # return 2D numpy array, shape=(n_atoms, 3)

    def get_coordinates_dict(self):
        # return a dict of encoded 2D numpy array

    @property
    def n_frames(self):
        return 2  # return number of frames
```


Combined
--------

```Python
class MyStructureTrajectory(nglview.Structure, nglview.Trajectory):
    ext = "pdb"  # or gro, cif, mol2, sdf
    params = {}  # loading options passed to NGL

    def get_structure_string(self):
        return "structure in the self.ext format"

    def get_coordinates(self, index):
        # return 2D numpy array, shape=(n_atoms, 3)

    def get_coordinates_dict(self):
        # return a dict of encoded 2D numpy array

    @property
    def n_frames(self):
        return 2  # return number of frames
```

License
=======

Generally MIT, see the LICENSE file for details.

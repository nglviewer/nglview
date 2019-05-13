import uuid
"""abstract base class.
"""

__all__ = ['Structure', 'Trajectory']


class Structure:
    """abstract base class
    """

    def __init__(self):
        self.ext = "pdb"
        self.params = {}
        self.id = str(uuid.uuid4())

    def get_structure_string(self):
        raise NotImplementedError()


class Trajectory:
    """abstract base class
    """

    def __init__(self):
        self.id = str(uuid.uuid4())
        self.shown = True

    def get_coordinates(self, index):
        raise NotImplementedError()

    @property
    def n_frames(self):
        raise NotImplementedError()


__doc__ = """
Extend NGLView classes

Structures
==========

```python
import nglview
struc = nglview.PdbIdStructure("3pqr")     # load file from RCSB PDB
view = nglview.NGLWidget(struc)            # create widget
view                                       # display widget
```


Trajectories
============

To enable trajectory access pass a second `Trajectory` argument to the widget
constructor or supply a combined `Structure`/`Trajectory` object as the first
argument.

Seperate `Structure` and `Trajectory` objects using `FileStructure` and
`SimpletrajStructure` (requires the [`simpletraj`](https://github.com/arose/simpletraj)
package):

```python
import nglview
struc = nglview.FileStructure(nglview.datafiles.GRO)
traj = nglview.SimpletrajStructure(nglview.datafiles.XTC)
nglview.NGLWidget(struc, traj)
```

Combined `Structure`/`Trajectory` object utilizing `MDTrajTrajectory` which
wraps a trajectory loaded with [MDTraj](http://mdtraj.org/):

```python
import nglview
import mdtraj
traj = mdtraj.load(nglview.datafiles.XTC, top=nglview.datafiles.GRO)
strucTraj = nglview.MDTrajTrajectory(traj)
nglview.NGLWidget(strucTraj)
```

The displayed frame can be changed by setting the `frame` property of the
widget instance `w`:

```python
view.frame = 100  # set to frame no 100
```


Interface classes
=================

You can create your own adaptors simply by following the interfaces for `Structure` and `Trajectory`,
which can also be combined into a single class.


Structure
---------

```python
class MyStructure(nglview.Structure):
    ext = "pdb"  # or gro, cif, mol2, sdf
    params = {}  # loading options passed to NGL
    def get_structure_string(self):
        return "structure in the self.ext format"
```


Trajectory
----------

```python
class MyTrajectory(nglview.Trajectory):
    def get_coordinates(self, index):
        # return 2D numpy array, shape=(n_atoms, 3)

    @property
    def n_frames(self):
        return 2  # return number of frames
```


Combined
--------

```python
class MyStructureTrajectory(nglview.Structure, nglview.Trajectory):
    ext = "pdb"  # or gro, cif, mol2, sdf
    params = {}  # loading options passed to NGL

    def get_structure_string(self):
        return "structure in the self.ext format"

    def get_coordinates(self, index):
        # return 2D numpy array, shape=(n_atoms, 3)

    def n_frames(self):
        # return total frames
```
"""

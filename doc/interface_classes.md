
# Extend NGLView classes

[Structures](#structures) | [Trajectories](#trajectories) | [Interface class](#interface-classes) | [Register your backend](#register-your-backend)

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
ngl_traj = nglview.MDTrajTrajectory(traj)
view = nglview.NGLWidget(ngl_traj)
view
```

You can also add more trajectories to the widget
```python
ngl_traj2 = nglview.MDTrajTrajectory(traj2)
view.add_trajectory(ngl_traj2)
```

Interface classes
=================

You can create your own adaptors simply by following the interfaces for `Structure` and `Trajectory`, which can also be combined into a single class.


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

Register your backend
=====================
```python
from nglview import register_backend

@register_backend(your_package_name)
class NewTrajectoryClass:
    def __init__(your_traj, *args, **kwargs):
        # define your own implementation here
    ...

# if you already register your class, you can add `your_traj` directly to `view`
# without creating `NewTrajectoryClass` instance.
view.add_trajectory(your_traj)
```

Further reading: [nglview/adaptor.py](https://github.com/arose/nglview/blob/master/nglview/adaptor.py)

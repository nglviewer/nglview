
![nglview](nglview.png)


An IPython widget to interactively view molecular structures and trajectories. Utilizes the embeddable [NGL Viewer](https://github.com/arose/ngl) for rendering.

Very much work in progress. Please contact me if you want to take part. Should work with Python 2 and 3. If you experience problems, please file an [issue](https://github.com/arose/nglview/issues).


Table of contents
=================

* [Installation](#installation)
* [Usage](#Usage)
* [Interface classes](#Interface classes)
* [Changelog](#changelog)
* [License](#license)


Installation
============

From PyPI:

    pip install nglview


Usage
=====

Open a notebook

    ipython notebook

and issue

```Python
import nglview
struc = nglview.PdbIdStructure( "3pqr" )  # load file from RCSB PDB
w = nglview.NGLWidget( struc )            # create widget
w                                         # display widget
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
struc = nglview.FileStructure( "mystruc.pdb" )
traj = nglview.SimpletrajStructure( "mytraj.xtc" )
nglview.NGLWidget( struc, traj )
```

Combined `Structure`/`Trajectory` object utilizing `MDTrajTrajectory` which
wraps a trajectory loaded with [MDTraj](http://mdtraj.org/):

```Python
import nglview
import mdtraj
t = mdtraj.load( "mytraj.xtc", top="mystruc.gro" )
strucTraj = nglview.MDTrajTrajectory( t )
nglview.NGLWidget( strucTraj )
```


Representations
---------------

Representations can be changed by overwriting the `representations` property
of the widget instance `w`. The available `type` and `params` are described
in the NGL Viewer [documentation](http://arose.github.io/ngl/doc).

```Python
w.representations = [
    { "type": "cartoon", "params": {
        "sele": "protein", "color": "residueindex"
    } },
    { "type": "ball+stick", "params": {
        "sele": "hetero"
    } }
]
```

The widget constructor also accepts a `representation` argument:

```Python
initial_repr = [
    { "type": "cartoon", "params": {
        "sele": "protein", "color": "sstruc"
    } }
]
nglview.NGLWidget( struc, representation=initial_repr )
```


Adaptors
--------

A number of adaptor classes are available to make structures and trajectories available to the widget.
They can support either the `Structure` (S) or the `Trajectory` (T) interface as well as both combined.

| Class                          | Description                                       | Interface |
|--------------------------------|---------------------------------------------------|-----------|
| `FileStructure( path )`        | Loads `path` from filesystem                      | S         |
| `PdbIdStructure( pdbid )`      | Fetches `pdbid` from RCSB PDB                     | S         |
| `SimpletrajTrajectory( path )` | Uses `simpletraj` to access trajectory at `path`  | T         |
| `MDTrajTrajectory( t )`        | Wraps `MDTraj` trajectory `t`                     | S and T   |


Interface classes
=================

You can create your own adaptors simply by following the interfaces for `Structure` and `Trajectory`, which can also be combined into a single class.


Structure
---------

```Python
class MyStructure(nglview.Structure):
    ext = "pdb"  # or gro, cif, mol2, sdf
    def get_structure_string( self ):
        return "structure in the self.ext format"
```


Trajectory
----------

```Python
class MyTrajectory(nglview.Trajectory):
    def get_coordinates_list( self, index ):
        # return list of coordinates for frame at given index
        return [ x1, y1, z1, x2, y2, z2 ]
```


Combined
--------

```Python
class MyStructureTrajectory(nglview.Structure, nglview.Trajectory):
    ext = "pdb"  # or gro, cif, mol2, sdf
    def get_structure_string( self ):
        return "structure in the self.ext format"
    def get_coordinates_list( self, index ):
        # return list of coordinates for frame at given index
        return [ x1, y1, z1, x2, y2, z2 ]
```


Changelog
=========

Version 0.3dev
--------------

*


Version 0.2
-----------

* MIGRATION: changed `get_string` to `get_structure_string`
* MIGRATION: changed `get_coordinates` to `get_coordinates_list`
* DOC: usage, interface classes
* ADD: MDTrajTrajectory adaptor
* CODE: added interface classes
* CODE: suggested packages; mdtraj, simpletraj


Version 0.1
-----------

* initial version, no release


License
=======

Generally MIT, see the LICENSE file for details.

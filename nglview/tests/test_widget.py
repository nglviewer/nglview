# adpated from Jupyter ipywidgets project.

# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from __future__ import print_function

import nose.tools as nt
import unittest
from numpy.testing import assert_equal as eq, assert_almost_equal as aa_eq
import numpy as np

from ipykernel.comm import Comm
import ipywidgets as widgets

from traitlets import TraitError
from ipywidgets import Widget


import pytraj as pt
import nglview as nv
from nglview.utils import PY2, PY3


#-----------------------------------------------------------------------------
# Utility stuff from ipywidgets tests
#-----------------------------------------------------------------------------

class DummyComm(Comm):
    comm_id = 'a-b-c-d'

    def open(self, *args, **kwargs):
        pass

    def send(self, *args, **kwargs):
        pass

    def close(self, *args, **kwargs):
        pass

_widget_attrs = {}
displayed = []
undefined = object()

def setup():
    _widget_attrs['_comm_default'] = getattr(Widget, '_comm_default', undefined)
    Widget._comm_default = lambda self: DummyComm()
    _widget_attrs['_ipython_display_'] = Widget._ipython_display_
    def raise_not_implemented(*args, **kwargs):
        raise NotImplementedError()
    Widget._ipython_display_ = raise_not_implemented

def teardown():
    for attr, value in _widget_attrs.items():
        if value is undefined:
            delattr(Widget, attr)
        else:
            setattr(Widget, attr, value)

#-----------------------------------------------------------------------------
# NGLView stuff
#-----------------------------------------------------------------------------

def test_add_repr_shortcut():
    view = nv.show_pytraj(pt.datafiles.load_tz2())
    assert isinstance(view, nv.NGLWidget), 'must be instance of NGLWidget'
    view.add_cartoon(color='residueindex')
    view.add_rope(color='red')

def test_remote_call():
    # how to test JS?
    view = nv.show_pytraj(pt.datafiles.load_tz2())
    view._remote_call('centerView', target='stage')

    fn = 'notebooks/tz2.pdb'
    kwargs = {'defaultRepresentation': True}
    view._remote_call('loadFile', target='stage', args=[fn,], kwargs=kwargs)

# @unittest.skip("mess up with scipy, skip mdtraj now")
def test_show_mdtraj():
    import mdtraj as md
    from mdtraj.testing import get_fn
    fn = nv.datafiles.PDB 
    traj = md.load(fn)
    view = nv.show_mdtraj(traj)

@unittest.skipUnless(PY2, "only test MDAnalysis with PY2")
def test_show_MDAnalysis():
    from MDAnalysis import Universe
    tn, fn = nv.datafiles.PDB, nv.datafiles.PDB
    u = Universe(fn, tn)
    view = nv.show_mdanalysis(u)

def test_show_parmed():
    import parmed as pmd
    fn = nv.datafiles.PDB 
    parm = pmd.load_file(fn)
    view = nv.show_parmed(parm)

def test_caching_bool():
    view = nv.show_pytraj(pt.datafiles.load_tz2())
    nt.assert_false(view.cache)
    view.caching()
    nt.assert_true(view.cache)
    view.uncaching()
    nt.assert_false(view.cache)

def test_encode_and_decode():
    xyz = np.arange(100).astype('f4')
    shape = xyz.shape

    b64_str = nv.encode_numpy(xyz)
    new_xyz = nv.decode_base64(b64_str, dtype='f4', shape=shape)
    aa_eq(xyz, new_xyz) 

def test_coordinates_meta():
    traj = pt.datafiles.load_tz2()
    view = nv.show_pytraj(traj)
    view.frame = 3
    aa_eq(view.coordinates, traj.xyz[3])
    
    _coordinates_meta = view._coordinates_meta
    nt.assert_in('data', _coordinates_meta)
    nt.assert_in('shape', _coordinates_meta)
    nt.assert_in('dtype', _coordinates_meta)
    nt.assert_equal(view._coordinates_meta['dtype'], 'f4')

    view.coordinates = traj[2].xyz
    aa_eq(view.coordinates, traj.xyz[2])

    data = view._coordinates_meta['data']
    shape = view._coordinates_meta['shape']
    dtype = view._coordinates_meta['dtype']
    aa_eq(view.coordinates, nv.decode_base64(data, dtype=dtype, shape=shape))

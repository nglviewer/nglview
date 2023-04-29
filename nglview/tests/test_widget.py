import gzip
import os
import sys
import time
import unittest
from functools import partial
from io import StringIO
from itertools import chain

import ipywidgets
import numpy as np
import pytest
import traitlets
from ipykernel.comm import Comm
from IPython.display import display
from ipywidgets import BoundedFloatText, Button, HBox, IntText, Layout, Widget
from mock import MagicMock, patch
from numpy.testing import assert_almost_equal as aa_eq
from traitlets import TraitError, link

import nglview as nv
from nglview import NGLWidget, interpolate, js_utils, widget_utils
from nglview.representation import RepresentationControl
from nglview.utils.py_utils import click, decode_base64, encode_base64, submit
# local
from utils import get_fn
from utils import repr_dict as REPR_DICT

try:
    import simpletraj
    has_simpletraj = True
except ImportError:
    has_simpletraj = False

try:
    import pytraj as pt
    has_pytraj = True
except ImportError:
    pt = None
    has_pytraj = False

try:
    import mdtraj as md
    has_mdtraj = True
except ImportError:
    has_mdtraj = False
try:
    import parmed as pmd
    has_parmed = True
except ImportError:
    has_parmed = False

try:
    import MDAnalysis
    has_MDAnalysis = True
except ImportError:
    has_MDAnalysis = False

try:
    import htmd
    has_HTMD = True
except ImportError:
    has_HTMD = False

try:
    import ase
    has_ase = True
except ImportError:
    has_ase = False

try:
    import pymatgen
    has_pymatgen = True
except ImportError:
    has_pymatgen = False

try:
    import Bio.PDB
    has_bio = True
except ImportError:
    has_bio = False

try:
    import qcelemental
    has_qcelemental = True
except ImportError:
    has_qcelemental = False


def get_simple_traj():
    return nv.SimpletrajTrajectory(nv.datafiles.XTC, nv.datafiles.GRO)


def default_view():
    view = nv.NGLWidget()
    view.add_trajectory(get_simple_traj())
    return view


#-----------------------------------------------------------------------------
# NGLView stuff
#-----------------------------------------------------------------------------

DEFAULT_REPR = [{
    'params': {
        'sele': 'polymer'
    },
    'type': 'cartoon'
}, {
    'params': {
        'sele': 'hetero OR mol'
    },
    'type': 'ball+stick'
}, {
    "type": "ball+stick",
    "params": {
        "sele": "not protein and not nucleic"
    }
}]


def _assert_dict_list_equal(listdict0, listdict1):
    for (dict0, dict1) in zip(listdict0, listdict1):
        for (key0, key1) in zip(sorted(dict0.keys()), sorted(dict1.keys())):
            assert key0 == key1
            assert dict0.get(key0) == dict1.get(key1)


def test_API_promise_to_have():

    # for Jupyter notebook extension
    nv._jupyter_nbextension_paths()

    view = nv.demo()

    # trigger _set_size
    with patch.object(view, '_remote_call') as mock_call:
        view.layout.width = '100px'
        view.layout.height = '500px'
        mock_call.assert_called_with('setSize',
                                     args=['', '500px'],
                                     target='Widget')

    # Structure
    structure = nv.Structure()
    structure.get_structure_string
    assert hasattr(structure, 'id')
    assert hasattr(structure, 'ext')
    assert hasattr(structure, 'params')

    # Widget
    nv.NGLWidget._set_coordinates

    nv.NGLWidget.add_component
    nv.NGLWidget.add_trajectory
    nv.NGLWidget._coordinates_dict
    nv.NGLWidget.set_representations
    nv.NGLWidget.clear
    nv.NGLWidget.center

    # add component
    view.add_component('rcsb://1tsu.pdb')
    view.add_pdbid('1tsu')

    # display
    js_utils.clean_error_output()
    display(view.player.widget_repr)
    view.player._display()
    view._display_image()

    # show
    try:
        nv.show_pdbid('1tsu')
    except:
        pass
    nv.show_url('https://dummy.pdb')
    # other backends will be tested in other sections

    # constructor
    ngl_traj = get_simple_traj()
    nv.NGLWidget(ngl_traj, parameters=dict(background_color='black'))
    nv.NGLWidget(ngl_traj, representations=[dict(type='cartoon', params={})])

    view.parameters
    view.camera
    view.camera = 'perspective'
    view._request_stage_parameters()
    view._ngl_repr_dict = REPR_DICT
    view._handle_repr_dict_changed(dict(new=dict(c0={})))

    # dummy
    class DummWidget():
        value = ''

    view.player.picked_widget = DummWidget()

    view._update_background_color(change=dict(new='blue'))
    tab = view.player._display()

    view.player.widget_repr = view.player._make_widget_repr()
    view._handle_n_components_changed(change=dict(new=2, old=1))
    view._handle_n_components_changed(change=dict(new=1, old=1))
    view._handle_n_components_changed(change=dict(new=1, old=0))
    view.on_loaded(change=dict(new=True))
    view.on_loaded(change=dict(new=False))

    view._first_time_loaded = False
    view
    view._first_time_loaded = True
    view
    view._init_gui = True
    view
    view._theme = 'dark'
    view

    view.display(gui=True, style='ngl')
    view.display(gui=False)
    view.display(gui=True, style='ipywidgets')
    view._set_sync_camera([view])
    view._set_unsync_camera([view])
    view._set_selection('.CA')
    view.color_by('atomindex')
    representations = [dict(type='cartoon', params=dict())]
    view.representations = representations
    repr_parameters = dict(opacity=0.3, params=dict())
    view.update_representation(parameters=repr_parameters)
    view._remove_representation()
    view.clear()
    view.add_representation('surface', selection='*', useWorker=True)
    view.add_representation('surface', selection='*', component=1)
    view.center()
    view._on_render_image(change=dict(new='xyz'))
    view.render_image()
    view.render_image(frame=2)
    view.download_image()

    assert view._dry_run(view._set_sync_camera,
                         [view])['methodName'] == 'setSyncCamera'

    msg = dict(type='request_frame', data=dict())
    view._handle_custom_msg(msg=msg, buffers=[])
    msg = dict(type='repr_parameters', data=dict(name='hello'))
    view._handle_custom_msg(msg=msg, buffers=[])
    view.loaded = True
    msg = dict(type='request_loaded', data=True)
    view._handle_custom_msg(msg=msg, buffers=[])
    view.loaded = False
    msg = dict(type='request_loaded', data=True)
    view._handle_custom_msg(msg=msg, buffers=[])
    msg = dict(type='all_reprs_info', data=REPR_DICT)
    view._handle_custom_msg(msg=msg, buffers=[])
    msg = dict(type='stage_parameters', data=dict())
    view._handle_custom_msg(msg=msg, buffers=[])
    # test negative frame (it will be set to self.count - 1)
    view.frame = -1
    msg = dict(type='request_frame', data=dict())
    # async_message
    msg  = {'type': 'async_message', 'data': 'ok'}
    view._handle_custom_msg(msg, [])
    # render_image
    r = view.render_image()
    msg = {'type': 'image_data', 'ID': r.model_id, 'data': b'YmxhIGJsYQ=='}
    view._handle_custom_msg(msg, [])
    view.loaded = True
    view.show_only([
        0,
    ])
    view._js_console()
    view._get_full_params()

    # iter
    for c in view:
        assert isinstance(c, nv.widget.ComponentViewer)


@unittest.skipUnless(has_pytraj, 'skip if not having pytraj')
@unittest.skipUnless(has_mdtraj, 'skip if not having mdtraj')
def test_add_trajectory():
    view = nv.NGLWidget(default=False)

    def update_coords(view=view):
        view.frame = 1000
        view.frame = 0

    p_traj = pt.load(nv.datafiles.TRR, nv.datafiles.PDB)
    view.add_trajectory(p_traj)
    m_traj = md.load(nv.datafiles.XTC, top=nv.datafiles.PDB)
    view.add_trajectory(m_traj)
    # trigger updating coordinates
    update_coords()
    assert len(view._coordinates_dict.keys()) == 2
    if has_MDAnalysis:
        from MDAnalysis import Universe
        mda_traj = Universe(nv.datafiles.PDB, nv.datafiles.TRR)
        view.add_trajectory(mda_traj)
        update_coords()
        assert len(view._coordinates_dict.keys()) == 3
    if has_HTMD:
        from htmd import Molecule
        htmd_traj = Molecule(nv.datafiles.PDB)
        htmd_traj.filter('protein')
        view.add_trajectory(htmd_traj)
        update_coords()
        if has_MDAnalysis:
            assert len(view._coordinates_dict.keys()) == 4
        else:
            assert len(view._coordinates_dict.keys()) == 3


def test_API_promise_to_have_add_more_backend():
    @nv.register_backend('dummy')
    class MyLovelyClass(nv.Structure, nv.Trajectory):
        pass

    assert 'dummy' in nv.BACKENDS


def test_handling_n_components_changed():
    view = nv.NGLWidget()
    n_traj = get_simple_traj()
    view.add_trajectory(n_traj)
    # fake updating n_components and _repr_dict from front-end
    view._ngl_repr_dict = REPR_DICT
    view.n_components = 1
    view.player.widget_repr = view.player._make_widget_repr()
    view.remove_component(n_traj.id)
    # fake updating n_components from front-end
    view._ngl_repr_dict = {'c0': {}}
    view.n_components = 0


def test_base_adaptor():
    # abstract base class
    def func_0():
        nv.Structure().get_structure_string()

    def func_1():
        nv.Trajectory().get_coordinates(1)

    def func_2():
        nv.Trajectory().n_frames

    pytest.raises(NotImplementedError, func_0)
    pytest.raises(NotImplementedError, func_1)
    pytest.raises(NotImplementedError, func_2)


def test_coordinates_dict():
    traj = get_simple_traj()
    view = nv.NGLWidget(traj)
    view.frame = 1
    coords = view._coordinates_dict[0]
    aa_eq(coords, traj.get_coordinates(1))

    # dummy
    view._send_binary = False
    view._coordinates_dict = {0: coords}
    # increase coverage for IndexError: make index=1000 (which is larger than n_frames)
    view.player.interpolate = True
    view._set_coordinates(1000)


def test_load_data():
    view = default_view()

    # load blob with ext
    blob = open(nv.datafiles.PDB).read()
    view._load_data(blob, ext='pdb')

    # raise if passing blob but does not provide ext
    with pytest.raises(ValueError):
        view._load_data(blob)

    # raise if passing dummy name
    with pytest.raises(NameError):
        view._load_data(hahahaha)

    # load PyTrajectory
    t0 = get_simple_traj()
    view._load_data(t0)

    # load current folder
    view._load_data(get_fn('tz2.pdb'))


def test_representations():
    view = default_view()
    view.representations = DEFAULT_REPR
    view.add_cartoon()
    representations_2 = DEFAULT_REPR[:]
    representations_2.append({'type': 'cartoon', 'params': {'sele': 'all'}})
    _assert_dict_list_equal(view.representations, representations_2)

    # accept dict too (to specify seperate reprs for different component
    def func():
        view.representations = {'0': MagicMock()}
    assert view._dry_run(func)['methodName'] == '_set_representation_from_repr_dict'

    # Representations
    # make fake params
    try:
        view._ngl_repr_dict = {'c0': {'0': {'parameters': {}}}}
    except (KeyError, TraitError):
        # in real application, we are not allowed to assign values
        pass

    view._ngl_repr_dict = REPR_DICT
    representation_widget = RepresentationControl(view, 0, 0)
    representation_widget
    representation_widget._on_parameters_changed(change=dict(new=dict()))


def test_representation_control():
    view = nv.demo()
    repr_control = view._display_repr()

    repr_control.name = 'surface'
    repr_control.name = 'cartoon'
    repr_control.repr_index = 1
    repr_control.component_index = 1


def test_add_repr_shortcut():
    view = default_view()
    assert isinstance(view, nv.NGLWidget), 'must be instance of NGLWidget'

    # add
    view.add_cartoon(color='residueindex')
    view.add_rope(color='red')

    # update
    view.update_cartoon(opacity=0.4)
    view.update_rope(coor='blue')

    # remove
    view.remove_cartoon()
    view.remove_rope()


def test_color_scheme():
    view = nv.demo()
    scheme = nv.color._ColorScheme([['red', '1-6'], ['yellow', '20-30']],
                                   'what')
    view.clear()
    view.add_cartoon(color=scheme)


def test_add_new_shape():
    view = nv.NGLWidget()
    sphere = ('sphere', [0, 0, 9], [1, 0, 0], 1.5)
    arrow = ('arrow', [1, 2, 7], [30, 3, 3], [1, 0, 1], 1.0)
    c0 = view._add_shape([sphere, arrow], name='my_shape')

    # Shape
    c1 = view.shape.add_arrow([1, 2, 7], [30, 3, 3], [1, 0, 1], 1.0)
    assert len(view._ngl_component_ids) == 2
    view.remove_component(c0)
    assert len(view._ngl_component_ids) == 1

    view.remove_component(c1)
    assert len(view._ngl_component_ids) == 0


def test_add_buffer():
    view = nv.NGLWidget()
    view
    kwargs = {
        "position": [0, 0, 0, 1, 1, 1],
        "color": [1, 0, 0, 255, 0, 0],
        "radius": [1., 2.]
    }

    view.shape.add_buffer('sphere', **kwargs)


def test_remote_call():
    # how to test JS?
    view = default_view()
    view._remote_call('centerView', target='stage')

    fn = 'notebooks/tz2.pdb'
    kwargs = {'defaultRepresentation': True}
    view._remote_call('loadFile', target='stage', args=[
        fn,
    ], kwargs=kwargs)


def test_download_image():
    """just make sure it can be called
    """
    view = default_view()
    view.download_image('myname.png', 2, False, False, True)


def test_show_structure_file():
    view = nv.show_structure_file(nv.datafiles.PDB)


def test_show_file():
    view = nv.show_file(nv.datafiles.PDB)


def test_show_text():
    text = open(nv.datafiles.PDB).read()
    nv.show_text(text)


@unittest.skipUnless(has_ase, 'skip if not having ase')
def test_show_ase():
    from ase import Atom, Atoms
    dimer = Atoms([Atom('X', (0, 0, 0)), Atom('X', (0, 0, 1))])
    dimer.set_positions([(1, 2, 3), (4, 5, 6.2)])
    nv.show_ase(dimer)


@unittest.skipUnless(has_pymatgen, 'skip if not having pymatgen')
def test_show_pymatgen():
    from pymatgen.core import Lattice, Structure
    lattice = Lattice.cubic(4.2)
    structure = Structure(lattice, ["Cs", "Cl"],
                             [[0, 0, 0], [0.5, 0.5, 0.5]])
    view = nv.show_pymatgen(structure)
    view


@unittest.skipUnless(has_qcelemental, 'skip if not having qcelemental')
def test_show_qcelemental():
    import qcelemental as qcel

    mol = qcel.models.Molecule.from_data("He 0 0 0") 
    view = nv.show_qcelemental(mol)
    view


@unittest.skipUnless(has_bio, 'skip if not having biopython')
def test_show_biopython():
    from Bio.PDB import PDBParser
    parser = PDBParser()
    structure = parser.get_structure('protein', nv.datafiles.PDB)
    nv.show_biopython(structure)


@unittest.skipUnless(has_simpletraj, 'skip if not having simpletraj')
def test_show_simpletraj():
    traj = nv.SimpletrajTrajectory(nv.datafiles.XTC, nv.datafiles.GRO)
    view = nv.show_simpletraj(traj)
    view
    view.frame = 3


@unittest.skipUnless(has_mdtraj, 'skip if not having mdtraj')
def test_show_mdtraj():
    import mdtraj as md
    traj = md.load(nv.datafiles.PDB)
    view = nv.show_mdtraj(traj)


@unittest.skipUnless(has_HTMD, 'skip if not having HTMD')
def test_show_htmd():
    from htmd import Molecule
    fn = nv.datafiles.PDB
    traj = Molecule(fn)
    view = nv.show_htmd(traj)
    # trigger updating cooridnates
    view.frame = 100
    index = 0
    view.frame = index
    xyz_htmd = np.squeeze(traj.coords[:, :, index])
    aa_eq(view._coordinates_dict[0], xyz_htmd)


@unittest.skipUnless(has_MDAnalysis, 'skip if not having MDAnalysis')
def test_show_MDAnalysis():
    from MDAnalysis import Universe
    tn, fn = nv.datafiles.PDB, nv.datafiles.PDB
    u = Universe(fn, tn)
    view = nv.show_mdanalysis(u)


@unittest.skipUnless(has_parmed, 'skip if not having ParmEd')
def test_show_parmed():
    import parmed as pmd
    fn = nv.datafiles.PDB
    parm = pmd.load_file(fn)
    view = nv.show_parmed(parm)

    ngl_traj = nv.ParmEdTrajectory(parm)
    ngl_traj.only_save_1st_model = False
    ngl_traj.get_structure_string()


def test_encode_and_decode():
    xyz = np.arange(100).astype('f4')
    shape = xyz.shape

    b64_str = encode_base64(xyz)
    new_xyz = decode_base64(b64_str, dtype='f4', shape=shape)
    aa_eq(xyz, new_xyz)


def test_structure_file():
    for fn in [get_fn('tz2.pdb'), nv.datafiles.GRO]:
        content = open(fn, 'r').read()
        fs1 = nv.FileStructure(fn)
        assert content == fs1.get_structure_string()

    # gz
    fn = get_fn('tz2_2.pdb.gz')
    fs2 = nv.FileStructure(fn)
    content = gzip.open(fn).read()
    assert content == fs2.get_structure_string()


def test_camelize_parameters():
    view = nv.NGLWidget()
    view.parameters = dict(background_color='black')
    assert 'backgroundColor' in view._parameters


def test_component_for_duck_typing():
    # FIXME: deprecate duck typing?
    # syntax looks ugly.
    view = NGLWidget()
    traj = nv.SimpletrajTrajectory(nv.datafiles.XTC, nv.datafiles.GRO)

    # add 3 components (trajectory is a component)
    view.add_component(get_fn('tz2.pdb'))
    view.add_component(get_fn('tz2_2.pdb.gz'))
    view.add_trajectory(traj)
    view.component_0.add_representation('cartoon')

    c0 = view[0]
    c1 = view[1]
    assert hasattr(view, 'component_0')
    assert hasattr(view, 'component_1')
    assert hasattr(view, 'trajectory_0')
    assert hasattr(view.trajectory_0, 'n_frames')
    assert hasattr(view.trajectory_0, 'get_coordinates')
    assert hasattr(view.trajectory_0, 'get_structure_string')

    c0.show()
    c0.hide()

    # 2 components left
    view.remove_component(c0.id)
    # c1 become 1st component
    assert not hasattr(view, 'component_2')
    assert len(view._ngl_component_ids) == 2

    # negative indexing
    assert view[0]._index == c1._index


def test_trajectory_show_hide_sending_cooridnates():
    view = NGLWidget()

    traj0 = get_simple_traj()
    traj1 = get_simple_traj()

    view.add_trajectory(traj0)
    view.add_trajectory(traj1)

    for traj in view._trajlist:
        assert traj.shown

    view.frame = 1

    def copy_coordinate_dict(view):
        # make copy to avoid memory free
        return {k: v.copy() for k, v in view._coordinates_dict.items()}

    coordinates_dict = copy_coordinate_dict(view)
    aa_eq(coordinates_dict[0], traj0.get_coordinates(1))
    aa_eq(coordinates_dict[1], traj1.get_coordinates(1))

    # hide 0
    view.hide([
        0,
    ])
    assert not view._trajlist[0].shown
    assert view._trajlist[1].shown

    # update frame so view can update its coordinates
    view.frame = 2
    coordinates_dict = copy_coordinate_dict(view)
    assert coordinates_dict[0].shape[0] == 0
    aa_eq(coordinates_dict[1], traj1.get_coordinates(2))

    # hide 0, 1
    view.hide([0, 1])
    assert not view._trajlist[0].shown
    assert not view._trajlist[1].shown
    view.frame = 3
    coordinates_dict = copy_coordinate_dict(view)
    assert coordinates_dict[0].shape[0] == 0
    assert coordinates_dict[1].shape[0] == 0

    # slicing, show only component 1
    view[1].show()
    view.frame = 0
    assert not view._trajlist[0].shown
    assert view._trajlist[1].shown
    coordinates_dict = copy_coordinate_dict(view)
    assert coordinates_dict[0].shape[0] == 0
    aa_eq(coordinates_dict[1], traj1.get_coordinates(0))

    # show all
    view[1].show()
    view[0].show()
    view.show(indices='all')
    view.show(indices=[
        0,
    ])
    view.show(indices=[0, 1])
    view.frame = 1
    assert view._trajlist[1].shown
    coordinates_dict = copy_coordinate_dict(view)
    aa_eq(coordinates_dict[0], traj0.get_coordinates(1))
    aa_eq(coordinates_dict[1], traj1.get_coordinates(1))

    # hide all
    view[1].hide()
    view[0].hide()
    view.frame = 2
    assert not view._trajlist[0].shown
    assert not view._trajlist[1].shown
    coordinates_dict = copy_coordinate_dict(view)
    assert coordinates_dict[0].shape[0] == 0
    assert coordinates_dict[1].shape[0] == 0


def test_existing_js_files():
    from glob import glob
    jsfiles = glob(os.path.join(os.path.dirname(nv.__file__), 'static', '*js'))
    mapfiles = glob(
        os.path.join(os.path.dirname(nv.__file__), 'static', '*map'))

    assert len(jsfiles) == 2
    assert len(mapfiles) == 1


def test_add_structure():
    view = nv.NGLWidget()
    with pytest.raises(ValueError):
        # raise if not is instance of nv.Structure
        view.add_structure(nv.datafiles.PDB)


def test_add_struture_then_trajectory():
    view = nv.show_structure_file(get_fn('tz2.pdb'))
    view.loaded = True
    traj = get_simple_traj()
    view.add_trajectory(traj)
    view.frame = 3
    coords = view._coordinates_dict[1].copy()
    expected = traj.get_coordinates(3)
    aa_eq(coords, expected)
    view.loaded = False
    view.add_trajectory(traj)


def test_loaded_attribute():
    traj = get_simple_traj()
    structure = nv.FileStructure(nv.datafiles.PDB)

    # False, empty constructor
    view = nv.NGLWidget()
    view.loaded = False
    view.add_structure(structure)
    view.add_trajectory(traj)
    view

    # False, constructor with a single Structure
    view = nv.NGLWidget(structure)
    view.loaded = False
    view.add_trajectory(traj)
    view

    # True
    view = nv.NGLWidget()
    view.loaded = True
    view.add_structure(structure)
    view.add_trajectory(traj)
    view

    # False then True, empty constructor
    view = nv.NGLWidget()
    view.loaded = False
    view.add_structure(structure)
    view.loaded = True
    view.add_trajectory(traj)
    view

    # False then True, constructor with a Trajectory
    view = nv.NGLWidget(traj)
    view.loaded = False
    view.add_structure(structure)
    view.loaded = True
    view.add_trajectory(traj)
    view


def test_player_simple():
    view = default_view()

    # dummy
    component_slider = ipywidgets.IntSlider()
    repr_slider = ipywidgets.IntSlider()

    # dummy test
    player = nv.player.TrajectoryPlayer(view)
    player.smooth()
    player.camera = 'perspective'
    player.camera = 'orthographic'
    player.frame
    player.frame = 10
    player.parameters = dict(step=2)
    player._display()
    player._make_button_center()
    w = player._make_widget_preference()
    w.children[0].value = 1.
    player.widget_preference = None
    w = player._make_widget_preference()
    w.children[0].value = 1.
    player._show_download_image()
    player._make_text_picked()
    player._refresh(component_slider, repr_slider)
    player._make_widget_repr()
    player._make_resize_notebook_slider()
    player._make_button_export_image()
    player._make_repr_playground()
    player._make_widget_picked()
    player._make_export_image_widget()
    player._make_general_box()
    player._update_padding()
    player._real_time_update = True
    player._make_widget_repr()
    player.widget_component_slider
    player.widget_repr_slider
    player._create_all_tabs()
    player._create_all_widgets()
    player.widget_tab = None
    player._create_all_widgets()
    player._simplify_repr_control()

    player._real_time_update = True
    player.widget_repr_slider.value = 0
    player.widget_repr_slider.value = 1
    slider_notebook = player._make_resize_notebook_slider()
    slider_notebook.value = 300

    player.widget_repr_name.value = 'surface'
    player.widget_repr_name.value = 'cartoon'


def test_player_submit_text():
    """ test_player_click_button """
    view = nv.demo(gui=True)
    submit(view.player._make_command_box())


def test_player_click_button():
    """ test_player_click_button """
    view = nv.demo(gui=True)
    view
    view._ngl_repr_dict = REPR_DICT
    view.player._create_all_widgets()
    view.player.widget_export_image = view.player._make_button_export_image()
    button_iter = chain.from_iterable([
        view.player.widget_repr_control_buttons.children,
        [
            view.player._show_download_image(),
            view.player._make_button_center(),
            view.player.widget_export_image.children[0].children[0],
            view.player.widget_repr_add.children[0],
        ],
        [
            w for w in view.player.widget_preference.children
            if isinstance(w, Button)
        ],
    ])
    for button in button_iter:
        click(button)


def test_player_link_to_ipywidgets():
    view = default_view()

    int_text = IntText(2)
    float_text = BoundedFloatText(40, min=10)
    HBox([int_text, float_text])
    link((int_text, 'value'), (view.player, 'step'))
    link((float_text, 'value'), (view.player, 'delay'))

    assert view.player.step == 2
    assert view.player.delay == 40

    float_text.value = 100
    assert view.player.delay == 100

    float_text.value = 0.00
    # we set min=10
    assert view.player.delay == 10


def test_player_interpolation():
    view = default_view()

    view.player.interpolate = True
    assert view.player.iparams.get('type') == 'linear'
    assert view.player.iparams.get('step') == 1


def test_player_picked():
    view = nv.demo()
    s = dict(x=3)
    view.player.widget_picked = view.player._make_text_picked()
    view.picked = s
    assert view.player.widget_picked.value == '{"x": 3}'


def test_widget_utils():
    box = HBox()
    i0 = IntText()
    i0._ngl_name = 'i0'
    i1 = IntText()
    i1._ngl_name = 'i1'
    box.children = [i0, i1]

    assert i0 is widget_utils.get_widget_by_name(box, 'i0')
    assert i1 is widget_utils.get_widget_by_name(box, 'i1')

    box.children = [i1, i0]
    assert i0 is widget_utils.get_widget_by_name(box, 'i0')
    assert i1 is widget_utils.get_widget_by_name(box, 'i1')

    assert widget_utils.get_widget_by_name(box, 'i100') is None
    assert widget_utils.get_widget_by_name(None, 'i100') is None


def test_adaptor_raise():
    with pytest.raises(ValueError):
        nv.FileStructure('hellotheredda.pdb')


def test_theme():
    from nglview import theme
    # FIXME: fill me


def test_player_click_tab():
    view = nv.demo()
    gui = view.player._display()
    assert isinstance(gui, ipywidgets.Tab)

    for i, child in enumerate(gui.children):
        try:
            gui.selected_index = i
            assert isinstance(child, ipywidgets.Box)
        except TraitError:
            pass


def test_interpolate():
    # dummy test
    traj = get_simple_traj()
    interpolate.linear(0, 0.4, traj, step=1)


def dummy_test_to_increase_coverage():
    nv.__version__


def test_viewer_control():
    view = nv.demo()
    view

    mat = [11, 12, 13, 14, 21, 22, 23, 24, 31, 32, 33, 34, 41, 42, 43, 44]

    vector = [0, 1, 2]

    view.control.align(mat)
    view.control.rotate(mat)
    view.control.translate(vector)
    view.control.apply_matrix(mat)
    view.control.center(vector)
    view.control.orient(mat)
    view.control.zoom(0.3)
    view.control.rotate(mat)
    view.control.spin(vector, 0.1)


def test_queuing_messages():
    view = nv.NGLWidget()
    view.add_component(nv.datafiles.PDB)
    view.download_image()
    view
    assert [f._method_name for f in view._ngl_displayed_callbacks_before_loaded] == \
           [
            'loadFile',
            '_downloadImage']
    assert [f['methodName'] for f in view._ngl_msg_archive] == \
           ['loadFile']

    # display 2nd time
    view
    assert [f['methodName'] for f in view._ngl_msg_archive] == \
           ['loadFile']


@patch('nglview.NGLWidget._unset_serialization')
def test_write_html(mock_unset):
    from nglview.color import ColormakerRegistry as cm
    from nglview.theme import ThemeManager
    import ipywidgets.embed as embed

    tm = ThemeManager()
    traj0 = get_simple_traj()
    traj1 = get_simple_traj()
    view = nv.NGLWidget()
    view.add_trajectory(traj0)
    view.add_trajectory(traj1)
    view.gui_style = 'ngl'
    view._gui_theme = 'dark'
    display(view)
    fp = StringIO()

    with patch.object(embed, 'embed_snippet') as mock_embed:
        nv.write_html(fp, [view], frame_range=(0, 3))
        mock_embed.assert_called_with([tm, cm, view])
    mock_unset.assert_called_with()
    assert len(view._ngl_coordinate_resource[0]) == 3
    assert len(view._ngl_coordinate_resource[1]) == 3

    # box
    with patch.object(embed, 'embed_snippet') as mock_embed:
        nv.write_html(fp, [HBox([view])], frame_range=(0, 3))
        # FIXME: assertion?


def test_trim_messages():
    view = nv.demo()
    view.remove_component(view[0])
    assert view._ngl_msg_archive == []
    view.add_component(nv.datafiles.ALA3)
    assert len(view._ngl_msg_archive) == 1
    assert view._ngl_msg_archive[0]['methodName'] == 'loadFile'

    view = nv.demo()
    c = view.add_component(nv.datafiles.ALA3)
    view.remove_component(c)
    assert len(view._ngl_msg_archive) == 1
    assert view._ngl_msg_archive[0]['methodName'] == 'loadFile'


def test_fullscreen():
    v = nv.demo()
    fs = nv.widget.Fullscreen(v, [v])
    fs.fullscreen()
    with patch.object(v, 'handle_resize'):
        fs._fullscreen_changed(MagicMock(new=True, old=False))
        assert v.handle_resize.called
    # just run the code
    fs._fullscreen_changed(MagicMock(new=True, old=False))
    fs._fullscreen_changed(MagicMock(new=False, old=True))

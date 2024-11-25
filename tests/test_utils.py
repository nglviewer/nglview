import gzip
import os

import pytest

import nglview
from nglview.utils import js_utils, py_utils
from nglview.utils.py_utils import (FileManager, _camelize, _camelize_dict,
                                    seq_to_string)
# local
# local
from utils import get_fn
from utils import repr_dict as repr_dict_example


def assert_equal(x, y):
    assert x == y


def test_get_name():
    fn = nglview.datafiles.PDB
    assert_equal(py_utils.get_name(object, name='hello'), 'hello')
    assert_equal(py_utils.get_name(nglview.FileStructure(fn)),
                 'nglview.adaptor.FileStructure')


def test_seq_to_string():
    assert_equal(seq_to_string([1, 2, 3]), '@1,2,3')
    assert_equal(seq_to_string('@1,2,3'), '@1,2,3')


def test_camelize():
    assert_equal(_camelize('remote_call'), 'remoteCall')
    assert_equal(_camelize('flat_shaded'), 'flatShaded')
    assert_equal(_camelize('near_clip'), 'nearClip')
    assert_equal(_camelize('radius_type'), 'radiusType')
    assert_equal(_camelize('aspect_ratio'), 'aspectRatio')
    assert_equal(_camelize('radius_segments'), 'radiusSegments')
    assert_equal(_camelize('sphere_detail'), 'sphereDetail')
    assert_equal(_camelize('contact_type'), 'contactType')
    assert_equal(_camelize('max_distance'), 'maxDistance')
    assert_equal(_camelize('max_angle'), 'maxAngle')
    assert_equal(_camelize('label_type'), 'labelType')
    assert_equal(_camelize('radial_segments'), 'radialSegments')
    assert_equal(_camelize('surface_type'), 'surfaceType')
    assert_equal(_camelize('probe_radius'), 'probeRadius')
    assert_equal(_camelize('scale_factor'), 'scaleFactor')

    # make sure this method does not change the camel
    assert_equal(_camelize('maxAngle'), 'maxAngle')
    assert_equal(_camelize('maxDistance'), 'maxDistance')
    assert_equal(_camelize('color'), 'color')


def test_dict():
    kwargs = dict(default_representation=True)
    kwargs2 = _camelize_dict(kwargs)
    assert 'defaultRepresentation' in kwargs2


def test_file_manager_use_url():
    fh = FileManager('rcsb://1tsu.pdb')
    assert fh.is_url
    assert_equal(fh.ext, 'pdb')
    assert not fh.is_compressed

    fh = FileManager('http://dummy.com/1tsu.pdb')
    assert fh.is_url
    assert_equal(fh.ext, 'pdb')

    fh = FileManager('http://dummy.com/1tsu.pdb.gz')
    assert fh.is_url
    assert_equal(fh.ext, 'pdb')
    assert fh.is_compressed


def test_file_not_use_filename():
    src = os.path.join(os.path.dirname(nglview.__file__), '__init__.py')
    fh = FileManager(src)
    assert not fh.is_compressed

    assert fh.ext.endswith('py')
    assert fh.is_filename

    with open(src) as src2:
        fh2 = FileManager(src2)
        assert not fh2.is_filename
        fh2_content = fh2.read()
        src2.seek(0)
        assert fh2_content == src2.read()


def test_file_current_folder():
    src = get_fn('tz2.pdb')
    fh = FileManager(src)
    assert fh.use_filename
    assert not fh.is_compressed

    assert fh.ext.endswith('pdb')
    assert fh.is_filename
    assert fh.read() == os.path.relpath(fh.src)

    with open(src) as src2:
        fh2 = FileManager(src2)
        assert not fh2.is_filename
        fh2_content = fh2.read()
        src2.seek(0)
        assert fh2_content == src2.read()

    with open(src) as src2:
        fh_content = fh.read()
        src2.seek(0)
        assert fh_content != src2.read()

    fh3 = FileManager(src)
    content = open(src, 'r').read()
    assert_equal(fh3.read(force_buffer=True), content)

    # blob
    fh4 = FileManager(content)
    assert_equal(fh4.read(force_buffer=True), content)


def test_file_gz():
    src = get_fn('tz2_2.pdb.gz')
    fh = FileManager(src)
    assert fh.use_filename
    assert fh.is_compressed
    assert_equal(fh.compressed_ext, 'gz')

    assert fh.ext.endswith('pdb')
    assert fh.is_filename

    with open(src, 'rb') as src2:
        fh2 = FileManager(src2, compressed=True)
        assert not fh2.is_filename
        assert fh.is_compressed

        def func():
            fh2.ext

        with pytest.raises(ValueError):
            func()

    with open(src, 'rb') as src3:
        fh3 = FileManager(src3, compressed=True, ext='pdb')
        assert fh3.ext.endswith('pdb')

    # specify compression
    fh4 = FileManager(src, compressed=True)
    assert fh4.is_compressed

    content = gzip.open(src).read()
    assert_equal(fh4.read(force_buffer=True), content)


def test_file_passing_blob():
    src = get_fn('tz2.pdb')
    blob = open(src).read()

    fm = FileManager(blob)
    assert not fm.is_filename
    with pytest.raises(ValueError):
        fm.ext


def test_file_passing_blob_from_gzip():
    import gzip
    src = get_fn('tz2_2.pdb.gz')
    blob = gzip.open(src).read()

    fm = FileManager(blob)
    assert not fm.is_filename

    with pytest.raises(ValueError):
        fm.ext


def test_get_repr_names_from_dict():
    assert_equal(py_utils.get_repr_names_from_dict(repr_dict_example, 0),
                 ['cartoon', 'base', 'ball+stick'])


def test_js_utils():
    js_utils.launch_qtconsole()
    js_utils.clean_empty_output_area()
    js_utils.clean_error_output()
    js_utils._set_ipython_cell()
    js_utils.ngl_demo()
    js_utils.init_funcs()
    js_utils.hide_toolbar()
    js_utils.show_toolbar()
    js_utils.execute('print("hello")')

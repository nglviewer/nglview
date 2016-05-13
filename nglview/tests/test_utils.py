from __future__ import print_function
import unittest
from nglview.utils import seq_to_string, _camelize, _camelize_dict, FileManager
import nose.tools as nt

def test_seq_to_string():
    nt.assert_equal(seq_to_string([1, 2, 3]), '@1,2,3')
    nt.assert_equal(seq_to_string('@1,2,3'), '@1,2,3')

def test_camelize():
    nt.assert_equal(_camelize('remote_call'), 'remoteCall')
    nt.assert_equal(_camelize('flat_shaded'), 'flatShaded')
    nt.assert_equal(_camelize('near_clip'), 'nearClip')
    nt.assert_equal(_camelize('radius_type'), 'radiusType')
    nt.assert_equal(_camelize('aspect_ratio'), 'aspectRatio')
    nt.assert_equal(_camelize('radius_segments'), 'radiusSegments')
    nt.assert_equal(_camelize('sphere_detail'), 'sphereDetail')
    nt.assert_equal(_camelize('contact_type'), 'contactType')
    nt.assert_equal(_camelize('max_distance'), 'maxDistance')
    nt.assert_equal(_camelize('max_angle'), 'maxAngle')
    nt.assert_equal(_camelize('label_type'), 'labelType')
    nt.assert_equal(_camelize('radial_segments'), 'radialSegments')
    nt.assert_equal(_camelize('surface_type'), 'surfaceType')
    nt.assert_equal(_camelize('probe_radius'), 'probeRadius')
    nt.assert_equal(_camelize('scale_factor'), 'scaleFactor')

    # make sure this method does not change the camel
    nt.assert_equal(_camelize('maxAngle'), 'maxAngle')
    nt.assert_equal(_camelize('maxDistance'), 'maxDistance')
    nt.assert_equal(_camelize('color'), 'color')

def test_dict():
    kwargs = dict(default_representation=True)
    kwargs2 = _camelize_dict(kwargs)
    nt.assert_true('defaultRepresentation' in kwargs2)

def test_file_not_current_or_subfolder():
    src = '../__init__.py'
    fh = FileManager(src)
    nt.assert_false(fh.current_or_subfolder)
    nt.assert_false(fh.compressed)

    nt.assert_true(fh.ext, 'py')
    nt.assert_true(fh.is_filename)

    with open(src) as src2:
        fh2 = FileManager(src2)
        nt.assert_false(fh2.is_filename)
        nt.assert_true(fh2.read(), src2.read())


def test_file_current_folder():
    src = 'data/tz2.pdb'
    fh = FileManager(src)
    nt.assert_true(fh.current_or_subfolder)
    nt.assert_false(fh.compressed)

    nt.assert_true(fh.ext, 'pdb')
    nt.assert_true(fh.is_filename)
    nt.assert_true(fh.read(), fh.src)

    with open(src) as src2:
        fh2 = FileManager(src2)
        nt.assert_false(fh2.is_filename)
        nt.assert_true(fh2.read(), src2.read())

    with open(src) as src2:
        nt.assert_true(fh.read(), src2.read())

def test_file_gz():
    src = 'data/tz2_2.pdb.gz'
    fh = FileManager(src)
    nt.assert_true(fh.current_or_subfolder)
    nt.assert_true(fh.compressed)

    nt.assert_true(fh.ext, 'pdb')
    nt.assert_true(fh.is_filename)

    with open(src, 'rb') as src2:
        fh2 = FileManager(src2, compressed=True)
        nt.assert_false(fh2.is_filename)
        nt.assert_true(fh.compressed)

        def func():
            fh2.ext
        nt.assert_raises(ValueError, func)

    with open(src, 'rb') as src3:
        fh3 = FileManager(src3, compressed=True, ext='pdb')
        nt.assert_true(fh3.ext, 'pdb')

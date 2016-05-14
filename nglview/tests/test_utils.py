from __future__ import print_function
import unittest
from nglview.utils import seq_to_string, _camelize, _camelize_dict, FileManager
import nose.tools as nt
import gzip

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

def test_file_not_use_filename():
    src = '../__init__.py'
    fh = FileManager(src)
    nt.assert_false(fh.use_filename)
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
    nt.assert_true(fh.use_filename)
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

    fh3 = FileManager(src)
    content = open(src, 'rb').read() 
    nt.assert_equal(fh3.read(force_buffer=True), content)

    # blob
    fh4 = FileManager(content)
    nt.assert_equal(fh4.read(force_buffer=True), content)

def test_file_gz():
    src = 'data/tz2_2.pdb.gz'
    fh = FileManager(src)
    nt.assert_true(fh.use_filename)
    nt.assert_true(fh.compressed)
    nt.assert_equal(fh.compressed_ext, 'gz')

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

    # specify compression
    fh4 = FileManager(src, compressed=True)
    nt.assert_true(fh4.compressed)

    content = gzip.open(src).read() 
    nt.assert_equal(fh4.read(force_buffer=True), content)

def test_file_passing_blob():
    src = 'data/tz2.pdb'
    blob = open(src).read()

    fm = FileManager(blob)
    nt.assert_false(fm.is_filename)
    nt.assert_raises(ValueError, lambda: fm.ext)

def test_file_passing_blob_from_gzip():
    import gzip
    src = 'data/tz2_2.pdb.gz'
    blob = gzip.open(src).read()

    fm = FileManager(blob)
    nt.assert_false(fm.is_filename)

    nt.assert_raises(ValueError, lambda: fm.ext)

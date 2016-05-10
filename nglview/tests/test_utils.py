#!/usr/bin/env python
from __future__ import print_function
import unittest
from nglview.utils import seq_to_string, _camelize
import nose.tools as nt

def test_seq_to_string():
    nt.assert_equal(seq_to_string([1, 2, 3]), '@1,2,3')

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

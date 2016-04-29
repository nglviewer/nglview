#!/usr/bin/env python
from __future__ import print_function
import unittest

def test_utils():
    from nglview.utils import seq_to_string
    assert seq_to_string([1, 2, 3]) == '@1,2,3', 'convert string'

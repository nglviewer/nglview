#!/usr/bin/env python
from __future__ import print_function
import unittest

class TestUtils(unittest.TestCase):
    '''testing utils
    '''

    def test_utils(self):
        from nglview.utils import seq_to_string
        assert seq_to_string([1, 2, 3]) == '@1,2,3', 'convert string'


if __name__ == "__main__":
    unittest.main()

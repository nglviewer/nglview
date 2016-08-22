from __future__ import print_function, absolute_import
import uuid

"""abstract base class.
"""

__all__ = ['Structure', 'Trajectory']

class Structure(object):
    """abstract base class
    """
    def __init__(self):
        self.ext = "pdb"
        self.params = {}
        self.id = str(uuid.uuid4())

    def get_structure_string(self):
        raise NotImplementedError()


class Trajectory(object):
    """abstract base class
    """
    def __init__(self):
        self.id = str(uuid.uuid4())
        self.shown = True

    def get_coordinates(self, index):
        raise NotImplementedError()

    @property
    def n_frames(self):
        raise NotImplementedError()

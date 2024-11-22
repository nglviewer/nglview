from ..adaptor import Trajectory, Structure
import numpy as np

def get_mocked_traj():
    class MockedTraj(Structure, Trajectory):
        def __init__(self):
            Structure.__init__(self)
            Trajectory.__init__(self)

        @property
        def n_frames(self):
            return 5

        def get_coordinates(self, frame):
            coordinates = [
                [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                [[0.1, 0.1, 0.1], [1.1, 0.1, 0.1], [0.1, 1.1, 0.1]],
                [[0.2, 0.2, 0.2], [1.2, 0.2, 0.2], [0.2, 1.2, 0.2]],
                [[0.3, 0.3, 0.3], [1.3, 0.3, 0.3], [0.3, 1.3, 0.3]],
                [[0.4, 0.4, 0.4], [1.4, 0.4, 0.4], [0.4, 1.4, 0.4]]
            ]
            return np.array(coordinates)[frame]

        def get_structure_string(self):
            return 'hello'

    return MockedTraj()
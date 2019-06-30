from .utils.py_utils import _camelize_dict


class Stage:
    """
    Try to mimic NGL.Stage
    """

    def __init__(self, view):
        self._view = view

    def set_parameters(self, **kwargs):
        self._view._remote_call('setParameters',
                                target='Stage',
                                kwargs=_camelize_dict(kwargs))

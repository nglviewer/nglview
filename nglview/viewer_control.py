class ViewerControl:
    """EXPERIMENTAL. Create viewer controls (rotate, zoom, translation, ...)

    Parameters
    ----------
    view : nglview.NGLWidget

    Notes
    -----
    Unstable feature

    Examples
    --------
    >>> import nglview as nv
    >>> from nglview.viewer_control import ViewerControl
    >>> view = nv.NGLWidget()
    >>> view # doctest: +SKIP
    >>> control = ViewerControl(view=view)
    >>> control.zoom(0.1)
    """

    def __init__(self, view):
        self.view = view

    def _call(self, funcname, *args):
        self.view._remote_call(funcname, target='viewerControls', args=args)

    def _view_xz_plane(self):
        self.view._remote_call('viewXZPlane', target='Widget')

    def align(self, basis):
        '''
        
        Parameters
        ----------
        basis : List[float], len=16
        '''
        self._call('align', basis)

    def apply_matrix(self, basis):
        '''
        
        Parameters
        ----------
        basis : List[float], len=16
        '''
        self._call('applyMatrix', basis)

    def center(self, vector):
        '''
        
        Parameters
        ----------
        vector: List[float], len=3
        '''
        self._call('center', vector)

    def orient(self, basis):
        '''
        
        Parameters
        ----------
        basis : List[float], len=16
        '''
        self._call('orient', basis)

    def rotate(self, basis):
        '''
        
        Parameters
        ----------
        basis : List[float], len=4
            quaternion
        '''
        self._call('rotate', basis)

    def translate(self, vector):
        '''
        
        Parameters
        ----------
        basis : List[float], len=3
        '''
        self._call('translate', vector)

    def spin(self, axis, angle):
        '''
        
        Parameters
        ----------
        axis: List[float], len=3
        angle : float
        '''
        self._call('spin', axis, angle)

    def zoom(self, delta):
        '''
        
        Parameters
        ----------
        delta : float
        '''
        self._call('zoom', delta)

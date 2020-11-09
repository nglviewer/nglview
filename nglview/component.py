from .utils import widget_utils


class ComponentViewer:
    """Convenient attribute for NGLWidget. See example below.

    Examples
    --------
    >>> view = nv.NGLWidget() # doctest: +SKIP
    ... view.add_trajectory(traj) # traj is a component 0
    ... c = view.add_component(filename) # component 1
    ... c.clear()
    ... c.add_cartoon()
    ... c.add_licorice()
    ... view.remove_component(c)
    """

    def __init__(self, view, id):
        self._view = view
        self._id = id
        widget_utils._add_repr_method_shortcut(self, self._view)
        self._borrow_attribute(self._view, [
            'clear_representations', '_remove_representations_by_name',
            '_update_representations_by_name', 'center_view', 'center',
            'clear', 'set_representations'
        ], ['get_structure_string', 'get_coordinates', 'n_frames'])

    @property
    def id(self):
        return self._id

    @property
    def _index(self):
        # FIXME: not use private attribute from `self._view`
        return self._view._ngl_component_ids.index(self._id)

    def set_coordinates(self, coordinates):
        """

        Parameters
        ----------
        coordinates : numpy.ndarray, shape=(3, n_atoms)
        """
        self._view.set_coordinates({self._index: coordinates})

    def hide(self):
        """set invisibility for given components (by their indices)
        """
        self._view._remote_call("setVisibility",
                                target='compList',
                                args=[
                                    False,
                                ],
                                kwargs={'component_index': self._index})
        traj = self._view._get_traj_by_id(self.id)
        if traj is not None:
            traj.shown = False

    def show(self):
        """set invisibility for given components (by their indices)
        """
        self._view._remote_call("setVisibility",
                                target='compList',
                                args=[
                                    True,
                                ],
                                kwargs={'component_index': self._index})

        traj = self._view._get_traj_by_id(self.id)
        if traj is not None:
            traj.shown = True

    def set_position(self, pos):
        """
        Parameters
        ----------
        pos: List-like of float, length = 3
        """
        self._call('setPosition', pos)

    def set_rotation(self, rot):
        """
        Parameters
        ----------
        rot: List-like of float, length = 3
        """
        self._call('setRotation', rot)

    def set_scale(self, scale):
        """
        Parameters
        ----------
        scale: float
        """
        self._call('setScale', scale)

    def add_representation(self, repr_type, selection='all', **kwargs):
        kwargs['component'] = self._index
        self._view.add_representation(repr_type=repr_type,
                                      selection=selection,
                                      **kwargs)

    def _borrow_attribute(self, view, attributes, trajectory_atts=None):
        from functools import partial
        from types import MethodType

        traj = view._get_traj_by_id(self.id)

        for attname in attributes:
            view_att = getattr(view, attname)
            setattr(self, '_' + attname, MethodType(view_att, view))
            self_att = partial(getattr(view, attname), component=self._index)
            setattr(self, attname, self_att)

        if traj is not None and trajectory_atts is not None:
            for attname in trajectory_atts:
                traj_att = getattr(traj, attname)
                setattr(self, attname, traj_att)

    def _call(self, method, *args, **kwargs):
        """

        >>> c = view.add_component('file.pdb') # doctest: +SKIP
        ... c._call('setPosition', [10, 20, 0])
        ... c._call("setRotation", [1, 2, 0])
        """
        kwargs2 = {}
        kwargs2.update(kwargs)
        kwargs2['component_index'] = self._index
        self._view._remote_call(method,
                                target='compList',
                                args=args,
                                kwargs=kwargs2)

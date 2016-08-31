
from . import default

class Shape(object):
    """TODO: doc

    Parameters
    ----------
    view : nglview.NGLWidget

    Notes
    -----
    Unstable feature

    Examples
    --------
    >>> import nglview as nv
    >>> view = nv.NGLWidget()
    >>> view
    >>> shape = nv.Shape(view=view)
    >>> # TODO: add example
    >>> shape.add_sphere(...)

    See also
    --------
    {ngl_url}
    """.format(ngl_url=default.NGL_BASE_URL)

    def __init__(self, view):
        self.view = view
        names = ['mesh', 'sphere', 'ellipsoid', 'cylinder', 'cone', 'arrow']
        self._make_func(names)

    def _make_func(self, names):
        from types import MethodType

        def make_func(name):
            def func(this, *args):
                args_with_name = [name, ] + list(args)
                self.add(*args_with_name)
            func.__doc__ = 'check `add` method'
            return func

        for name in names:
            func_name = 'add_' + name
            func = make_func(name)
            setattr(self, func_name, MethodType(func, self))

    def add(self, *args):
        """

        Examples
        --------
        # TODO : add me

        See also
        --------
        {ngl_url}
        """.format(ngl_url=default.NGL_BASE_URL)

        self.view._add_shape([args,])

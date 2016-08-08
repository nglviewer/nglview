
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
    >>> shape.add_sphere(from=[0, 0, 0], to=[10, 10, 10], color=[1, 0, 0], radius=1.5])

    See also
    --------
    http://arose.github.io/ngl/api/dev/Shape.html
    """

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
        >>> shape = nv.Shape(view)

        See also
        --------
        http://arose.github.io/ngl/api/dev/Shape.html
        """
        self.view._add_shape([args,])

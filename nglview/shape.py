from .utils.py_utils import _update_url

SHAPE_EXAMPLES = {
    'mesh':
    """
# add_mesh(position, color, index, normal, name)
>>> shape.add_mesh(
        [ 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1 ],
        [ 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0 ]
)
""",
    'sphere':
    """
# add_sphere(position, color, radius, name)
>>> shape.add_sphere([0, 0, 9], [1, 0, 0], 1.5)
""",
    'ellipsoid':
    """
# add_ellipsoid(position, color, radius, majorAxis, minorAxis, name)
>>> shape.add_ellipsoid([ 6, 0, 0], [ 1, 0, 0 ], 1.5, [ 3, 0, 0 ], [ 0, 2, 0 ])
""",
    'cylinder':
    """
# add_cylinder(position1, position2, color, radius, name)
>>> shape.add_cylinder( [ 0, 2, 7 ], [ 0, 0, 9 ], [ 1, 1, 0 ], 0.5 )
""",
    'cone':
    """
# add_cone(position1, position2, color, radius, name) 
>>> shape.add_cone( [ 0, 2, 7 ], [ 0, 3, 3 ], [ 1, 1, 0 ], 1.5 )
""",
    'arrow':
    """
# add_arrow(position1, position2, color, radius, name)
>>> shape.add_arrow( [ 0, 2, 7 ], [ 0, 0, 9 ], [ 1, 1, 0 ], 0.5 )
""",
    'label':
    """
# add_text(position, color, size, text)
>>> shape.add_text( [ 10, -2, 4 ], [ 0.2, 0.5, 0.8 ], 0.5, "Hello" )
""",
    'text':
    """
# add_text(position, color, size, text)
>>> shape.add_text( [ 10, -2, 4 ], [ 0.2, 0.5, 0.8 ], 0.5, "Hello" )
"""
}


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
    """

    def __init__(self, view):
        self.view = view
        self._make_func(SHAPE_EXAMPLES.keys())

    def _make_func(self, names):
        from types import MethodType

        def make_func(name):
            def func(this, *args):
                args_with_name = [
                    name,
                ] + list(args)
                self.add(*args_with_name)

            func.__doc__ = SHAPE_EXAMPLES[name]
            return func

        for name in names:
            func_name = 'add_' + name
            func = make_func(name)
            setattr(self, func_name, MethodType(func, self))

    @_update_url
    def add(self, *args):
        """

        Examples
        --------
        >>> view.shape.add('text', [0, 4, -1], [0.2, 0.5, 0.8], 2.5, 'Meow')

        See also
        --------
        {ngl_url}
        """

        self.view._add_shape([
            args,
        ])

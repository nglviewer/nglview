import time

from ..parameters import REPRESENTATION_NAME_PAIRS


def wait(widget, attribute='value', timeout=5):
    """ EXPERIMENTAL. Require `ipython_blocking` package.

    Block further code execution until `attribute` of `widget` is updated.
    """
    from ipython_blocking import CaptureExecution
    c = CaptureExecution()
    with c:
        t0 = time.time()
        while True or time.time() - t0 > timeout:
            attr = getattr(widget, attribute)
            if attr:
                break
            c.step()
    return widget


def get_widget_by_name(box, widget_name):

    if hasattr(box, '_ngl_children'):
        children = box._ngl_children
    elif hasattr(box, 'children'):
        children = box.children
    else:
        children = None

    if children is not None:
        for widget in children:
            if hasattr(widget,
                       '_ngl_name') and widget._ngl_name == widget_name:
                return widget
    return None


def _add_repr_method_shortcut(self, other):
    from types import MethodType

    def make_func_add(rep):
        """return a new function object
        """

        def func(this, selection='all', **kwargs):
            """
            """
            self.add_representation(repr_type=rep[1],
                                    selection=selection,
                                    **kwargs)

        func.__doc__ = """Shortcut for `add_representation` method

        Examples
        --------
        >>> view.add_{name}()
        >>> # is equal to
        >>> view.add_representation('{name}')
        """.format(name=rep[0])
        return func

    def make_func_remove(rep):
        """return a new function object
        """

        def func(this, **kwargs):
            """
            """
            self._remove_representations_by_name(repr_name=rep[1], **kwargs)

        return func

    def make_func_update(rep):
        """return a new function object
        """

        def func(this, **kwargs):
            """
            """
            self._update_representations_by_name(repr_name=rep[1], **kwargs)

        return func

    for rep in REPRESENTATION_NAME_PAIRS:
        for make_func, root_fn in [(make_func_add, 'add'),
                                   (make_func_update, 'update'),
                                   (make_func_remove, 'remove')]:
            func = make_func(rep)
            fn = '_'.join((root_fn, rep[0]))
            setattr(self, fn, MethodType(func, other))

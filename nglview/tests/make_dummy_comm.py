from ipykernel.comm import Comm
from ipywidgets import Widget

#-----------------------------------------------------------------------------
# Utility stuff from ipywidgets tests
#
# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

#-----------------------------------------------------------------------------


class DummyComm(Comm):
    comm_id = 'a-b-c-d'

    def open(self, *args, **kwargs):
        pass

    def send(self, *args, **kwargs):
        pass

    def close(self, *args, **kwargs):
        pass


_widget_attrs = {}
displayed = []
undefined = object()


def setup():
    _widget_attrs['_comm_default'] = getattr(Widget, '_comm_default',
                                             undefined)
    Widget._comm_default = lambda self: DummyComm()
    _widget_attrs['_ipython_display_'] = Widget._ipython_display_

    def raise_not_implemented(*args, **kwargs):
        raise NotImplementedError()

    Widget._ipython_display_ = lambda _: _


def teardown():
    for attr, value in _widget_attrs.items():
        if value is undefined:
            delattr(Widget, attr)
        else:
            setattr(Widget, attr, value)

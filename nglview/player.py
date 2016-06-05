from ipywidgets import DOMWidget
from traitlets import Int, Bool, Dict, Float
from traitlets import observe

class TrajectoryPlayer(DOMWidget):
    step = Int().tag(sync=True)
    sync_frame = Bool().tag(sync=True)
    delay = Float().tag(sync=True)
    parameters = Dict().tag(sync=True)

    def __init__(self, view, step=1, delay=0.1, sync_frame=False):
        self._view = view
        self.step = step
        self.sync_frame = sync_frame
        self.delay = delay

    @property
    def frame(self):
        return self._view.frame

    @frame.setter
    def frame(self, value):
        self._view.frame = value

    @property
    def count(self):
        return self._view.count

    @observe('sync_frame')
    def update_sync_frame(self, change):
        value = change['new']
        if value:
            self._view._set_sync_frame()
        else:
            self._view._set_unsync_frame()

    @observe("delay")
    def update_delay(self, change):
        delay = change['new']
        self._view._set_delay(delay)

    @observe('parameters')
    def update_parameters(self, change):
        params = change['new']
        self.sync_frame = params.get("sync_frame", self.sync_frame)
        self.delay = params.get("delay", self.delay)
        self.step = params.get("step", self.step)

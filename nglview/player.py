from ipywidgets import DOMWidget
from traitlets import Int, Bool, Dict, Float

class TrajectoryPlayer(DOMWidget):
    step = Int() 
    _sync_frame = Bool()
    _delay = Float()
    _params = Dict()

    def __init__(self, view, step=1, delay=0.1, sync_frame=False):
        self._view = view
        self.step = step
        self._sync_frame = sync_frame
        self._delay = delay

        self._params = dict(sync_frame=self.sync_frame,
                            delay=self.delay,
                            step=self.step)

    @property
    def frame(self):
        return self._view.frame

    @frame.setter
    def frame(self, value):
        return self._view.frame = value

    @property
    def count(self):
        return self._view.count

    @property
    def sync_frame(self):
        return self._sync_frame

    @sync_frame.setter
    def sync_frame(self, value):
        if value:
            self._view._set_sync_frame()
            self._sync_frame = True
        else:
            self._view._set_unsync_frame()
            self._sync_frame = False

    @property
    def delay(self):
        return self._delay

    @delay.setter
    def delay(self, delay):
        self._delay = delay
        self._view._set_delay(delay)

    @property
    def parameters(self):
        return self._params

    @parameters.setter
    def parameters(self, params):
        self.sync_frame = params.get("sync_frame", self.sync_frame)
        self.delay = params.get("delay", self.delay)
        self.step = params.get("step", self.step)

        self._params = dict(sync_frame=self.sync_frame,
                            delay=self.delay,
                            step=self.step)

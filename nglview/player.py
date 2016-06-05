from ipywidgets import DOMWidget, IntText, BoundedFloatText, HBox, VBox, Checkbox
from traitlets import Int, Bool, Dict, Float
from traitlets import observe, link

class TrajectoryPlayer(DOMWidget):
    # should set default values here different from desired defaults
    # so `observe` can be triggered
    step = Int(0).tag(sync=True)
    sync_frame = Bool(True).tag(sync=True)
    delay = Float(0.0).tag(sync=True)
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

    def _display(self):
        int_text = IntText(self.step, description='step')
        float_txt = BoundedFloatText(self.delay, description='delay', min=0.001)
        checkbox_sync_frame = Checkbox(self.sync_frame, description='sync_frame')

        link((int_text, 'value'), (self, 'step'))
        link((float_txt, 'value'), (self, 'delay'))
        link((checkbox_sync_frame, 'value'), (self, 'sync_frame'))

        return VBox([int_text, float_txt, checkbox_sync_frame])

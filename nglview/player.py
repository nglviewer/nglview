from ipywidgets import (DOMWidget, IntText, FloatText, HBox, VBox, Checkbox,
                        ColorPicker, FloatSlider,
                        Dropdown)
                        
from traitlets import Int, Bool, Dict, Float, CaselessStrEnum
from traitlets import observe, link

class TrajectoryPlayer(DOMWidget):
    # should set default values here different from desired defaults
    # so `observe` can be triggered
    step = Int(0).tag(sync=True)
    sync_frame = Bool(True).tag(sync=True)
    interpolate = Bool(False).tag(sync=False)
    delay = Float(0.0).tag(sync=True)
    parameters = Dict().tag(sync=True)
    iparams = Dict().tag(sync=False)
    _t_interpolation = Float().tag(sync=False)
    _type_interpolation = CaselessStrEnum(['linear', 'spline']).tag(sync=False)

    def __init__(self, view, step=1, delay=100, sync_frame=False, min_delay=40):
        self._view = view
        self.step = step
        self.sync_frame = sync_frame
        self.delay = delay
        self.min_delay = min_delay
        self._t_interpolation = 0.5
        self._type_interpolation = 'linear'
        self.iparams = dict(t=self._t_interpolation, step=1, type=self._type_interpolation)

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

    @observe('_t_interpolation')
    def _t_interpolation_changed(self, change):
        self.iparams['t'] = change['new']

    @observe('_type_interpolation')
    def _t_interpolation_changed(self, change):
        self.iparams['type'] = change['new']

    def _display(self):
        step_text = IntText(self.step, description='step')
        delay_text = FloatText(value=self.delay, description='delay')
        checkbox_interpolate = Checkbox(self.interpolate, description='interpolate')
        bg_color = ColorPicker(value='white', description='background_color')
        t_interpolation = FloatSlider(value=0.5, min=0, max=1.0, step=0.1)
        type_iterpolation = Dropdown(value=self._type_interpolation,
                options=['linear', 'spline'], description='interpolation type')

        link((step_text, 'value'), (self, 'step'))
        link((delay_text, 'value'), (self, 'delay'))
        link((checkbox_interpolate, 'value'), (self, 'interpolate'))
        link((t_interpolation, 'value'), (self, '_t_interpolation'))
        link((type_iterpolation, 'value'), (self, '_type_interpolation'))
        link((bg_color, 'value'), (self._view, 'background'))

        return VBox([step_text, delay_text, bg_color,
                     checkbox_interpolate,
                     t_interpolation,
                     type_iterpolation])

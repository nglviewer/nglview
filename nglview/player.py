# TODO: reorg
import ipywidgets
from ipywidgets import (DOMWidget, IntText, FloatText, HBox, VBox, Checkbox,
                        ColorPicker, IntSlider, FloatSlider,
                        Dropdown,
                        interactive)

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
    _interpolation_t = Float().tag(sync=False)
    _iterpolation_type = CaselessStrEnum(['linear', 'spline']).tag(sync=False)
    spin = Bool(False).tag(sync=False)
    _spin_x = Int(1).tag(sync=False)
    _spin_y = Int(0).tag(sync=False)
    _spin_z = Int(0).tag(sync=False)
    _spin_speed = Float(0.005).tag(sync=False)
    camera = CaselessStrEnum(['perspective', 'orthographic'],
        default_value='perspective').tag(sync=False)

    def __init__(self, view, step=1, delay=100,
                 sync_frame=False, min_delay=40):
        self._view = view
        self.step = step
        self.sync_frame = sync_frame
        self.delay = delay
        self.min_delay = min_delay
        self._interpolation_t = 0.5
        self._iterpolation_type = 'linear'
        self.iparams = dict(
            t=self._interpolation_t,
            step=1,
            type=self._iterpolation_type)

    @observe('camera')
    def on_camera_changed(self, change):
        camera_type = change['new']
        self._view._remote_call("setParameters",
                                target='Stage',
                                kwargs=dict(cameraType=camera_type))

    @observe('camera')
    def on_camera_changed(self, change):
        camera_type = change['new']
        self._view._remote_call("setParameters",
                target='Stage',
                kwargs=dict(cameraType=camera_type))

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

    @observe('_interpolation_t')
    def _interpolation_t_changed(self, change):
        self.iparams['t'] = change['new']

    @observe('_iterpolation_type')
    def _interpolation_t_changed(self, change):
        self.iparams['type'] = change['new']

    @observe('spin')
    def on_spin_changed(self, change):
        self.spin = change['new']
        if self.spin:
            self._view._set_spin([self._spin_x, self._spin_y, self._spin_z],
                                 self._spin_speed)
        else:
            # stop
            self._view._set_spin(None, None)

    @observe('_spin_x')
    def on_spin_x_changed(self, change):
        self._spin_x = change['new']
        if self.spin:
            self._view._set_spin([self._spin_x, self._spin_y, self._spin_z],
                                 self._spin_speed)

    @observe('_spin_y')
    def on_spin_y_changed(self, change):
        self._spin_y = change['new']
        if self.spin:
            self._view._set_spin([self._spin_x, self._spin_y, self._spin_z],
                                 self._spin_speed)

    @observe('_spin_z')
    def on_spin_z_changed(self, change):
        self._spin_z = change['new']
        if self.spin:
            self._view._set_spin([self._spin_x, self._spin_y, self._spin_z],
                                 self._spin_speed)

    @observe('_spin_speed')
    def on_spin_speed_changed(self, change):
        self._spin_speed = change['new']
        if self.spin:
            self._view._set_spin([self._spin_x, self._spin_y, self._spin_z],
                                 self._spin_speed)

    def _display(self):
        step_slide = IntSlider(
            value=self.step,
            min=-100,
            max=100,
            description='step')
        delay_text = IntSlider(
            value=self.delay,
            min=10,
            max=1000,
            description='delay')
        checkbox_interpolate = Checkbox(
            self.interpolate, description='interpolate')
        checkbox_spin = Checkbox(self.spin, description='spin')
        spin_x_slide = IntSlider(
            self._spin_x,
            min=-1,
            max=1,
            description='spin_x')
        spin_y_slide = IntSlider(
            self._spin_y,
            min=-1,
            max=1,
            description='spin_y')
        spin_z_slide = IntSlider(
            self._spin_z,
            min=-1,
            max=1,
            description='spin_z')
        spin_speed_slide = FloatSlider(
            self._spin_speed,
            min=0,
            max=0.2,
            step=0.001,
            description='spin speed')
        bg_color = ColorPicker(value='white', description='background_color')
        # t_interpolation = FloatSlider(value=0.5, min=0, max=1.0, step=0.1)
        interpolation_type = Dropdown(value=self._iterpolation_type,
                                      options=['linear', 'spline'], description='interpolation type')

        camera_type = Dropdown(value=self.camera,
                               options=['perspective', 'orthographic'], description='camera')

        camera_type = Dropdown(value=self.camera,
                options=['perspective', 'orthographic'], description='camera')

        link((step_slide, 'value'), (self, 'step'))
        link((delay_text, 'value'), (self, 'delay'))
        link((checkbox_interpolate, 'value'), (self, 'interpolate'))
        # link((t_interpolation, 'value'), (self, '_interpolation_t'))
        link((interpolation_type, 'value'), (self, '_iterpolation_type'))
        link((camera_type, 'value'), (self, 'camera'))
        link((bg_color, 'value'), (self._view, 'background'))

        # spin
        link((checkbox_spin, 'value'), (self, 'spin'))
        link((spin_x_slide, 'value'), (self, '_spin_x'))
        link((spin_y_slide, 'value'), (self, '_spin_y'))
        link((spin_z_slide, 'value'), (self, '_spin_z'))
        link((spin_speed_slide, 'value'), (self, '_spin_speed'))

        ibox = HBox([checkbox_interpolate, interpolation_type])
        v0 = VBox([step_slide,
                   delay_text,
                   bg_color,
                   ibox,
                   camera_type])

        v1 = VBox([checkbox_spin,
                   spin_x_slide,
                   spin_y_slide,
                   spin_z_slide,
                   spin_speed_slide])

        genbox = HBox([v0, v1])
        prefbox = self._show_preference()
        tab = ipywidgets.Tab([genbox, prefbox])
        tab.set_title(0, 'General')
        tab.set_title(1, 'Speed')
        return tab

    def _show_preference(self):
        def func(pan_speed=0.8,
                 rotate_speed=2,
                 zoom_speed=1.2,
                 clip_dist=10):
            self._view.parameters = dict(
                panSpeed=pan_speed,
                rotateSpeed=rotate_speed,
                zoomSpeed=zoom_speed,
                clipDist=clip_dist)

        return interactive(func,
                  pan_speed=(0, 10, 0.1),
                  rotate_speed=(0, 10, 1),
                  zoom_speed=(0, 10, 1),
                  clip_dist=(0, 200, 5))

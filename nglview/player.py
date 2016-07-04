# TODO: reorg
# simplify code
import json
import ipywidgets
from IPython.display import display, Javascript
from ipywidgets import (DOMWidget, IntText, FloatText,
                        Box, HBox, VBox, Checkbox,
                        ColorPicker, IntSlider, FloatSlider,
                        Dropdown,
                        Button,
                        Text, Textarea,
                        interactive,
                        Tab)

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
    _render_params = Dict().tag(sync=False)

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
        self._render_params = dict(factor=4,
                                   antialias=True,
                                   trim=False,
                                   transparent=False)
        self.picked_widget = self._add_text_picked()
        self.repr_widget = self._add_text_repr_widget()

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
                                      options=['linear', 'spline'])

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

        qtconsole_button = self._add_button_qtconsole()

        ibox = HBox([checkbox_interpolate, interpolation_type])
        center_button = self._add_button_center()
        render_button = self._show_download_image()
        center_render_hbox = HBox([center_button, render_button])

        v0_left = VBox([step_slide,
                   delay_text,
                   bg_color,
                   ibox,
                   camera_type,
                   center_render_hbox,
                   qtconsole_button,
                   ])

        spinbox= VBox([checkbox_spin,
                   spin_x_slide,
                   spin_y_slide,
                   spin_z_slide,
                   spin_speed_slide])

        genbox = HBox([v0_left, ])
        prefbox = self._show_preference()
        themebox = Box([self._add_button_theme(), self._add_button_reset_theme()])
        hidebox = Box([])
        help_url = self._show_website()

        picked_box = HBox([self.picked_widget,])
        repr_box= HBox([self.repr_widget, self._add_repr_sliders()])

        extrabox = Tab([repr_box, picked_box])
        extrabox.set_title(0, 'repr')
        extrabox.set_title(1, 'picked atom')

        export_image_box = HBox([self._add_button_export_image()])

        tab = Tab([genbox, spinbox, prefbox, themebox, 
            hidebox, help_url, extrabox,
            export_image_box])

        tab.set_title(0, 'General')
        tab.set_title(1, 'Spin')
        tab.set_title(2, 'Speed')
        tab.set_title(3, 'Theme')
        tab.set_title(4, 'Hide')
        tab.set_title(5, 'Help')
        tab.set_title(6, 'Extra')
        tab.set_title(7, 'Image')
        return tab

    def _add_button_center(self):
        button = Button(description='Center')
        def on_click(button):
            self._view.center()
        button.on_click(on_click)
        return button

    def _add_button_theme(self):
        button = Button(description='Oceans16')
        def on_click(button):
            from nglview import theme
            display(theme.oceans16())
        button.on_click(on_click)
        return button

    def _add_button_reset_theme(self):
        from nglview.jsutils import js_clean_empty_output_area
        button = Button(description='Default')
        def on_click(button):
            display(Javascript('$("#nglview_style").remove()'))
            display(Javascript(js_clean_empty_output_area))
        button.on_click(on_click)
        return button

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

    def _show_download_image(self):
        # "interactive" does not work for True/False in ipywidgets 4 yet.
        button = Button(description='Screenshot')
        def on_click(button):
            self._view.download_image()
        button.on_click(on_click)
        return button

    def _make_button_url(self, url, description):
        from nglview.jsutils import js_open_url_template
        button = Button(description=description)

        def on_click(button):
            display(Javascript(js_open_url_template.format(url=url)))

        button.on_click(on_click)
        return button

    def _show_website(self):
        buttons = [self._make_button_url(url, description) for url, description in
            [("'http://arose.github.io/nglview/latest/'", "nglview"),
            ("'http://arose.github.io/ngl/api/dev/'", "NGL"),
            ("'http://arose.github.io/ngl/api/dev/tutorial-selection-language.html'", "Selection"),
            ("'http://arose.github.io/ngl/api/dev/tutorial-molecular-representations.html'", "Representation")]
        ]
        return HBox(buttons)

    def _add_button_qtconsole(self):
        from nglview.jsutils import js_launch_qtconsole
        button = Button(description='qtconsole')

        def on_click(button):
            display(Javascript(js_launch_qtconsole))
        button.on_click(on_click)
        return button

    def _add_text_picked(self):
        ta = Textarea(value=json.dumps(self._view.picked), description='Picked atom')
        return ta

    def _add_text_repr_widget(self):
        button_info = Button(description='Refresh', tooltip='Get representation info')
        button_update = Button(description='Update', tooltip='Update representation')
        bbox = HBox([button_info, button_update])
        repr_name = Text(value='', description='repr_name')
        component_slider = IntSlider(value=0, description='cindex')
        repr_slider = IntSlider(value=0, description='rindex')
        ta = Textarea(value='', description='rinfo')

        def on_click_info(button):
            self._view._request_repr_parameters(component=int(component_slider.value),
                                                repr_index=int(repr_slider.value))
        button_info.on_click(on_click_info)

        def on_click_update(button):
            parameters = json.loads(ta.value.replace("False", "false").replace("True", "true"))
            self._view.update_representation(component=int(component_slider.value),
                                             repr_index=int(repr_slider.value),
                                             **parameters)
            self._view._request_update_reprs()
        button_update.on_click(on_click_update)

        def update_slide_info(change):
            self._view._request_repr_parameters(component=int(component_slider.value),
                                                repr_index=int(repr_slider.value))
        repr_slider.observe(update_slide_info, names='value')
        component_slider.observe(update_slide_info, names='value')

        # NOTE: if you update below list, make sure to update _add_repr_sliders
        # or refactor
        return VBox([bbox, repr_name, component_slider, repr_slider, ta])

    def _add_repr_sliders(self):
        repr_checkbox = Checkbox(value=False, description='repr slider')

        vbox = VBox([repr_checkbox])

        def create_widget(change):
            if change['new']:
                # repr_name
                # TODO: correctly upate name
                name = self.repr_widget.children[1].value
                component_slider = self.repr_widget.children[2]
                repr_slider = self.repr_widget.children[3]
                widget = self._view._display_repr(component=int(component_slider.value),
                                         repr_index=int(repr_slider.value),
                                         name=name)
                vbox.children = [repr_checkbox, widget]
            else:
                vbox.children = [repr_checkbox, ]
        repr_checkbox.observe(create_widget, names='value')
        return vbox

    def _add_button_export_image(self):
        slider_factor = IntSlider(value=4, min=1, max=10, description='scale')
        checkbox_antialias = Checkbox(value=True, description='antialias')
        checkbox_trim = Checkbox(value=False, description='trim')
        checkbox_transparent = Checkbox(value=False, description='transparent')
        filename_text = Text(value='Screenshot', description='Filename')

        button = Button(description='Export Image')

        def on_click(button):
            self._view.download_image(factor=slider_factor.value,
                    antialias=checkbox_antialias.value,
                    trim=checkbox_trim.value,
                    transparent=checkbox_transparent.value,
                    filename=filename_text.value)

        button.on_click(on_click)

        return VBox([button,
            filename_text,
            slider_factor,
            checkbox_antialias,
            checkbox_trim,
            checkbox_transparent])

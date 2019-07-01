# TODO: reorg
# simplify code
import json
import time
import uuid
from collections import defaultdict

from IPython.display import Javascript, display
from ipywidgets import (Box, Button, Checkbox, ColorPicker, Dropdown,
                        FloatSlider, FloatText, HBox, IntSlider, IntText,
                        Label, Layout, Play, Tab, Text, Textarea, ToggleButton,
                        VBox, interactive, jslink)
from traitlets import (Any, Bool, CaselessStrEnum, Dict, Float, HasTraits, Int,
                       link, observe)
import traitlets

from . import default
from .layout import (_make_autofit, _make_box_layout, _make_delay_tab,
                     _relayout, _relayout_master, make_form_item_layout)
from .parameters import REPRESENTATION_NAMES
from .utils import js_utils


def _dry_run(v, func, *args, **kwargs):
    # Get remote call message.
    msg = []

    def _another_call(*args, **kwargs):
        msg.append(v._get_remote_call_msg(*args, **kwargs))

    old = v._remote_call
    setattr(v, '_remote_call', _another_call)
    func(*args, **kwargs)
    setattr(v, '_remote_call', old)

    return msg and msg[-1] or None  # FIXME: all messages?


def _cature_msg(f):
    def wrap(*args, **kwargs):
        self = args[0]
        v = self._view
        msg_dict = v._player_dict
        w = f(*args, **kwargs)
        # FIXME: add code here
        # if isinstance(w, Button):
        #     msg_dict[w.model_id] = _dry_run(v, w.click)
        return w
    return wrap


class TrajectoryPlayer(HasTraits):
    # should set default values here different from desired defaults
    # so `observe` can be triggered
    step = Int(0)
    interpolate = Bool(False)
    delay = Float(0.0)
    parameters = Dict()
    iparams = Dict()
    _interpolation_t = Float()
    _iterpolation_type = CaselessStrEnum(['linear', 'spline'])
    camera = CaselessStrEnum(['perspective', 'orthographic'],
                             default_value='perspective')
    _render_params = Dict()
    _real_time_update = Bool(False)

    widget_tab = Any(None)
    widget_repr = Any(None)
    widget_repr_parameters = Any(None)
    widget_quick_repr = Any(None)
    widget_general = Any(None)
    widget_picked = Any(None)
    widget_preference = Any(None)
    widget_extra = Any(None)
    widget_help = Any(None)
    widget_export_image = Any(None)
    widget_component_slider = Any(None)
    widget_repr_slider = Any(None)
    widget_repr_choices = Any(None)
    widget_repr_control_buttons = Any(None)
    widget_repr_add = Any(None)
    widget_accordion_repr_parameters = Any(None)
    widget_repr_parameters_dialog = Any(None)
    widget_repr_name = Any(None)
    widget_component_dropdown = Any(None)
    widget_player = Any(None)
    widget_background = Any(None)
    widget_camera = Any(None)
    btn_center = Any(None)

    def __init__(self, view, step=1, delay=100, min_delay=40):
        self._view = view
        self.step = step
        self.delay = delay
        self.min_delay = min_delay
        self._iplayer = None
        self._interpolation_t = 0.5
        self._iterpolation_type = 'linear'
        self.iparams = dict(t=self._interpolation_t,
                            step=1,
                            type=self._iterpolation_type)
        self._render_params = dict(factor=4,
                                   antialias=True,
                                   trim=False,
                                   transparent=False)

        self._widget_names = [w for w in dir(self) if w.startswith('wiget_')]

    def _update_padding(self, padding=default.DEFAULT_PADDING):
        widget_collection = [
            self.widget_general, self.widget_repr, self.widget_preference,
            self.widget_repr_parameters, self.widget_help, self.widget_extra,
            self.widget_picked
        ]
        for widget in widget_collection:
            if widget is not None:
                widget.layout.padding = padding

    def _create_all_widgets(self, thread=False):

        def make():
            self.widget_tab = self._display()

            old_index = self.widget_tab.selected_index
            for index, _ in enumerate(self.widget_tab.children):
                self.widget_tab.selected_index = index
                if thread:
                    time.sleep(0.1)

            self.widget_tab.selected_index = old_index
        if thread:
            self._view._run_on_another_thread(make)
        else:
            make()

    def smooth(self):
        self.interpolate = True

    @observe('delay')
    def _on_delay(self, change):
        if self.widget_player:
            self.widget_player.interval = change['new']

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
    def max_frame(self):
        return self._view.max_frame

    @observe('parameters')
    def update_parameters(self, change):
        params = change['new']
        self.delay = params.get("delay", self.delay)
        self.step = params.get("step", self.step)

    @observe('_interpolation_t')
    def _interpolation_t_changed(self, change):
        self.iparams['t'] = change['new']

    def _display(self):
        if self.widget_tab is None:
            box_factory = [(self._make_general_box, 'General'),
                           (self._make_widget_repr, 'Representation'),
                           (self._make_widget_preference, 'Preference'),
                           (self._make_extra_box, 'Extra')]

            tab = _make_delay_tab(box_factory, selected_index=0)
            # tab = _make_autofit(tab)
            tab.layout.align_self = 'center'
            tab.layout.align_items = 'stretch'

            self.widget_tab = self._view._igui = tab

        return self.widget_tab

    def _make_widget_tab(self):
        return self._display()

    @_cature_msg
    def _make_button_center(self):
        self.btn_center = Button(description=' Center', icon='fa-bullseye')

        @self.btn_center.on_click
        def on_click(_):
            self._view.center()

        return self.btn_center

    def _make_widget_preference(self, width='100%'):
        def make_func():
            parameters = self._view._ngl_full_stage_parameters

            def func(pan_speed=parameters.get('panSpeed', 0.8),
                     rotate_speed=parameters.get('rotateSpeed', 2),
                     zoom_speed=parameters.get('zoomSpeed', 1.2),
                     clip_dist=parameters.get('clipDist', 10),
                     camera_fov=parameters.get('cameraFov', 40),
                     clip_far=parameters.get('clipFar', 100),
                     clip_near=parameters.get('clipNear', 0),
                     fog_far=parameters.get('fogFar', 100),
                     fog_near=parameters.get('fogNear', 50),
                     impostor=parameters.get('impostor', True),
                     light_intensity=parameters.get('lightIntensity', 1),
                     quality=parameters.get('quality', 'medium'),
                     sample_level=parameters.get('sampleLevel', 1)):

                self._view.parameters = dict(panSpeed=pan_speed,
                                             rotateSpeed=rotate_speed,
                                             zoomSpeed=zoom_speed,
                                             clipDist=clip_dist,
                                             clipFar=clip_far,
                                             clipNear=clip_near,
                                             cameraFov=camera_fov,
                                             fogFar=fog_far,
                                             fogNear=fog_near,
                                             impostor=impostor,
                                             lightIntensity=light_intensity,
                                             quality=quality,
                                             sampleLevel=sample_level)

            return func

        def make_widget_box():
            widget_sliders = interactive(make_func(),
                                         pan_speed=(0, 10, 0.1),
                                         rotate_speed=(0, 10, 1),
                                         zoom_speed=(0, 10, 1),
                                         clip_dist=(0, 200, 5),
                                         clip_far=(0, 100, 1),
                                         clip_near=(0, 100, 1),
                                         camera_fov=(15, 120, 1),
                                         fog_far=(0, 100, 1),
                                         fog_near=(0, 100, 1),
                                         light_intensity=(0, 10, 0.02),
                                         quality=['low', 'medium', 'high'],
                                         sample_level=(-1, 5, 1))

            for child in widget_sliders.children:
                if isinstance(child, (IntSlider, FloatSlider)):
                    child.layout.width = default.DEFAULT_SLIDER_WIDTH
            return widget_sliders

        if self.widget_preference is None:
            widget_sliders = make_widget_box()
            reset_button = Button(description='Reset')
            widget_sliders.children = [
                reset_button,
            ] + list(widget_sliders.children)

            @reset_button.on_click
            def on_click(reset_button):
                self._view.parameters = self._view._original_stage_parameters
                self._view._ngl_full_stage_parameters = self._view._original_stage_parameters
                widget_sliders.children = [
                    reset_button,
                ] + list(make_widget_box().children)

            self.widget_preference = _relayout_master(widget_sliders,
                                                      width=width)
        return self.widget_preference

    def _show_download_image(self):
        # "interactive" does not work for True/False in ipywidgets 4 yet.
        button = Button(description=' Screenshot', icon='fa-camera')

        @button.on_click
        def on_click(button):
            self._view.download_image()

        return button

    def _make_text_picked(self):
        ta = Textarea(value=json.dumps(self._view.picked),
                      description='Picked atom')
        ta.layout.width = '300px'
        _dry_run(self._view,
                 self._view._on_picked,
                 change={
                     'old': '',
                     'new': ''
                 })
        return ta

    def _refresh(self, component_slider, repr_slider):
        """update representation and component information
        """
        self._view._request_repr_parameters(component=component_slider.value,
                                            repr_index=repr_slider.value)
        self._view._update_repr_dict()
        self._view._handle_repr_dict_changed(change=dict(
            new=self._view._ngl_repr_dict))

    @_cature_msg
    def _make_button_repr_control(self, component_slider, repr_slider,
                                  repr_selection):
        button_refresh = Button(description=' Refresh',
                                tooltip='Get representation info',
                                icon='fa-refresh')
        button_center_selection = Button(description=' Center',
                                         tooltip='center selected atoms',
                                         icon='fa-bullseye')
        button_center_selection._ngl_name = 'button_center_selection'
        button_hide = Button(description=' Hide',
                             icon='fa-eye-slash',
                             tooltip='Hide/Show current representation')
        button_remove = Button(description=' Remove',
                               icon='fa-trash',
                               tooltip='Remove current representation')

        @button_refresh.on_click
        def on_click_refresh(button):
            self._refresh(component_slider, repr_slider)

        @button_center_selection.on_click
        def on_click_center(center_selection):
            self._view.center(selection=repr_selection.value,
                              component=component_slider.value)

        @button_hide.on_click
        def on_click_hide(button_hide):
            component = component_slider.value
            repr_index = repr_slider.value

            if button_hide.description == 'Hide':
                hide = True
                button_hide.description = 'Show'
            else:
                hide = False
                button_hide.description = 'Hide'

            self._view._remote_call('setVisibilityForRepr',
                                    target='Widget',
                                    args=[component, repr_index, not hide])

        @button_remove.on_click
        def on_click_remove(button_remove):
            self._view._remove_representation(component=component_slider.value,
                                              repr_index=repr_slider.value)
            self._view._request_repr_parameters(
                component=component_slider.value, repr_index=repr_slider.value)

        bbox = _make_autofit(
            HBox([
                button_refresh, button_center_selection, button_hide,
                button_remove
            ]))
        for b in bbox.children:
            if hasattr(b, 'click'):
                _dry_run(self._view, b.click)
        return bbox

    def _make_widget_repr(self):
        self.widget_repr_name = Text(value='', description='representation')
        self.widget_repr_name._ngl_name = 'repr_name_text'
        repr_selection = Text(value=' ', description='selection')
        repr_selection._ngl_name = 'repr_selection'
        repr_selection.width = self.widget_repr_name.width = default.DEFAULT_TEXT_WIDTH

        max_n_components = max(self._view.n_components - 1, 0)
        self.widget_component_slider = IntSlider(value=0,
                                                 max=max_n_components,
                                                 min=0,
                                                 description='component')
        self.widget_component_slider._ngl_name = 'component_slider'

        cvalue = ' '
        self.widget_component_dropdown = Dropdown(value=cvalue,
                                                  options=[
                                                      cvalue,
                                                  ],
                                                  description='component')
        self.widget_component_dropdown._ngl_name = 'component_dropdown'

        self.widget_repr_slider = IntSlider(value=0,
                                            description='representation',
                                            width=default.DEFAULT_SLIDER_WIDTH)
        self.widget_repr_slider._ngl_name = 'repr_slider'
        self.widget_repr_slider.visible = True

        self.widget_component_slider.layout.width = default.DEFAULT_SLIDER_WIDTH
        self.widget_repr_slider.layout.width = default.DEFAULT_SLIDER_WIDTH
        self.widget_component_dropdown.layout.width = self.widget_component_dropdown.max_width = default.DEFAULT_TEXT_WIDTH

        # turn off for now
        self.widget_component_dropdown.layout.display = 'none'
        self.widget_component_dropdown.description = ''

        # self.widget_accordion_repr_parameters = Accordion()
        self.widget_accordion_repr_parameters = Tab()
        self.widget_repr_parameters = self._make_widget_repr_parameters(
            self.widget_component_slider, self.widget_repr_slider,
            self.widget_repr_name)
        self.widget_accordion_repr_parameters.children = [
            self.widget_repr_parameters, Box()
        ]
        self.widget_accordion_repr_parameters.set_title(0, 'Parameters')
        self.widget_accordion_repr_parameters.set_title(1, 'Hide')
        self.widget_accordion_repr_parameters.selected_index = 1

        checkbox_reprlist = Checkbox(value=False, description='reprlist')
        checkbox_reprlist._ngl_name = 'checkbox_reprlist'
        self.widget_repr_choices = self._make_repr_name_choices(
            self.widget_component_slider, self.widget_repr_slider)
        self.widget_repr_choices._ngl_name = 'reprlist_choices'

        self.widget_repr_add = self._make_add_widget_repr(
            self.widget_component_slider)

        def on_update_checkbox_reprlist(change):
            self.widget_repr_choices.visible = change['new']

        checkbox_reprlist.observe(on_update_checkbox_reprlist, names='value')

        def on_repr_name_text_value_changed(change):
            name = change['new'].strip()
            old = change['old'].strip()

            should_update = (self._real_time_update and old and name
                             and name in REPRESENTATION_NAMES
                             and name != change['old'].strip())

            if should_update:
                component = self.widget_component_slider.value
                repr_index = self.widget_repr_slider.value
                self._view._remote_call(
                    'setRepresentation',
                    target='Widget',
                    args=[change['new'], {}, component, repr_index])
                self._view._request_repr_parameters(component, repr_index)

        def on_component_or_repr_slider_value_changed(change):
            self._view._request_repr_parameters(
                component=self.widget_component_slider.value,
                repr_index=self.widget_repr_slider.value)
            self.widget_component_dropdown.options = tuple(
                self._view._ngl_component_names)

            if self.widget_accordion_repr_parameters.selected_index >= 0:
                self.widget_repr_parameters.name = self.widget_repr_name.value
                self.widget_repr_parameters.repr_index = self.widget_repr_slider.value
                self.widget_repr_parameters.component_index = self.widget_component_slider.value

        def on_repr_selection_value_changed(change):
            if self._real_time_update:
                component = self.widget_component_slider.value
                repr_index = self.widget_repr_slider.value
                self._view._set_selection(change['new'],
                                          component=component,
                                          repr_index=repr_index)

        def on_change_component_dropdown(change):
            choice = change['new']
            if choice and choice.strip():  # bool(" ") is True
                self.widget_component_slider.value = self._view._ngl_component_names.index(
                    choice)

        self.widget_component_dropdown.observe(on_change_component_dropdown,
                                               names='value')

        self.widget_repr_slider.observe(
            on_component_or_repr_slider_value_changed, names='value')
        self.widget_component_slider.observe(
            on_component_or_repr_slider_value_changed, names='value')
        self.widget_repr_name.observe(on_repr_name_text_value_changed,
                                      names='value')
        repr_selection.observe(on_repr_selection_value_changed, names='value')

        self.widget_repr_control_buttons = self._make_button_repr_control(
            self.widget_component_slider, self.widget_repr_slider,
            repr_selection)

        blank_box = Box([Label("")])

        all_kids = [
            self.widget_repr_control_buttons, blank_box, self.widget_repr_add,
            self.widget_component_dropdown, self.widget_repr_name,
            repr_selection, self.widget_component_slider,
            self.widget_repr_slider, self.widget_repr_choices,
            self.widget_accordion_repr_parameters
        ]

        vbox = VBox(all_kids)

        self._view._request_repr_parameters(
            component=self.widget_component_slider.value,
            repr_index=self.widget_repr_slider.value)

        self.widget_repr = _relayout_master(vbox, width='100%')

        self._refresh(self.widget_component_slider, self.widget_repr_slider)

        setattr(self.widget_repr, "_saved_widgets", [])
        for _box in self.widget_repr.children:
            if hasattr(_box, 'children'):
                for kid in _box.children:
                    self.widget_repr._saved_widgets.append(kid)

        return self.widget_repr

    def _make_widget_repr_parameters(self,
                                     component_slider,
                                     repr_slider,
                                     repr_name_text=None):
        name = repr_name_text.value if repr_name_text is not None else ' '
        widget = self._view._display_repr(component=component_slider.value,
                                          repr_index=repr_slider.value,
                                          name=name)
        widget._ngl_name = 'repr_parameters_box'
        return widget

    def _make_button_export_image(self):
        slider_factor = IntSlider(value=4, min=1, max=10, description='scale')
        checkbox_antialias = Checkbox(value=True, description='antialias')
        checkbox_trim = Checkbox(value=False, description='trim')
        checkbox_transparent = Checkbox(value=False, description='transparent')
        filename_text = Text(value='Screenshot', description='Filename')
        delay_text = FloatText(value=1,
                               description='delay (s)',
                               tooltip='hello')

        start_text, stop_text, step_text = (IntText(value=0,
                                                    description='start'),
                                            IntText(value=self._view.max_frame+1,
                                                    description='stop'),
                                            IntText(value=1,
                                                    description='step'))

        start_text.layout.max_width = stop_text.layout.max_width = step_text.layout.max_width \
                = filename_text.layout.max_width = delay_text.layout.max_width = default.DEFAULT_TEXT_WIDTH

        button_movie_images = Button(description='Export Images')

        def download_image(filename):
            self._view.download_image(factor=slider_factor.value,
                                      antialias=checkbox_antialias.value,
                                      trim=checkbox_trim.value,
                                      transparent=checkbox_transparent.value,
                                      filename=filename)

        @button_movie_images.on_click
        def on_click_images(button_movie_images):
            for i in range(start_text.value, stop_text.value, step_text.value):
                self._view.frame = i
                time.sleep(delay_text.value)
                download_image(filename=filename_text.value + str(i))
                time.sleep(delay_text.value)

        vbox = VBox([
            button_movie_images,
            start_text,
            stop_text,
            step_text,
            delay_text,
            filename_text,
            slider_factor,
            checkbox_antialias,
            checkbox_trim,
            checkbox_transparent,
        ])

        form_items = _relayout(vbox, make_form_item_layout())
        form = Box(form_items, layout=_make_box_layout())
        # form = _relayout_master(vbox)
        return form

    def _make_resize_notebook_slider(self):
        resize_notebook_slider = IntSlider(min=300,
                                           max=2000,
                                           description='resize notebook')

        def on_resize_notebook(change):
            width = change['new']
            self._view._remote_call('resizeNotebook',
                                    target='Widget',
                                    args=[
                                        width,
                                    ])

        resize_notebook_slider.observe(on_resize_notebook, names='value')
        return resize_notebook_slider

    def _make_add_widget_repr(self, component_slider):
        dropdown_repr_name = Dropdown(options=REPRESENTATION_NAMES,
                                      value='cartoon')
        repr_selection = Text(value='*', description='')
        repr_button = Button(description='Add',
                             tooltip="""Add representation.
        You can also hit Enter in selection box""")
        repr_button.layout = Layout(width='auto', flex='1 1 auto')

        dropdown_repr_name.layout.width = repr_selection.layout.width = default.DEFAULT_TEXT_WIDTH

        def on_click_or_submit(button_or_text_area):
            self._view.add_representation(
                selection=repr_selection.value.strip(),
                repr_type=dropdown_repr_name.value,
                component=component_slider.value)

        repr_button.on_click(on_click_or_submit)
        repr_selection.on_submit(on_click_or_submit)
        add_repr_box = HBox([repr_button, dropdown_repr_name, repr_selection])
        add_repr_box._ngl_name = 'add_repr_box'

        return add_repr_box

    def _make_repr_playground(self):
        vbox = VBox()
        children = []

        rep_names = REPRESENTATION_NAMES[:]
        excluded_names = ['ball+stick', 'distance']
        for name in excluded_names:
            rep_names.remove(name)

        repr_selection = Text(value='*')
        repr_selection.layout.width = default.DEFAULT_TEXT_WIDTH
        repr_selection_box = HBox([Label('selection'), repr_selection])
        setattr(repr_selection_box, 'value', repr_selection.value)

        for index, name in enumerate(rep_names):
            button = ToggleButton(description=name)

            def make_func():
                def on_toggle_button_value_change(change, button=button):
                    selection = repr_selection.value
                    new = change['new']  # True/False
                    if new:
                        self._view.add_representation(button.description,
                                                      selection=selection)
                    else:
                        self._view._remove_representations_by_name(
                            button.description)

                return on_toggle_button_value_change

            button.observe(make_func(), names='value')
            children.append(button)

        pd = self._view._player_dict
        pd['widget_quick_repr'] = defaultdict(dict)
        for k in children:
            pd['widget_quick_repr'][k.model_id] = {}

            def click_true(k):
                k.value = True

            pd['widget_quick_repr'][k.model_id][True] = _dry_run(self._view,
                                                                 click_true,
                                                                 k=k)

            def click_false(k):
                k.value = False

            pd['widget_quick_repr'][k.model_id][False] = _dry_run(self._view,
                                                                  click_false,
                                                                  k=k)

        button_clear = Button(description='clear',
                              button_style='info',
                              icon='fa-eraser')

        @button_clear.on_click
        def on_clear(button_clear):
            self._view.clear()
            for kid in children:
                # unselect
                kid.value = False

        if hasattr(button_clear, 'click'):
            _dry_run(self._view, button_clear.click)

        vbox.children = children + [repr_selection, button_clear]
        _make_autofit(vbox)
        self.widget_quick_repr = vbox
        return self.widget_quick_repr

    def _make_repr_name_choices(self, component_slider, repr_slider):
        repr_choices = Dropdown(options=[
            " ",
        ])

        def on_chosen(change):
            repr_name = change.get('new', " ")
            try:
                repr_index = repr_choices.options.index(repr_name)
                repr_slider.value = repr_index
            except ValueError:
                pass

        repr_choices.observe(on_chosen, names='value')
        repr_choices.layout.width = default.DEFAULT_TEXT_WIDTH

        self.widget_repr_choices = repr_choices
        return self.widget_repr_choices

    def _make_widget_picked(self):
        self.widget_picked = self._make_text_picked()
        picked_box = HBox([
            self.widget_picked,
        ])
        return _relayout_master(picked_box, width='75%')

    def _make_export_image_widget(self):
        if self.widget_export_image is None:
            self.widget_export_image = HBox([self._make_button_export_image()])
        return self.widget_export_image

    def _make_extra_box(self):
        if self.widget_extra is None:
            extra_list = [(self._make_widget_picked, 'Picked'),
                          (self._make_repr_playground, 'Quick'),
                          (self._make_export_image_widget, 'Image'),
                          (self._make_command_box, 'Command')]

            self.widget_extra = _make_delay_tab(extra_list, selected_index=0)
        return self.widget_extra

    def _make_general_box(self):
        if self.widget_general is None:
            step_slide = IntSlider(value=self.step,
                                   min=-100,
                                   max=100,
                                   description='step')
            delay_text = IntSlider(value=self.delay,
                                   min=10,
                                   max=1000,
                                   description='delay')
            toggle_button_interpolate = ToggleButton(
                self.interpolate,
                description='Smoothing',
                tooltip='smoothing trajectory')
            link((toggle_button_interpolate, 'value'), (self, 'interpolate'))

            self.widget_background = background_color_picker = ColorPicker(
                value=self._view._ngl_full_stage_parameters.get(
                    'backgroundColor', 'white'),
                description='background')
            self.widget_camera = camera_type = Dropdown(
                value=self.camera,
                options=['perspective', 'orthographic'],
                description='camera')

            link((step_slide, 'value'), (self, 'step'))
            link((delay_text, 'value'), (self, 'delay'))
            link((toggle_button_interpolate, 'value'), (self, 'interpolate'))
            link((camera_type, 'value'), (self, 'camera'))
            link((background_color_picker, 'value'),
                 (self._view, 'background'))

            center_button = self._make_button_center()
            render_button = self._show_download_image()
            center_render_hbox = _make_autofit(
                HBox([
                    toggle_button_interpolate,
                    center_button,
                    render_button,
                ]))

            v0_left = VBox([
                step_slide,
                delay_text,
                background_color_picker,
                camera_type,
                center_render_hbox,
            ])

            v0_left = _relayout_master(v0_left, width='100%')
            self.widget_general = v0_left
        return self.widget_general

    def _ipython_display_(self, **kwargs):
        display(self._display())

    def _make_command_box(self):
        widget_text_command = Text()

        @widget_text_command.on_submit
        def _on_submit_command(_):
            command = widget_text_command.value
            js_utils.execute(command)
            widget_text_command.value = ''

        return widget_text_command

    def _create_all_tabs(self):
        tab = self._display()
        for index, _ in enumerate(tab.children):
            # trigger ceating widgets
            tab.selected_index = index

        self.widget_extra = self._make_extra_box()
        for index, _ in enumerate(self.widget_extra.children):
            self.widget_extra.selected_index = index

    def _simplify_repr_control(self):
        for widget in self.widget_repr._saved_widgets:
            if not isinstance(widget, Tab):
                widget.layout.display = 'none'
        self.widget_repr_choices.layout.display = 'flex'
        self.widget_accordion_repr_parameters.selected_index = 0

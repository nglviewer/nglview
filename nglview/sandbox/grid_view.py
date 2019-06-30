from uuid import uuid4

from ipywidgets import Button, GridBox, HBox, Layout  # ipywidgets >= 7.3

from nglview import NGLWidget
from nglview.utils import js_utils

_code_set_size = """
var ww = window.outerWidth
var wh = window.outerHeight
var vw = ww / %s + 'px'
var vh = wh / %s + 50 + 'px'
var pmodel =  this.model.widget_manager.get_model('%s')

pmodel.then(function(model){
    var key = Object.keys(model.views)[0];
    model.views[key].then(function(v){
        v.children_views.views.forEach(function(pv){
             pv.then(function(v){
                 v.setSize(vw, vh);
             })
        })
    })
})
"""

_code_esc_callback = """
document.onkeydown = function(event){
    if (event.keyCode === 27){
        pmodel.then(function(model){
            var key = Object.keys(model.views)[0];
            model.views[key].then(function(v){
                v.children_views.views.forEach(function(pv){
                    pv.then(function(pvv){
                        pvv.setSize('400px', '300px');
                    })
                })
            })
        })
    }
}
"""


class GridBoxNGL(GridBox):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._ngl_class_id = None

    def add_class(self, class_id):
        super().add_class(class_id)
        self._ngl_class_id = class_id
        return self

    def fullscreen(self, js_code=None):
        code = (f"""
            var ele = document.getElementsByClassName('{self._ngl_class_id}')[0];
            ele.requestFullscreen();
        """)
        n_rows = len(self.children) // self._n_columns
        n_columns = self._n_columns
        code_fullscreen = _code_set_size % (n_columns, n_rows,
                                            self.model_id) + _code_esc_callback
        self.children[0]._execute_js_code(code_fullscreen)

        if js_code is not None:
            code += js_code
        js_utils.run(code)


class GridBoxViewAndPlayer(GridBoxNGL):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._view = self.children[0]
        self._player = self._view.player
        self._player._create_all_widgets()
        self._add_button_fullscreen()

    def _add_button_fullscreen(self):
        btn = Button(description='Fullscreen', icon='fullscreen')

        def on_click(b):
            self.fullscreen()

        btn.on_click(on_click)

        self._player.widget_general.children = list(
            self._player.widget_general.children) + [btn]


# def _sync_camera_pair(v0, v1):
#     v0._set_sync_camera(v0
#     def on_change_camera(change):
#         new = change['new']
#         owner = change['owner']
#         if owner == v0 and v0._ngl_focused:
#             v1.control.orient(v0._camera_orientation)
#         elif owner == v1 and v1._ngl_focused:
#             v0.control.orient(v1._camera_orientation)
#
#     v0.observe(on_change_camera, '_camera_orientation')
#     v1.observe(on_change_camera, '_camera_orientation')

# def _sync_all(views):
#     from itertools import combinations
#     views = [v for v in views if isinstance(v, NGLWidget)]
#     for v0, v1 in combinations(views, 2):
#         _sync_camera_pair(v0, v1)


def _sync_all(views):
    views = {v for v in views if isinstance(v, NGLWidget)}
    for v in views:
        v._set_sync_camera(views)


def grid_view(views,
              n_columns,
              grid_class=GridBoxViewAndPlayer,
              fullscreen=False,
              sync_camera=False,
              **grid_kwargs):
    if sync_camera:
        _sync_all(views)

    grid_template_columns = f"{str(100 / n_columns)}% " * n_columns
    grid_kwargs = grid_kwargs or {
        'grid_template_columns': grid_template_columns,
        'grid_spacing': '0px 0px'
    }
    box = grid_class(views, layout=Layout(**grid_kwargs))
    class_id = f"nglview-grid-{str(uuid4())}"
    box.add_class(class_id)
    box._n_columns = n_columns

    if fullscreen:

        def on_displayed(box):
            box.fullscreen()

        box.on_displayed(on_displayed)
    return box


def fullscreen_mode(view):
    player = view.player
    player._create_all_widgets()
    b = player._display()
    b.layout.width = '400px'
    b.layout.align_self = 'stretch'
    bb = GridBox([view, b], layout=Layout(grid_template_columns='70% 30%'))
    class_id = f'nglview-{uuid4()}'
    bb.add_class(class_id)

    btn = Button(description='Fullscreen', icon='fullscreen')

    def on_click(b):
        view._execute_js_code("""
        var ele = document.getElementsByClassName("%s")[0];

        this.stage.viewer.container.style.height = window.screen.height + 'px'
        console.log(window.innerHeight)
        this.stage.toggleFullscreen(ele);
        """ % class_id)

    btn.on_click(on_click)
    player.widget_general.children = list(
        player.widget_general.children) + [btn]
    return bb

from nglview.utils import js_utils
from ipywidgets import GridBox, Layout # ipywidgets >= 7.3
from uuid import uuid4


_code_fullscreen = """
var ww = window.outerWidth
var wh = window.outerHeight
var vw = ww / %s + 'px'
var vh = wh / %s + 50 + 'px'
var pmodel =  this.model.widget_manager.get_model('%s')
console.log(pmodel);

pmodel.then(function(model){
    var key = Object.keys(model.views)[0];
    model.views[key].then(function(v){
        console.log(v.el.style)
        console.log(v.el.style)
        v.children_views.views.forEach(function(pv){
             pv.then(function(v){
                 v.setSize(vw, vh);
             })
        })
    })
})

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
        """           )
        n_rows = len(self.children) // self._n_columns
        n_columns = self._n_columns
        code_fullscreen = _code_fullscreen % (n_columns, n_rows, self.model_id)
        self.children[0]._execute_js_code(code_fullscreen)

        if js_code is not None:
            code += js_code
        js_utils.run(code)


def _sync_camera_pair(v0, v1):
    def on_change_camera(change):
        new = change['new']
        owner = change['owner']
        if owner == v0 and v0._ngl_focused:
            v1.control.orient(v0._camera_orientation)
        elif owner == v1 and v1._ngl_focused:
            v0.control.orient(v1._camera_orientation)

    v0.observe(on_change_camera, '_camera_orientation')
    v1.observe(on_change_camera, '_camera_orientation')


def _sync_all(views):
    from itertools import combinations
    for v0, v1 in combinations(views, 2):
        _sync_camera_pair(v0, v1)


def grid_view(views, n_columns, fullscreen=False, sync_camera=True):
    if sync_camera:
        _sync_all(views)

    grid_template_columns = f"{str(100 / n_columns)}% " * n_columns
    box = GridBoxNGL(views,
            layout=Layout(
                grid_template_columns=grid_template_columns,
                grid_spacing='0px 0px',
                ))
    class_id = f"nglview-grid-{str(uuid4())}"
    box.add_class(class_id)
    box._n_columns = n_columns

    if fullscreen:
        def on_displayed(box):
            box.fullscreen()

        box.on_displayed(on_displayed)
    return box

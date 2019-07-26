var Jupyter
var widgets = require("@jupyter-widgets/base")
var NGL = require('ngl')
require('./gui.js') // NGL.SidebarWidget


var SidebarModel = widgets.DOMWidgetModel.extend({
    defaults: function(){
        return _.extend(widgets.DOMWidgetModel.prototype.defaults(), {
            _model_name: 'SidebarModel',
            _model_module: 'nglview-js-widgets',
            _model_module_version: require("../package.json").version,
            _view_name: "SidebarView",
            _view_module: "nglview-js-widgets",
            _view_module_version: require("../package.json").version,
        });
    }
})

var SidebarView = widgets.DOMWidgetView.extend({
    render: function() {
        this.sidebar = new NGL.SidebarWidget()
        this.sidebar.container
        this.el.appendChild(this.sidebar.container.dom)
    },

    createView: function(){
        var target_model_id = this.model.get("_target_model_id")
        return this.model.widget_manager.get_model(target_model_id).then((model) =>{
            var key = Object.keys(model.views)[0]
            return model.views[key].then((view) => {
                return new NGL.SidebarWidget(view.stage)
            })
        })
    },
})


module.exports = {
    'SidebarModel': SidebarModel,
    'SidebarView': SidebarView,
};

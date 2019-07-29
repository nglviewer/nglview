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
        this.sidebar = undefined
    },

    setStage: function(stage){
        // stage: NGL.Stage
        if (this.sidebar){
            this.sidebar.dispose()
            this.sidebar = undefined
        }
        this.sidebar = NGL.SidebarWidget(stage)
        this.el.appendChild(this.sidebar.dom)
    }
})


module.exports = {
    'SidebarModel': SidebarModel,
    'SidebarView': SidebarView,
};

import * as widgets from "@jupyter-widgets/base"
var NGL = require('ngl')
require('./gui.js') // NGL.SidebarWidget


export
class SidebarModel extends widgets.DOMWidgetModel{
    defaults(){
        return {
            ...super.defaults(),
            _model_name: 'SidebarModel',
            _model_module: 'nglview-js-widgets',
            _model_module_version: require("../package.json").version,
            _view_name: "SidebarView",
            _view_module: "nglview-js-widgets",
            _view_module_version: require("../package.json").version,
        }
    }
}

export
class SidebarView extends widgets.DOMWidgetView{
    render: function() {
        this.sidebar = undefined
    }

    setStage: function(stage, view=undefined){
        // stage: NGL.Stage
        // view: NGLView
        if (this.sidebar){
            this.sidebar.dispose()
            this.sidebar = undefined
        }
        this.sidebar = NGL.SidebarWidget(stage, view)
        this.el.appendChild(this.sidebar.dom)
    }
}

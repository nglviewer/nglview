import * as _ from 'underscore'
import { BaseView } from "./base"
import * as widgets from "@jupyter-widgets/base"


export
class ThemeManagerModel extends widgets.DOMWidgetModel {
    defaults(){
        return _.extend(widgets.DOMWidgetModel.prototype.defaults(), {
            _model_name: 'ThemeManagerModel',
            _model_module: 'nglview-js-widgets',
            _model_module_version: require("../package.json").version,
            _view_name: "ThemeManagerView",
            _view_module: "nglview-js-widgets",
            _view_module_version: require("../package.json").version,
        })
    }
}


export
class ThemeManagerView extends BaseView {

    render(){
        super.render()
        if (this.isEmbeded()){
            console.log("Embed mode for ThemeManagerView")
            this.handleEmbed()
        }
    }

    handleEmbed(){
        // should be called later by NGLView
        // to make sure theme is setup fisrt.
        this.on_msg(this.model.get("_msg_ar").pop())
    }

    setTheme(cssContent: string){
        var ele = document.getElementById('nglview_style')
        if (ele){
            document.head.removeChild(ele)
        }
        ele = document.createElement('style')
        ele.id = 'nglview_style'
        ele.appendChild(document.createTextNode(cssContent))
        document.head.appendChild(ele)
    }
}

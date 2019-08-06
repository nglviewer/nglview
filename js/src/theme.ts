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
        this.displayed.then(() => this.handleThemeChanged())
    }

    handleThemeChanged(){
        console.log("ThemeManagerView: handleThemeChanged")
        this.setTheme(this.model.get("_theme_css"))
    }

    handleEmbed(){
        this.handleThemeChanged()
    }

    setTheme(cssContent: string){
        if (cssContent == undefined){
            return
        }
        var ele = document.getElementById('nglview_style')
        if (ele){
            document.head.removeChild(ele)
        }
        ele = document.createElement('style')
        ele.id = 'nglview_style'
        ele.appendChild(document.createTextNode(cssContent))
        document.head.appendChild(ele)
        this.send({'type': 'call_method', 'data': 'handle_resize'})
    }
}

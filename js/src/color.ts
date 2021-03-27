import * as _ from 'underscore'
import * as NGL from "ngl"
import { BaseView } from "./base"
import * as widgets from "@jupyter-widgets/base"


export
class ColormakerRegistryModel extends widgets.DOMWidgetModel {
    defaults(){
        return _.extend(widgets.DOMWidgetModel.prototype.defaults(), {
            _model_name: 'ColormakerRegistryModel',
            _model_module: 'nglview-js-widgets',
            _model_module_version: require("../package.json").version,
            _view_name: "ColormakerRegistryView",
            _view_module: "nglview-js-widgets",
            _view_module_version: require("../package.json").version,
        })
    }
}


export
class ColormakerRegistryView extends BaseView {

    render(){
        super.render()
        if (this.isEmbeded()){
            console.log("Embed mode for ColormakerRegistryView")
            this.handleEmbed()
        }
    }

    handleEmbed(){
        // should be called by NGLView later
        this.model.get("_msg_ar").forEach(msg =>{
            this.on_msg(msg)
        })
    }

    addSelectionScheme(label, args){
        var id = NGL.ColormakerRegistry.addSelectionScheme(args, label)
        this._updateId(id, label)
    }

    addSelectionSchemeOriginal(label, args){
        var id = NGL.ColormakerRegistry.addSelectionScheme(args, label);
        var scheme = NGL.ColormakerRegistry.userSchemes[id];
        NGL.ColormakerRegistry.removeScheme(id);
        // hard code the scheme ID
        NGL.ColormakerRegistry.add(label, scheme);
    }

    addScheme(label, func_str){
        var func = Function("return " + func_str)()
        var id = NGL.ColormakerRegistry.addScheme(function(params){
            this.atomColor = func
        })
        this._updateId(id, label)
    }

    _updateId(oldId, newId){
        var scheme = NGL.ColormakerRegistry.userSchemes[oldId]
        NGL.ColormakerRegistry.add(newId, scheme)
        NGL.ColormakerRegistry.removeScheme(oldId)
    }
}

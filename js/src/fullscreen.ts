import * as widgets from "@jupyter-widgets/base"
import * as NGL from 'ngl'
import { BaseView } from "./base"


export
class FullscreenModel extends widgets.DOMWidgetModel{
    defaults(){
        return {
                ...super.defaults(),
                _model_name: 'FullscreenModel',
                _model_module: 'nglview-js-widgets',
                _model_module_version: require("../package.json").version,
                _view_name: "FullscreenView",
                _view_module: "nglview-js-widgets",
                _view_module_version: require("../package.json").version,
            }
        }
}

export 
class ABC{
}


export
class FullscreenView extends BaseView{
    stage : any
    render() {
        this.stage = new NGL.Stage()
        var that = this
        this.model.on("msg:custom", function(msg){
            that.on_msg(msg)
        })
        this.handleSignals()
    }

    fullscreen(model_id){
        var that = this
        this.model.widget_manager.get_model(model_id).then((model) =>{
            var key = Object.keys(model.views)[0]
            model.views[key].then((view) => {
                that.stage.toggleFullscreen(view.el)
            })
        })
    }

    handleSignals(){
        var that = this
        this.stage.signals.fullscreenChanged.add(function (isFullscreen) {
            that.model.set("_is_fullscreen", isFullscreen)
            that.touch()
        })
    }

    executeCode(code){
        eval(code);
    }

    on_msg(msg){
        if ('executeCode' in msg){
            this.executeCode(msg.executeCode)
        }
    }
}

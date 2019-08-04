import * as widgets from "@jupyter-widgets/base"
import * as NGL from "ngl"


export
class BaseView extends widgets.DOMWidgetView {

    render(){
        this.handleMessage();
        this.displayed.then(() =>{
            this.model.set("_ready", true)
            this.touch()
        })
    }

    executeCode(code){
        eval(code);
    }

    handleMessage(){
        this.model.on("msg:custom", function(msg){
            this.on_msg(msg)
        }.bind(this))
    }

    on_msg(msg){
        if (msg.type == 'callMethod'){
            this[msg.methodName].apply(this, msg.args, msg.kwargs)
        }
    }
}

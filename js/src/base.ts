import * as widgets from "@jupyter-widgets/base"
var NGL = require("ngl") // not sure why I can not import it


export
class BaseView extends widgets.DOMWidgetView {

    render(){
    if (this.isEmbeded()){
            // embed mode: let NGLView call the handleEmbed directly
            return
        }
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
            console.log("from BaseView", msg)
            this[msg.methodName].apply(this, msg.args, msg.kwargs)
        }
    }

    isEmbeded(){
        return (this.model.comm == undefined)
    }
}

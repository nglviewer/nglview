var widgets = require("@jupyter-widgets/base")
var NGL = require("ngl")


var BaseView = widgets.DOMWidgetView.extend({

    render: function(){
        this.handleMessage();
        this.displayed.then(() =>{
            this.model.set("_ready", true)
            this.touch()
        })
    },

    executeCode: function(code){
        eval(code);
    },

    handleMessage: function(){
        this.model.on("msg:custom", function(msg){
            this.on_msg(msg)
        }.bind(this))
    },

    on_msg: function(msg){
        if (msg.type == 'callMethod'){
            this[msg.methodName].apply(this, msg.args, msg.kwargs)
        }
    }
})

module.exports = {
    "BaseView": BaseView
}

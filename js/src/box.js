var widgets = require("@jupyter-widgets/base")

var GridBoxNGLModel = widgets.DOMWidgetModel.extend({
    defaults: function(){
        return _.extend(widgets.DOMWidgetModel.prototype.defaults(), {
            _model_name: 'GridBoxNGLModel',
            _model_module: 'nglview-js-widgets',
            _model_module_version: require("../package.json").version,
            _view_name: "GridBoxNGLView",
            _view_module: "nglview-js-widgets",
            _view_module_version: require("../package.json").version,
        });
    }
})

var GridBoxNGLView = widgets.DOMWidgetView.extend({
    render: function() {
        this.model.on("msg:custom", function(msg){
            this.on_msg(msg)
        })
    },

    execute_code: function(code){
        eval(code);
    },

    on_msg: function(msg){
        if ('execute_code' in msg){
            this.execute_code(msg.execute_code)
        }
    }

});

module.exports = {
    'GridBoxNGLView': GridBoxNGLView,
    'GridBoxNGLModel': GridBoxNGLModel,
};

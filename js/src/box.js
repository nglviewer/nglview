var widgets = require("@jupyter-widgets/base")
var NGL = require("ngl")

var GridBoxNGLModel = widgets.GridBoxModel.extend({
    defaults: function(){
        return _.extend(widgets.GridBoxModel.prototype.defaults(), {
            _model_name: 'GridBoxNGLModel',
            _model_module: 'nglview-js-widgets',
            _model_module_version: require("../package.json").version,
            _view_name: "GridBoxNGLView",
            _view_module: "nglview-js-widgets",
            _view_module_version: require("../package.json").version,
        });
    }
})

var GridBoxNGLView = widgets.GridBoxView.extend({
    render: function() {
        this.stage = new NGL.Stage()
        widgets.GridBoxView.prototype.render.call(this)
        var that = this
        this.model.on("msg:custom", function(msg){
            that.on_msg(msg)
        })
        this.handleResize()
        this.handleSignals()
    },

    handleSignals: function(){
        var that = this
        this.stage.signals.fullscreenChanged.add(function (isFullscreen) {
            this.handleResize()
            if (!isFullscreen){
                this.el.style.height = '300px' // FIXME: record previous height?
            }
        })
    },

    handleResize: function(){
        var that = this
        this.children_views.views.forEach((view)  => {
            view.then((view) => {
                var box = that.el.getBoundingClientRect()
                var w = box.width / 2
                w = w + 'px'
                var h = box.height + 'px'
                view.setSize(w, h)
            })
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

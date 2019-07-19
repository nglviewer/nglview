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
        this.displayed.then(() =>{
            var size = this.getSize()
            // FIXME: got "0px" for height after displaying
            console.log("size")
            console.log(size)
            this.model.set("width", size[0])
            this.touch()
            this.model.set("height", '300px')
            this.touch()
        })
    },

    handleSignals: function(){
        var that = this
        this.stage.signals.fullscreenChanged.add(function (isFullscreen) {
            that.handleResize()
            if (!isFullscreen){
                console.log("not isFullscreen")
                that.setSize(
                    that.model.get("width"),
                    that.model.get("height"))
            }
        })
    },

    getSize: function(){
        var box = this.el.getBoundingClientRect()
        return [box.width + 'px', box.height + 'px']
    },

    setSize: function(w, h){
        // px
        console.log('setSize')
        console.log(w, h)
        this.el.style.width = w
        this.el.style.height = h
        this.handleResize()
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

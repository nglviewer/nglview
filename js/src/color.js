var NGL = require("ngl")
var BaseView = require("./base").BaseView
var widgets = require("@jupyter-widgets/base")


var ColormakerRegistryModel = widgets.DOMWidgetModel.extend({
    defaults: function(){
        return _.extend(widgets.DOMWidgetModel.prototype.defaults(), {
            _model_name: 'ColormakerRegistryModel',
            _model_module: 'nglview-js-widgets',
            _model_module_version: require("../package.json").version,
            _view_name: "ColormakerRegistryView",
            _view_module: "nglview-js-widgets",
            _view_module_version: require("../package.json").version,
        });
    }
})


var ColormakerRegistryView = BaseView.extend({
    addSelectionScheme: function(label, args){
        console.log('label', + label)
        console.log(args)
        var id = NGL.ColormakerRegistry.addSelectionScheme(args, label)
        var scheme = NGL.ColormakerRegistry.userSchemes[id]
        // hard code the scheme ID
        NGL.ColormakerRegistry.removeScheme(id)
        NGL.ColormakerRegistry.add(label, scheme)
        return label
    },

    addScheme: function(func_str){
        var func = eval(func_str)
        var schemeId = NGL.ColormakerRegistry.addScheme(func)
    }

    removeScheme: function(schemeId){
        NGL.ColormakerRegistry.removeScheme(schemeId)
    }
})


module.exports = {
    ColormakerRegistryView: ColormakerRegistryView,
    ColormakerRegistryModel: ColormakerRegistryModel
}

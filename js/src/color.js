var _ = require('underscore')
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
        var id = NGL.ColormakerRegistry.addSelectionScheme(args, label)
        this._updateId(id, label)
    },

    addSelectionSchemeOriginal: function(label, args){
        var id = NGL.ColormakerRegistry.addSelectionScheme(args, label);
        var scheme = NGL.ColormakerRegistry.userSchemes[id];
        NGL.ColormakerRegistry.removeScheme(id);
        // hard code the scheme ID
        NGL.ColormakerRegistry.add(label, scheme);
    },

    addScheme: function(label, func_str){
        var func = Function("return " + func_str)()
        console.log(func)
        var id = NGL.ColormakerRegistry.addScheme(function(params){
            this.atomColor = func
        })
        this._updateId(id, label)
    },

    _updateId: function(oldId, newId){
        var scheme = NGL.ColormakerRegistry.userSchemes[oldId]
        console.log(oldId, scheme)
        NGL.ColormakerRegistry.add(newId, scheme)
        NGL.ColormakerRegistry.removeScheme(oldId)
    },

    removeScheme: function(schemeId){
        NGL.ColormakerRegistry.removeScheme(schemeId)
    },
})


module.exports = {
    ColormakerRegistryView: ColormakerRegistryView,
    ColormakerRegistryModel: ColormakerRegistryModel
}

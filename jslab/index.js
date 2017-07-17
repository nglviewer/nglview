
var nglview_js_widgets = require('nglview-js-widgets');
var jupyterlab_widgets = require('@jupyter-widgets/jupyterlab-manager');


module.exports = {
  id: 'nglview-js-widgets',
  requires: [jupyterlab_widgets.INBWidgetExtension],
  activate: function(app, widgets) {
      widgets.registerWidget({
          name: 'nglview-js-widgets',
          version: nglview_js_widgets.version,
          exports: nglview_js_widgets
      });
    },
  autoStart: true
};

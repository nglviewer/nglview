var nglview_js_widgets = require('./index');
var base = require('@jupyter-widgets/base');


module.exports = {
  id: 'nglview-js-widgets',
  requires: [base.IJupyterWidgetRegistry],
  activate: function(app, widgets) {
      widgets.registerWidget({
          name: 'nglview-js-widgets',
          version: nglview_js_widgets.version,
          exports: nglview_js_widgets
      });
    },
  autoStart: true
};

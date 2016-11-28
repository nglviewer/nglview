var nglview = require('./index');

var jupyterlab_widgets = require('@jupyterlab/nbwidgets');

/**
 * The widget manager provider.
 */
module.exports = {
  id: 'jupyter.extensions.nglview',
  requires: [jupyterlab_widgets.INBWidgetExtension],
  activate: function(app, widgets) {
      widgets.registerWidget({
          name: 'nglview',
          version: nglview.version,
          exports: nglview
      });
    },
  autoStart: true
};

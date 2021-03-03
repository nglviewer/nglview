import _ from 'underscore';
import NGL from 'ngl';
import { BaseView } from './base';
import widgets from '@jupyter-widgets/base';
import pkg = require('../package.json');

export class ColormakerRegistryModel extends widgets.DOMWidgetModel {
  defaults() {
    return _.extend(widgets.DOMWidgetModel.prototype.defaults(), {
      _model_name: 'ColormakerRegistryModel',
      _model_module: 'nglview-js-widgets',
      _model_module_version: pkg.version,
      _view_name: 'ColormakerRegistryView',
      _view_module: 'nglview-js-widgets',
      _view_module_version: pkg.version
    });
  }
}

export class ColormakerRegistryView extends BaseView {
  render() {
    super.render();
    if (this.isEmbeded()) {
      console.log('Embed mode for ColormakerRegistryView');
      this.handleEmbed();
    }
  }

  handleEmbed() {
    // should be called by NGLView later
    this.model.get('_msg_ar').forEach(msg => {
      this.on_msg(msg);
    });
  }

  addSelectionScheme(label, args) {
    const id = NGL.ColormakerRegistry.addSelectionScheme(args, label);
    this._updateId(id, label);
  }

  addSelectionSchemeOriginal(label, args) {
    const id = NGL.ColormakerRegistry.addSelectionScheme(args, label);
    const scheme = NGL.ColormakerRegistry.userSchemes[id];
    NGL.ColormakerRegistry.removeScheme(id);
    // hard code the scheme ID
    NGL.ColormakerRegistry.add(label, scheme);
  }

  addScheme(label, func_str) {
    const func = Function('return ' + func_str)();
    const id = NGL.ColormakerRegistry.addScheme(function (params) {
      this.atomColor = func;
    });
    this._updateId(id, label);
  }

  _updateId(oldId, newId) {
    const scheme = NGL.ColormakerRegistry.userSchemes[oldId];
    NGL.ColormakerRegistry.add(newId, scheme);
    NGL.ColormakerRegistry.removeScheme(oldId);
  }
}

import widgets from '@jupyter-widgets/base';
import NGL from 'ngl';
import { BaseView } from './base';
import pkg = require('../package.json');

export class FullscreenModel extends widgets.DOMWidgetModel {
  defaults() {
    return {
      ...super.defaults(),
      _model_name: 'FullscreenModel',
      _model_module: 'nglview-js-widgets',
      _model_module_version: pkg.version,
      _view_name: 'FullscreenView',
      _view_module: 'nglview-js-widgets',
      _view_module_version: pkg.version
    };
  }
}

export class ABC {}

export class FullscreenView extends BaseView {
  stage: any;
  render() {
    this.stage = new NGL.Stage();
    const that = this;
    this.model.on('msg:custom', msg => {
      that.on_msg(msg);
    });
    this.handleSignals();
  }

  fullscreen(model_id) {
    const that = this;
    this.model.widget_manager.get_model(model_id).then(model => {
      const key = Object.keys(model.views)[0];
      model.views[key].then(view => {
        that.stage.toggleFullscreen(view.el);
      });
    });
  }

  handleSignals() {
    const that = this;
    this.stage.signals.fullscreenChanged.add(isFullscreen => {
      that.model.set('_is_fullscreen', isFullscreen);
      that.touch();
    });
  }

  executeCode(code) {
    eval(code);
  }

  on_msg(msg) {
    if ('executeCode' in msg) {
      this.executeCode(msg.executeCode);
    }
  }
}

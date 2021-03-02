import {
  JupyterFrontEnd,
  JupyterFrontEndPlugin
} from '@jupyterlab/application';

/**
 * Initialization data for the nglview-js-widgets extension.
 */
const extension: JupyterFrontEndPlugin<void> = {
  id: 'nglview-js-widgets:plugin',
  autoStart: true,
  activate: (app: JupyterFrontEnd) => {
    console.log('JupyterLab extension nglview-js-widgets is activated!');
  }
};

export default extension;

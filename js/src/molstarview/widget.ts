import * as widgets from '@jupyter-widgets/base';
import * as _ from 'underscore';
import { PluginConfig } from 'molstar/lib/mol-plugin/config';
import { createPluginUI } from 'molstar/lib/mol-plugin-ui';
import * as molStructure from 'molstar/lib/mol-plugin-state/actions/structure';
import { BuiltInTrajectoryFormat } from 'molstar/lib/mol-plugin-state/formats/trajectory';
import { PluginCommands } from 'molstar/lib/mol-plugin/commands';
// require('molstar/lib/mol-plugin-ui/skin/light.scss'); // FIXME: loader issue for labextension building.
import * as representation from "./representation";

// import { basicSpec } from "./ui"

// See example.py for the kernel counterpart to this file.


// Custom Model. Custom widgets models must at least provide default values
// for model attributes, including
//
//  - `_view_name`
//  - `_view_module`
//  - `_view_module_version`
//
//  - `_model_name`
//  - `_model_module`
//  - `_model_module_version`
//
//  when different from the base class.

// When serialiazing the entire widget state for embedding, only values that
// differ from the defaults will be specified.
var MolstarModel = widgets.DOMWidgetModel.extend({
    defaults: _.extend(widgets.DOMWidgetModel.prototype.defaults(), {
        _model_name: 'MolstarModel',
        _view_name: 'MolstarView',
        _model_module: 'molstarview-widget',
        _view_module: 'molstarview-widget',
        _model_module_version: '0.1.0', // FIXME: require('../package.json').version,
        _view_module_version: '0.1.0', // FIXME: require('../package.json').version,
    })
});


// Custom View. Renders the widget model.
var MolstarView = widgets.DOMWidgetView.extend({
    // Defines how the widget gets rendered into the DOM
    async render() {
        this.handleMessage();
        this.displayed.then(async () => {
            await this.initializeDisplay();
            if (this.model.comm == undefined) {
                this.handleEmbed();
            }
            await this.finalizeDisplay();
        });
    },

    async initializeDisplay() {
        this.setupContainer();
        this.plugin = await createPluginUI(this.container);
        this._focused = false;
        await this.checkLeaderView();
    },

    setupContainer() {
        const container = document.createElement('div');
        container.style.width = '800px';
        container.style.height = '600px';
        this.el.appendChild(container);
        this.container = container;
    },

    async checkLeaderView() {
        console.log('Find a leader view');
        this.isLeader = this.model.views.length < 2;
        var hasLeader = false;

        for (var k in this.model.views) {
            var view = await this.model.views[k];
            if (view.isLeader) {
                hasLeader = true;
                break;
            }
        }

        if (!hasLeader) {
            for (var k in this.model.views) {
                var view = await this.model.views[k];
                view.isLeader = true;
                hasLeader = true;
                break;
            }
        }

        if (!this.isLeader) {
            for (var k in this.model.views) {
                var view = await this.model.views[k];
                if (view.isLeader) {
                    var data = await view.plugin.state.getSnapshot();
                    await this.plugin.state.setSnapshot(data);
                    break;
                }
            }
        }
    },

    handleSignals() {
        var that = this;
        this.container.addEventListener('mouseover', function(e: any) {
            that._focused = 1;
            e; // linter
            that.mouseOverDisplay('block');
        }, false);

        this.container.addEventListener('mouseout', function(e: any) {
            that._focused = 0;
            e; // linter
            that.mouseOverDisplay('none');
        }, false);
    },

    async finalizeDisplay() {
        this.send({
            'type': 'request_loaded',
            'data': true
        });
    },

    // from molstar: https://github.com/molstar/molstar/blob/d1e17785b8404eec280ad04a6285ad9429c5c9f3/src/apps/viewer/app.ts#L219-L223
    async loadStructureFromData(
        data: string | number[],
        format: BuiltInTrajectoryFormat,
        preset?: any,
        options?: { dataLabel?: string }
    ) {
        const _data = await this.plugin.builders.data.rawData({ data, label: options?.dataLabel });
        const trajectory = await this.plugin.builders.structure.parseTrajectory(_data, format);

        if (preset) {
            console.log("Calling loadStructureFromData with preset", preset);
            await this.plugin.builders.structure.hierarchy.applyPreset(trajectory, preset);
        } else {
            console.log('Calling loadStructureFromData without preset');
            await this.plugin.builders.structure.createModel(trajectory);
        }
    },

    // from molstar: https://github.com/molstar/molstar/blob/d1e17785b8404eec280ad04a6285ad9429c5c9f3/src/apps/viewer/app.ts#L219-L223
    // this method is taken from the Viewer class
    loadPdb(pdb: any) {
        const params = molStructure.DownloadStructure.createDefaultParams(this.plugin.state.data.root.obj, this.plugin);
        const provider = this.plugin.config.get(PluginConfig.Download.DefaultPdbProvider);
        return this.plugin.runTask(this.plugin.state.data.applyAction(molStructure.DownloadStructure, {
            source: {
                name: 'pdb',
                params: {
                    provider: {
                        id: pdb,
                        server: {
                            name: provider,
                            params: molStructure.PdbDownloadProvider[provider as keyof typeof molStructure.PdbDownloadProvider].defaultValue
                        }
                    },
                    options: { ...params.source.params.options },
                }
            }
        }));
    },

    executeCode(code: any) {
        eval(code);
    },

    on_msg(msg: any) {
        if (msg.type == 'call_method') {
            var new_args = msg.args.slice();
            new_args.push(msg.kwargs);

            switch (msg.target) {
                case 'Widget':
                    var func = this[msg.methodName];
                    if (func) {
                        func.apply(this, new_args);
                    } else {
                        console.log('Can not create func for ' + msg.methodName);
                    }
                    break;
            }
        } else if (msg.type == 'binary_single') {
            this.handleBinaryMessage(msg);
        }
    },

    handleBinaryMessage(msg: any) {
        var coordinateMeta = msg.data;
        var coordinates;
        var keys = Object.keys(coordinateMeta);

        for (var i = 0; i < keys.length; i++) {
            var traj_index = keys[i];
            coordinates = new Float32Array(msg.buffers[i].buffer);
            if (coordinates.byteLength > 0) {
                this.updateCoordinates(coordinates, traj_index);
            }
        }
    },

    handleEmbed() {
        var snaphShot = this.model.get("molstate");
        this.setState(snaphShot);
    },

    handleMessage() {
        this.model.on("msg:custom", (msg: any) => {
            this.on_msg(msg);
        }, this);

        if (this.model.comm) {
            this.model.comm.on_msg((msg: any) => {
                var buffers = msg.buffers;
                var content = msg.content.data.content;
                if (buffers.length && content) {
                    content.buffers = buffers;
                }
                this.model._handle_comm_msg.call(this.model, msg);
            });
        }
    },

    updateCoordinates(coordinates: any, modelIndex: any) {
        var component = 0; // FIXME
        if (coordinates && typeof component != 'undefined') {
            // FIXME: update
        }
    },

    exportImage(modelId: any) {
        this.plugin.helpers.viewportScreenshot.getImageDataUri().then((data: string) => {
            data = data.replace("data:image/png;base64,", "");
            var msg = { "type": "exportImage", "data": data, "model_id": modelId };
            this.send(msg);
        });
    },

    downloadState() {
        PluginCommands.State.Snapshots.DownloadToFile(this.plugin, { type: 'json' });
    },

    async getState() {
        if (this.isLeader) {
            var data = this.plugin.state.getSnapshot();
            this.model.set("molstate", data);
            this.touch();
        }
    },

    async setState(data: any) {
        await this.plugin.state.setSnapshot(data);
    },

    addRepresentation(params: any, modelIndex: any) {
        representation.addRepresentation(this.plugin, params, modelIndex);
    },

    removeRepresentation(modelIndex: any) {
        var st = this.plugin.managers.structure.hierarchy.current.structures[modelIndex];
        this.plugin.managers.structure.component.removeRepresentations(st.components);
    },

    resetCamera() {
        PluginCommands.Camera.Reset(this.plugin, {});
    },

    setCamera(params: any) {
        var durationMs = 0.0;
        this.plugin.canvas3d.requestCameraReset({ durationMs, params });
    },

    getCamera() {
        var snapshot = this.plugin.canvas3d.camera.getSnapshot();
        this.send({ "type": "getCamera", "data": snapshot });
    },

    syncCamera() {
        var that = this;
        if (that._synced_model_ids.length > 0 && that._focused) {
            that._synced_model_ids.forEach(async function(mid: any) {
                var model = await that.model.widget_manager.get_model(mid);
                for (var k in model.views) {
                    var view = await model.views[k];
                    if (view !== that) {
                        view.setCamera(that.plugin.canvas3d.camera.getSnapshot());
                    }
                }
            });
        }
    },
});

module.exports = {
    'MolstarModel': MolstarModel,
    'MolstarView': MolstarView
};
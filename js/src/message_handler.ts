
import { NGLView } from "./widget_ngl";

export class MessageHandler {
    view: NGLView;

    constructor(view: NGLView) {
        this.view = view;
    }

    async on_msg(msg) {
        if (msg.type === 'call_method') {
            await this.handleCallMethod(msg);
        } else if (msg.type === 'base64_single') {
            this.handleBase64Single(msg.data);
        } else if (msg.type === 'binary_single') {
            this.handleBinarySingle(msg);
        } else if (msg.type === 'movie_image_data') {
            this.handleMovieMaking(msg.render_params);
        } else if (msg.type === 'get') {
            this.handleGetRequest(msg.data);
        }
    }

    async handleCallMethod(msg) {
        var index, component, func, stage;
        var new_args = msg.args.slice();
        new_args.push(msg.kwargs);

        if (msg.methodName === 'addRepresentation' && msg.reconstruc_color_scheme) {
            msg.kwargs.color = this.view.addColorScheme(msg.kwargs.color, msg.kwargs.color_label);
        }

        switch (msg.target) {
            case 'Stage':
                await this.handleStageMethod(msg, new_args);
                break;
            case 'Viewer':
                this.view.stage.viewer[msg.methodName].apply(this.view.stage.viewer, new_args);
                break;
            case 'viewerControls':
                this.view.stage.viewerControls[msg.methodName].apply(this.view.stage.viewerControls, new_args);
                break;
            case 'compList':
                component = this.view.stage.compList[msg.component_index];
                component[msg.methodName].apply(component, new_args);
                break;
            case 'Widget':
                this.view[msg.methodName].apply(this.view, new_args);
                break;
            case 'Representation':
                component = this.view.stage.compList[msg.component_index];
                var repr = component.reprList[msg.repr_index];
                repr[msg.methodName].apply(repr, new_args);
                break;
            default:
                console.log('Unknown target: ' + msg.target);
                break;
        }
    }

    async handleStageMethod(msg, new_args) {
        var stage_func = this.view.stage[msg.methodName];
        if (msg.methodName === 'removeComponent') {
            var component = this.view.stage.compList[msg.args[0]];
            this.view.stage.removeComponent(component);
        } else if (msg.methodName === 'loadFile') {
            if (this.view.model.views.length > 1 && msg.kwargs && msg.kwargs.defaultRepresentation) {
                msg.kwargs.defaultRepresentation = false;
            }
            await this.view._handleStageLoadFile(msg);
        } else {
            stage_func.apply(this.view.stage, new_args);
        }
    }

    handleBase64Single(coordinatesDict) {
        var keys = Object.keys(coordinatesDict);
        for (var i = 0; i < keys.length; i++) {
            var traj_index = parseInt(keys[i], 10);
            var coordinates = this.view.decode_base64(coordinatesDict[traj_index]);
            if (coordinates && coordinates.byteLength > 0) {
                this.view.updateCoordinates(coordinates, traj_index);
            }
        }
    }

    handleBinarySingle(msg) {
        var coordinateMeta = msg.data;
        var keys = Object.keys(coordinateMeta);
        for (var i = 0; i < keys.length; i++) {
            var traj_index = parseInt(keys[i], 10);
            var coordinates = new Float32Array(msg.buffers[i].buffer);
            if (coordinates.byteLength > 0) {
                this.view.updateCoordinates(coordinates, traj_index);
            }
        }
        if (msg.movie_making) {
            this.handleMovieMaking(msg.render_params);
        }
    }

    handleGetRequest(data) {
        if (data === 'camera') {
            this.view.send(JSON.stringify(this.view.stage.viewer.camera));
        } else if (data === 'parameters') {
            this.view.send(JSON.stringify(this.view.stage.parameters));
        } else {
            console.log("Number of components", this.view.stage.compList.length);
            console.log("ngl_view_id", this.view.ngl_view_id);
        }
    }

    async handleMovieMaking(render_params) {
        if (this.view.ngl_view_id == this.view.get_last_child_id()) {
            var blob = await this.view.stage.makeImage(render_params);
            var reader = new FileReader();
            var arr_str;
            reader.onload = function () {
                arr_str = (reader.result as string).replace("data:image/png;base64,", "");
                this.view.send({
                    "data": arr_str,
                    "type": "movie_image_data",
                });
                this.view.send({ 'type': 'async_message', 'data': 'ok' });
            }.bind(this);
            reader.readAsDataURL(blob);
        }
    }
}
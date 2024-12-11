import { NGLView } from "./widget_ngl";
import * as NGL from 'ngl';
import * as _ from 'lodash';

export class EmbedHandler {
    view: NGLView;

    constructor(view: NGLView) {
        this.view = view;
    }

    isEmbeded() {
        return (this.view.model.get("_ngl_serialize") || (this.view.model.comm == undefined));
    }

    async handleEmbed() {
        var ngl_msg_archive = this.view.model.get("_ngl_msg_archive");
        var ngl_stage_params = this.view.model.get('_ngl_full_stage_parameters');
        const camera_orientation = this.view.model.get("_camera_orientation");

        if (
            Object.keys(ngl_stage_params).length === 0
            && camera_orientation.length === 0
        ) {
            console.log("No state stored; initializing embedded widget for the first time.");
            for (const msg of ngl_msg_archive) {
                await this.view.on_msg(msg);
            }
            return;
        }

        var loadfile_list = [];

        _.each(ngl_msg_archive, (msg: any) => {
            if (msg.methodName == 'loadFile') {
                if (msg.kwargs && msg.kwargs.defaultRepresentation) {
                    msg.kwargs.defaultRepresentation = false;
                }
                loadfile_list.push(this.view._getLoadFilePromise(msg));
            }
        });

        await Promise.all(loadfile_list);
        this.view.stage.setParameters(ngl_stage_params);
        this.view.set_camera_orientation(camera_orientation);
        this.view.touch();

        if (this.view.model.comm === undefined) {
            var ngl_coordinate_resource = this.view.model.get("_ngl_coordinate_resource");
            var n_frames = ngl_coordinate_resource['n_frames'] || 1;
            this.view.model.set("max_frame", n_frames - 1);
            this.view.touch();
            var model = await this.view.getPlayerModel();
            var pmodel = model.get("children")[0];
            this.view.listenTo(pmodel, "change:value", () => {
                this.view.updateCoordinatesFromDict(ngl_coordinate_resource, pmodel.get("value"));
            });
        }

        for (const msg of this.view.model.get("_ngl_msg_archive")) {
            if (msg.fire_embed) {
                await this.view.on_msg(msg);
            }
        }

        this.view._set_representation_from_repr_dict(this.view.model.get("_ngl_repr_dict"));
        this.view.handleResize();
    }

    _handleEmbedBeforeStage() {
        var ngl_color_dict = this.view.model.get("_ngl_color_dict");
        var label;
        for (label in ngl_color_dict) {
            if (!NGL.ColormakerRegistry.hasScheme(label)) {
                this.view.addColorScheme(ngl_color_dict[label], label);
            }
        }
    }
}
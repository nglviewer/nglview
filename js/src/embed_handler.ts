
import { NGLView } from "./widget_ngl";
import * as NGL from 'ngl';

export class EmbedHandler {
    view: NGLView;

    constructor(view: NGLView) {
        this.view = view;
    }

    isEmbeded() {
        return (this.view.model.get("_ngl_serialize") || (this.view.model.comm == undefined));
    }

    async handleEmbed() {
        var that = this.view;
        var ngl_msg_archive = that.model.get("_ngl_msg_archive");
        var ngl_stage_params = that.model.get('_ngl_full_stage_parameters');
        const camera_orientation = that.model.get("_camera_orientation");

        if (
            Object.keys(ngl_stage_params).length === 0
            && camera_orientation.length === 0
        ) {
            console.log("No state stored; initializing embedded widget for the first time.");
            for (const msg of ngl_msg_archive) {
                await that.on_msg(msg);
            }
            return;
        }

        var loadfile_list = [];

        _.each(ngl_msg_archive, function (msg: any) {
            if (msg.methodName == 'loadFile') {
                if (msg.kwargs && msg.kwargs.defaultRepresentation) {
                    msg.kwargs.defaultRepresentation = false;
                }
                loadfile_list.push(that._getLoadFilePromise(msg));
            }
        });

        await Promise.all(loadfile_list);
        that.stage.setParameters(ngl_stage_params);
        that.set_camera_orientation(camera_orientation);
        that.touch();

        if (that.model.comm === undefined) {
            var ngl_coordinate_resource = that.model.get("_ngl_coordinate_resource");
            var n_frames = ngl_coordinate_resource['n_frames'] || 1;
            that.model.set("max_frame", n_frames - 1);
            that.touch();
            var model = await that.getPlayerModel();
            var pmodel = model.get("children")[0];
            that.listenTo(pmodel, "change:value", function () {
                that.updateCoordinatesFromDict(ngl_coordinate_resource, pmodel.get("value"));
            });
        }

        for (const msg of that.model.get("_ngl_msg_archive")) {
            if (msg.fire_embed) {
                await that.on_msg(msg);
            }
        }

        that._set_representation_from_repr_dict(that.model.get("_ngl_repr_dict"));
        that.handleResize();
    }

    _handleEmbedBeforeStage() {
        var that = this.view;
        var ngl_color_dict = that.model.get("_ngl_color_dict");
        var label;
        for (label in ngl_color_dict) {
            if (!NGL.ColormakerRegistry.hasScheme(label)) {
                that.addColorScheme(ngl_color_dict[label], label);
            }
        }
    }
}
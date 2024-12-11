import * as NGL from "ngl";
import { NGLView } from "./widget_ngl";
import * as $ from "jquery";
import { WidgetView } from "@jupyter-widgets/base";

export class StageManager {
    stage: NGL.Stage;
    view: NGLView;

    constructor(view: NGLView) {
        this.view = view;
    }

    async createStage() {
        var stage_params = {
            ...this.view.model.get("_ngl_full_stage_parameters")
        };
        if (!("backgroundColor" in stage_params)) {
            stage_params["backgroundColor"] = "white";
        }
        this.stage = new NGL.Stage(undefined);
        this.view.$container = $(this.stage.viewer.container);
        this.view.$el.append(this.view.$container);
        this.stage.setParameters(stage_params);
        this.view.$container = $(this.stage.viewer.container);
        this.view.handleResizable();
        this.view.ngl_view_id = this.view.uuid;
        this.view.touch();
        var width = this.view.model.get("_view_width") || this.view.$el.parent().width() + "px";
        var height = this.view.model.get("_view_height") || "300px";
        this.setSize(width, height);
        this.view.uiManager.createFullscreenBtn();
        this.view.uiManager.createIPlayer();
        this.view.GUIStyleChanged();

        this.view.$container.resizable(
            "option", "maxWidth", this.view.$el.parent().width()
        );
        if (this.view.embedHandler.isEmbeded()) {
            console.log("Embed mode for NGLView");
            this.view.embedHandler.handleEmbed();
        } else {
            this.view.requestUpdateStageParameters();
            const viewsDict = await this.view.model.views;
            const views = await Promise.all(Object.values(viewsDict));
            if (views.length === 1) {
                this.view.serialize_camera_orientation();
            } else {
                this.view.set_camera_orientation(this.view.model.get("_camera_orientation"));
            }
        }
    }

    handleResize() {
        var width = this.view.$el.width();
        var height = this.view.$el.height() + "px";
        if (this.view.stage_widget) {
            width = width - $(this.view.stage_widget.sidebar.dom).width();
        }
        width = width + "px";
        this.setSize(width, height);
    }

    setSize(width: string, height: string) {
        this.stage.viewer.container.style.width = width;
        this.stage.viewer.container.style.height = height;
        this.stage.handleResize();
    }
}
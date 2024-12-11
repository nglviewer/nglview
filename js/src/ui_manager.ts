import { StageWidget } from "./gui";
import { StageManager } from "./stage_manager";
import { NGLView } from "./widget_ngl";
import * as NGL from 'ngl';

export class UIManager {
    view: NGLView;

    constructor(view: NGLView) {
        this.view = view;
    }

    async createFullscreenBtn() {
        this.view.btn_pview_fullscreen = this.view.createView("_ibtn_fullscreen");
        var view = await this.view.btn_pview_fullscreen;
        var stage: NGL.Stage = this.view.stage;

        var pe = view.el;
        pe.style.position = 'absolute';
        pe.style.zIndex = 100;
        pe.style.top = '5%';
        pe.style.right = '5%';
        pe.style.opacity = '0.7';
        pe.style.width = '35px';
        pe.style.background = 'white';
        pe.style.opacity = '0.3';
        pe.style.display = 'none';
        pe.onclick = function () {
            this.view.stage.toggleFullscreen();
        }.bind(this);
        stage.viewer.container.append(view.el);
        stage.signals.fullscreenChanged.add((isFullscreen) => {
            if (isFullscreen) {
                view.model.set("icon", "compress");
            } else {
                view.model.set("icon", "expand");
            }
        });
    }

    async createIPlayer() {
        this.view.player_pview = this.view.createView("_iplayer");
        var view = await this.view.player_pview;
        var pe = view.el;
        pe.style.position = 'absolute';
        pe.style.zIndex = 100;
        pe.style.bottom = '5%';
        pe.style.left = '10%';
        pe.style.opacity = '0.7';
        this.view.stage.viewer.container.append(view.el);
        pe.style.display = 'none';
    }

    async createImageBtn() {
        this.view.image_btn_pview = this.view.createView("_ibtn_image");
        var view = await this.view.image_btn_pview;
        var pe = view.el;
        pe.style.position = 'absolute';
        pe.style.zIndex = 100;
        pe.style.top = '5%';
        pe.style.right = '10%';
        pe.style.opacity = '0.7';
        pe.style.width = '35px';
        this.view.stage.viewer.container.append(view.el);
    }

    async createGUI() {
        this.view.pgui_view = this.view.createView("_igui");
        var view = await this.view.pgui_view;
        var pe = view.el;
        pe.style.position = 'absolute';
        pe.style.zIndex = 100;
        pe.style.top = '5%';
        pe.style.right = '10%';
        pe.style.width = '300px';
        this.view.stage.viewer.container.append(view.el);
    }

    createNglGUI() {
        this.view.stage_widget = new StageWidget(this.view);
    }
}
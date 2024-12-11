import { NGLModel, NGLView } from "./widget_ngl";
import * as $ from 'jquery'

export class EventHandler {
    view: NGLView;

    constructor(view: NGLView) {
        this.view = view;
    }

    handleSignals = () => {
        var container = this.view.stage.viewer.container;
        container.addEventListener('mouseover', (e) => {
            this.view._ngl_focused = 1;
            e; // linter
            this.view.mouseOverDisplay('block');
        }, false);

        container.addEventListener('mouseout', (e) => {
            this.view._ngl_focused = 0;
            e; // linter
            this.view.mouseOverDisplay('none');
        }, false);

        container.addEventListener('contextmenu', (e) => {
            e.stopPropagation();
            e.preventDefault();
        }, true);

        this.view.stage.signals.componentAdded.add((component) => {
            this.view.comp_uuids.push(component.uuid);
            var len = this.view.stage.compList.length;
            this.view.model.set("n_components", len);
            this.view.touch();
            var comp = this.view.stage.compList[len - 1];
            comp.signals.representationRemoved.add(() => {
                this.view.request_repr_dict();
            });
            comp.signals.representationAdded.add((repr) => {
                this.view.request_repr_dict();
                repr.signals.parametersChanged.add(() => {
                    console.log("repr.parametersChanged");
                    this.view.request_repr_dict();
                });
            });
        }, this.view);

        this.view.stage.signals.componentRemoved.add(async (component) => {
            var cindex = this.view.comp_uuids.indexOf(component.uuid);
            this.view.comp_uuids.splice(cindex, 1);
            var n_components = this.view.stage.compList.length;
            this.view.model.set("n_components", n_components);
            this.view.touch();
            console.log('componentRemoved', component, component.uuid);

            var pviews = [];
            for (var k in this.view.model.views) {
                pviews.push(this.view.model.views[k]);
            }

            var views = await Promise.all(pviews);
            console.log(views);
            var update_backend = false;
            for (var k in views) {
                var view = views[k];
                if ((view.uuid != this.view.uuid) && (view.stage.compList.length > n_components)) {
                    view.stage.removeComponent(view.stage.compList[cindex]);
                    update_backend = true;
                }
            }
            if (update_backend) {
                console.log("should update backend");
                this.view.send({ "type": "removeComponent", "data": cindex });
            }
        }, this.view);

        this.view.stage.signals.parametersChanged.add(() => {
            this.view.requestUpdateStageParameters();
        }, this.view);

        this.view.stage.viewerControls.signals.changed.add(() => {
            setTimeout(() => {
                this.view.serialize_camera_orientation();
            }, 100);

            var m = this.view.stage.viewerControls.getOrientation();
            if (this.view._synced_model_ids.length > 0 && this.view._ngl_focused == 1) {
                this.view._synced_model_ids.forEach(async (mid) => {
                    var model = await this.view.model.widget_manager.get_model(mid) as NGLModel;
                    for (var k in model.views) {
                        var view = await model.views[k] as NGLView;
                        if (view.uuid != this.view.uuid) {
                            view.stage.viewerControls.orient(m);
                        }
                    }
                });
            }
        });
    }

    handleMessage = () => {
        this.view.model.on("msg:custom", (msg) => {
            this.view.on_msg(msg);
        }, this.view);

        if (this.view.model.comm) {
            this.view.model.comm.on_msg((msg) => {
                var buffers = msg.buffers;
                var content = msg.content.data.content;
                if (buffers.length && content) {
                    content.buffers = buffers;
                }
                this.view.model._handle_comm_msg.call(this.view.model, msg);
            });
        }
    }

    handlePicking = () => {
        this.view.$pickingInfo = $("<div></div>")
            .css("position", "absolute")
            .css("top", "5%")
            .css("left", "3%")
            .css("background-color", "white")
            .css("padding", "2px 5px 2px 5px")
            .css("opacity", "0.7")
            .appendTo(this.view.$container);

        this.view.stage.signals.clicked.add((pd) => {
            if (pd) {
                this.view.model.set('picked', {}); //refresh signal
                this.view.touch();

                var pd2 = {} as any;
                var pickingText = "";
                if (pd.atom) {
                    pd2.atom1 = pd.atom.toObject();
                    pd2.atom1.name = pd.atom.qualifiedName();
                    pickingText = "Atom: " + pd2.atom1.name;
                } else if (pd.bond) {
                    pd2.bond = pd.bond.toObject();
                    pd2.atom1 = pd.bond.atom1.toObject();
                    pd2.atom1.name = pd.bond.atom1.qualifiedName();
                    pd2.atom2 = pd.bond.atom2.toObject();
                    pd2.atom2.name = pd.bond.atom2.qualifiedName();
                    pickingText = "Bond: " + pd2.atom1.name + " - " + pd2.atom2.name;
                }
                if (pd.instance) pd2.instance = pd.instance;

                var n_components = this.view.stage.compList.length;
                for (var i = 0; i < n_components; i++) {
                    var comp = this.view.stage.compList[i];
                    if (comp.uuid == pd.component.uuid) {
                        pd2.component = i;
                    }
                }

                this.view.model.set('picked', pd2);
                this.view.touch();

                this.view.$pickingInfo.text(pickingText);
            }
        });
    }

    async mouseOverDisplay(type) {
        if (this.view.btn_pview_fullscreen) {
            var btn = await this.view.btn_pview_fullscreen
            btn.el.style.display = type
            if (this.view.stage_widget) {
                // If NGL's GUI exists, use its fullscreen button.
                btn.el.style.display = 'none'
            }
        }

        if (this.view.player_pview) {
            var v = await this.view.player_pview
            v.el.style.display = type
            // Need to check if max_frame is available (otherwise NaN)
            // https://github.com/jupyter-widgets/ipywidgets/issues/2485
            if (!this.view.model.get("max_frame") || (this.view.model.get("max_frame") == 0)) {
                // always hide if there's no trajectory.
                v.el.style.display = 'none'
            }
        }
    }

    parametersChanged() {
        var _parameters = this.view.model.get("_parameters");
        this.view.setParameters(_parameters);
    }

    GUIStyleChanged() {
        var style = this.view.model.get("gui_style");
        if (style === 'ngl') {
            this.view.createNglGUI();
        } else {
            if (this.view.stage_widget) {
                this.view.stage_widget.dispose()
                this.view.stage_widget = undefined
                this.view.$container.resizable("enable")
                var width = this.view.$el.parent().width() + "px";
                var height = this.view.$el.parent().height() + "px";
                this.view.setSize(width, height);
            }
        }
    }
}
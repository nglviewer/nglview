import { NGLView } from "./widget_ngl";

export class EventHandler {
    view: NGLView;

    constructor(view: NGLView) {
        this.view = view;
    }

    handleSignals = () => {
        var container = this.view.stage.viewer.container;
        var that = this.view;
        container.addEventListener('mouseover', (e) => {
            that._ngl_focused = 1;
            e; // linter
            that.mouseOverDisplay('block');
        }, false);

        container.addEventListener('mouseout', (e) => {
            that._ngl_focused = 0;
            e; // linter
            that.mouseOverDisplay('none');
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
                that.request_repr_dict();
            });
            comp.signals.representationAdded.add((repr) => {
                that.request_repr_dict();
                repr.signals.parametersChanged.add(() => {
                    console.log("repr.parametersChanged");
                    that.request_repr_dict();
                });
            });
        }, this.view);

        this.view.stage.signals.componentRemoved.add(async (component) => {
            var that = this.view;
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
                if ((view.uuid != that.uuid) && (view.stage.compList.length > n_components)) {
                    view.stage.removeComponent(view.stage.compList[cindex]);
                    update_backend = true;
                }
            }
            if (update_backend) {
                console.log("should update backend");
                that.send({ "type": "removeComponent", "data": cindex });
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
            if (that._synced_model_ids.length > 0 && that._ngl_focused == 1) {
                that._synced_model_ids.forEach(async (mid) => {
                    var model = await that.model.widget_manager.get_model(mid);
                    for (var k in model.views) {
                        var view = await model.views[k];
                        if (view.uuid != that.uuid) {
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
}
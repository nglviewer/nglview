import { NGLView } from "./widget_ngl";

export class RepresentationHandler {
    view: NGLView;

    constructor(view: NGLView) {
        this.view = view;
    }

    requestReprParameters(component_index, repr_index) {
        var comp = this.view.stage.compList[component_index];
        var repr = comp.reprList[repr_index];
        var msg = repr.repr.getParameters();

        if (msg) {
            msg['name'] = repr.name;
            this.view.send({
                'type': 'repr_parameters',
                'data': msg
            });
        }
    }

    request_repr_dict() {
        var repr_dict = this.getReprDictFrontEnd();
        this.view.send({
            'type': 'request_repr_dict',
            'data': repr_dict,
        });
        if (this.view._synced_repr_model_ids.length > 0) {
            this.view._synced_repr_model_ids.forEach(async (mid) => {
                var model = await this.view.model.widget_manager.get_model(mid);
                for (var k in model.views) {
                    var view = await model.views[k];
                    if (view.uuid != this.view.uuid) {
                        view._set_representation_from_repr_dict(repr_dict);
                    }
                }
            });
        }
    }

    getReprDictFrontEnd() {
        var repr_dict = {};
        var n_components = this.view.stage.compList.length;
        for (var i = 0; i < n_components; i++) {
            var comp = this.view.stage.compList[i];
            repr_dict[i] = {};
            var msgi = repr_dict[i];
            for (var j = 0; j < comp.reprList.length; j++) {
                var repr = comp.reprList[j];
                msgi[j] = {};
                msgi[j]['type'] = repr.name;
                msgi[j]['params'] = repr.repr.getParameters();
            }
        }
        return repr_dict;
    }

    syncReprForAllViews() {
        var repr_dict_backend = this.view.model.get("_ngl_repr_dict");
        var repr_dict_frontend = this.getReprDictFrontEnd();
        if (JSON.stringify(repr_dict_frontend) !== JSON.stringify(repr_dict_backend)) {
            this._set_representation_from_repr_dict(repr_dict_backend);
        }
    }

    async syncReprWithMe() {
        var repr_dict = this.getReprDictFrontEnd();
        for (var k in this.view.model.views) {
            var v = await this.view.model.views[k];
            if (v.uuid != this.view.uuid) {
                v._set_representation_from_repr_dict(repr_dict);
            }
        }
        this.request_repr_dict();
    }

    setSyncRepr(model_ids) {
        this.view._synced_repr_model_ids = model_ids;
    }

    set_representation_from_backend() {
        var repr_dict = this.view.model.get('_ngl_repr_dict');
        this._set_representation_from_repr_dict(repr_dict);
    }

    _set_representation_from_repr_dict(repr_dict) {
        var compList = this.view.stage.compList;
        if (compList.length > 0) {
            for (var index in repr_dict) {
                var comp = compList[index];
                comp.removeAllRepresentations();
                var reprlist = repr_dict[index];
                for (var j in reprlist) {
                    var repr = reprlist[j];
                    if (repr) {
                        comp.addRepresentation(repr.type, repr.params);
                    }
                }
            }
        }
    }

    setVisibilityForRepr(component_index, repr_index, value) {
        var comp = this.view.stage.compList[component_index];
        var repr = comp.reprList[repr_index];
        repr.setVisibility(value);
    }

    removeRepresentation(component_index, repr_index) {
        var comp = this.view.stage.compList[component_index];
        comp.removeRepresentation(comp.reprList[repr_index]);
    }

    removeRepresentationsByName(repr_name, component_index) {
        var comp = this.view.stage.compList[component_index];
        comp.reprList = comp.reprList.filter(repr => repr.name !== repr_name);
    }

    updateRepresentationForComponent(repr_index, component_index, params) {
        var comp = this.view.stage.compList[component_index];
        var repr = comp.reprList[repr_index];
        repr.setParameters(params);
    }

    updateRepresentationsByName(repr_name, component_index, params) {
        var comp = this.view.stage.compList[component_index];
        comp.reprList.forEach(repr => {
            if (repr.name === repr_name) {
                repr.setParameters(params);
            }
        });
    }
}
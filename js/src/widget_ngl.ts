import * as widgets from "@jupyter-widgets/base"
import * as NGL from "ngl"
import * as $ from 'jquery'
import * as _ from 'underscore'
import "./lib/signals.min.js"
import "./lib/colorpicker.min.js"
import "./ui/ui.js"
import "./ui/ui.extra.js"
import "./ui/ui.ngl.js"
import 'jquery-ui/ui/widgets/dialog' // FIXME: remove?
import 'jquery-ui/themes/base/all.css' // FIXME: remove?

import { StageWidget } from "./gui"
import { FullscreenModel, FullscreenView } from "./fullscreen"
import { ColormakerRegistryModel, ColormakerRegistryView } from "./color"
import { ThemeManagerModel, ThemeManagerView } from "./theme"
import { MolstarModel, MolstarView } from "./molstarview/widget"
import { StageManager } from "./stage_manager"
import { EventHandler } from "./event_handler"
import { UIManager } from "./ui_manager"
import { EmbedHandler } from "./embed_handler"
import { RepresentationHandler } from "./representation_handler"
import { MessageHandler } from "./message_handler"


// From NGL
// http://www.broofa.com/Tools/Math.uuid.htm
const chars = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'.split('')
const uuid = new Array(36)

function generateUUID() {
    let rnd = 0
    let r

    for (let i = 0; i < 36; i++) {
        if (i === 8 || i === 13 || i === 18) {
            uuid[i] = '-'
        } else if (i === 14) {
            uuid[i] = '4'
        } else {
            if (rnd <= 0x02) rnd = 0x2000000 + (Math.random() * 0x1000000) | 0
            r = rnd & 0xf
            rnd = rnd >> 4
            uuid[i] = chars[(i === 19) ? (r & 0x3) | 0x8 : r]
        }
    }

    return uuid.join('')
}


export class NGLModel extends widgets.DOMWidgetModel {
    defaults() {
        return _.extend(widgets.DOMWidgetModel.prototype.defaults(), {
            _model_name: 'NGLModel',
            _model_module: 'nglview-js-widgets',
            _model_module_version: require("../package.json").version,
            _view_name: "NGLView",
            _view_module: "nglview-js-widgets",
            _view_module_version: require("../package.json").version,
        });
    }
}

export class NGLView extends widgets.DOMWidgetView {
    stage: NGL.Stage;
    stageManager: StageManager;
    eventHandler: EventHandler;
    uiManager: UIManager;
    embedHandler: EmbedHandler;
    representationHandler: RepresentationHandler;
    messageHandler: MessageHandler;
    model: NGLModel;
    player_pview: Promise<any>;
    btn_pview_fullscreen: Promise<any>;
    image_btn_pview: Promise<any>;
    pgui_view: Promise<any>;
    $container: any;
    ngl_view_id: string;
    uuid: string;
    stage_widget: StageWidget;
    _ngl_focused: number;
    comp_uuids: string[];
    _synced_model_ids: string[];
    _synced_repr_model_ids: string[];
    atomColor: any;
    $pickingInfo: any;


    constructor(options: any) {
        super(options);
        this.stageManager = new StageManager(this);
        this.eventHandler = new EventHandler(this);
        this.uiManager = new UIManager(this);
        this.embedHandler = new EmbedHandler(this);
        this.representationHandler = new RepresentationHandler(this);
        this.messageHandler = new MessageHandler(this);
    }

    render() {
        this.beforeDisplay();
        this.displayed.then(() => {
            this.stageManager.createStage();
            this.stage = this.stageManager.stage;
            this.eventHandler.handlePicking();
            this.eventHandler.handleSignals();
            this.eventHandler.handleMessage();
            this.finalizeDisplay();
        });
    }

    beforeDisplay() {
        this.model.on("change:_parameters", this.parametersChanged, this);
        this.model.on("change:gui_style", this.GUIStyleChanged, this);
        this.model.set('_ngl_version', NGL.Version);
        this._ngl_focused = 0;
        this.uuid = generateUUID();
        this.stage_widget = undefined;
        this.comp_uuids = [];
        this._synced_model_ids = this.model.get("_synced_model_ids");
        this._synced_repr_model_ids = this.model.get("_synced_repr_model_ids");

        if (this.embedHandler.isEmbeded()) {
            this.embedHandler._handleEmbedBeforeStage();
        }
    }

    finalizeDisplay() {
        // for callbacks from Python
        // must be after initializing NGL.Stage
        this.send({
            'type': 'request_loaded',
            'data': true
        })
        var state_params = this.stage.getParameters();
        this.model.set('_ngl_original_stage_parameters', state_params);
        this.touch();
        if (!this.embedHandler.isEmbeded() && this.stage.compList.length < this.model.get("n_components")) {
            // only call this in notebook to avoid calling handleEmbed twice in embeded mode.
            this.embedHandler.handleEmbed()
        }
        var ngl_view_ids = this.model.get("_ngl_view_id")
        ngl_view_ids.push(this.ngl_view_id)
        this.send({ "type": "updateIDs", "data": ngl_view_ids })

        // FIXME: Why below doesn't update _ngl_view_id in backend?
        // this.model.set("_ngl_view_id", ngl_view_ids)
        // this.touch()
    }

    async mouseOverDisplay(type) {
        if (this.btn_pview_fullscreen) {
            var btn = await this.btn_pview_fullscreen;
            btn.el.style.display = type;
            if (this.stage_widget) {
                // If NGL's GUI exists, use its fullscreen button.
                btn.el.style.display = 'none';
            }
        }

        if (this.player_pview) {
            var v = await this.player_pview;
            v.el.style.display = type;
            // Need to check if max_frame is available (otherwise NaN)
            // https://github.com/jupyter-widgets/ipywidgets/issues/2485
            if (!this.model.get("max_frame") || (this.model.get("max_frame") == 0)) {
                // always hide if there's no trajectory.
                v.el.style.display = 'none';
            }
        }
    }

    updateNGLTheme(css_content) {
        var ele = document.getElementById("nglview_style")
        if (ele != undefined) {
            document.head.removeChild(ele)
        }
        var style = document.createElement("style")
        style.id = 'nglview_style'
        style.type = 'text/css'
        style.appendChild(document.createTextNode(css_content))
        document.head.appendChild(style)
    }

    serialize_camera_orientation() {
        var m = this.stage.viewerControls.getOrientation();
        this.model.set('_camera_orientation', m.elements);
        this.touch();
    }

    set_camera_orientation(orientation) {
        if (orientation.length > 0) {
            this.stage.viewerControls.orient(orientation);
            this.serialize_camera_orientation();
        }
    }

    executeCode(code) {
        eval(code);
    }

    updateCoordinatesFromDict(cdict, frame_index) {
        // update coordinates for given "index"
        // cdict = Dict[int, List[base64]]
        var keys = Object.keys(cdict).filter(k => (k !== 'n_frames'));

        for (var i = 0; i < keys.length; i++) {
            var traj_index = parseInt(keys[i], 10);
            var coordinates = this.decode_base64(cdict[traj_index][frame_index]);
            if (coordinates && coordinates.byteLength > 0) {
                this.updateCoordinates(coordinates, traj_index);
            }
        }
    }

    requestFrame() {
        this.send({
            'type': 'request_frame',
            'data': 'frame'
        });
    }

    requestUpdateStageParameters() {
        var updated_params = this.stage.getParameters();
        this.model.set('_ngl_full_stage_parameters', updated_params);
        this.touch();
    }

    requestReprParameters(component_index, repr_index) {
        this.representationHandler.requestReprParameters(component_index, repr_index);
    }

    request_repr_dict() {
        this.representationHandler.request_repr_dict();
    }

    getReprDictFrontEnd() {
        return this.representationHandler.getReprDictFrontEnd();
    }

    syncReprForAllViews() {
        this.representationHandler.syncReprForAllViews();
    }

    async syncReprWithMe() {
        await this.representationHandler.syncReprWithMe();
    }

    setSyncRepr(model_ids) {
        this.representationHandler.setSyncRepr(model_ids);
    }

    setSyncCamera(model_ids) {
        this._synced_model_ids = model_ids;
    }

    viewXZPlane() {
        var m = new NGL.Matrix4().makeRotationX(Math.PI / 2);
        var q = new NGL.Quaternion().setFromRotationMatrix(m);
        this.stage.viewerControls.rotate(q);
    }

    set_representation_from_backend() {
        this.representationHandler.set_representation_from_backend();
    }

    _set_representation_from_repr_dict(repr_dict) {
        this.representationHandler._set_representation_from_repr_dict(repr_dict);
    }

    async createView(trait_name) {
        // Create a view for the model with given `trait_name`
        // e.g: in backend, 'view.<trait_name>`
        console.log("Creating view for model " + trait_name);
        var manager = this.model.widget_manager;
        var model_id = this.model.get(trait_name).replace("IPY_MODEL_", "");
        return await manager.create_view(await manager.get_model(model_id))
    }

    getPlayerModel() {
        // return a Promise
        var model_id = this.model.get("_iplayer").replace("IPY_MODEL_", "");
        return this.model.widget_manager.get_model(model_id)
    }

    async createIPlayer() {
        this.uiManager.createIPlayer();
    }

    async createImageBtn() {
        this.uiManager.createImageBtn();
    }

    async createFullscreenBtn() {
        this.uiManager.createFullscreenBtn();
    }

    async createGUI() {
        this.uiManager.createGUI();
    }

    createNglGUI() {
        this.stage_widget = new StageWidget(this)
    }


    setVisibilityForRepr(component_index, repr_index, value) {
        this.representationHandler.setVisibilityForRepr(component_index, repr_index, value);
    }

    removeRepresentation(component_index, repr_index) {
        this.representationHandler.removeRepresentation(component_index, repr_index);
    }

    removeRepresentationsByName(repr_name, component_index) {
        this.representationHandler.removeRepresentationsByName(repr_name, component_index);
    }

    updateRepresentationForComponent(repr_index, component_index, params) {
        this.representationHandler.updateRepresentationForComponent(repr_index, component_index, params);
    }

    updateRepresentationsByName(repr_name, component_index, params) {
        this.representationHandler.updateRepresentationsByName(repr_name, component_index, params);
    }

    setColorByResidue = (colors, component_index, repr_index) => {
        var repr = this.stage.compList[component_index].reprList[repr_index];
        var schemeId = NGL.ColormakerRegistry.addScheme((params) => {
            this.atomColor = (atom) => {
                var color = colors[atom.residueIndex];
                return color;
            };
            params; // to pass eslint; ack;
        });
        repr.setColor(schemeId);
    }

    addShape(name, shapes) {
        // shapes: List[Tuple[str, ...]]
        // e.g: [('sphere', ...), ('cone', ...)]
        var shape = new NGL.Shape(name);
        var shape_dict = {
            'sphere': shape.addSphere,
            'ellipsoid': shape.addEllipsoid,
            'cylinder': shape.addCylinder,
            'cone': shape.addCone,
            'mesh': shape.addMesh,
            'arrow': shape.addArrow,
            'text': shape.addText,
            'label': shape.addText,
            'tetrahedron': shape.addTetrahedron,
            'octahedron': shape.addOctahedron,
            'torus': shape.addTorus
        };
        for (var i = 0; i < shapes.length; i++) {
            var shapes_i = shapes[i]
            var shape_type = shapes_i[0];
            var params = shapes_i.slice(1, shapes_i.length);
            // e.g params = ('sphere', [ 0, 0, 9 ], [ 1, 0, 0 ], 1.5)

            var func = shape_dict[shape_type];
            func.apply(shape, params);
            // shape.func(params);
        }
        var shapeComp = this.stage.addComponentFromObject(shape);
        if (shapeComp) {
            shapeComp.addRepresentation("buffer", {});
        }
    }

    addBuffer(name, kwargs) {
        var class_dict = {
            "arrow": NGL.ArrowBuffer,
            "box": NGL.BoxBuffer,
            "cone": NGL.ConeBuffer,
            "cylinder": NGL.CylinderBuffer,
            "ellipsoid": NGL.EllipsoidBuffer,
            "octahedron": NGL.OctahedronBuffer,
            "sphere": NGL.SphereBuffer,
            "text": NGL.TextBuffer,
            "tetrahedron": NGL.TetrahedronBuffer,
            "torus": NGL.TorusBuffer
        };

        var params = {};
        for (var key in kwargs) {
            params[key] = new Float32Array(kwargs[key]);
        }
        var shape = new NGL.Shape("shape");
        var buffer_class = class_dict[name];
        var buffer = new buffer_class(params);
        shape.addBuffer(buffer);
        var shapeComp = this.stage.addComponentFromObject(shape);
        if (shapeComp) {
            shapeComp.addRepresentation("buffer", {});
        }
    }

    async replaceStructure(structure) {
        var blob = new Blob([structure.data], { type: "text/plain" });
        var stage = this.stage;
        var params = structure.params || {};
        params.ext = structure.ext;
        params.defaultRepresentation = false;
        var comp = this.stage.compList[0];
        var representations = comp.reprList.slice();
        var old_orientation = this.stage.viewerControls.getOrientation();
        var component = await this.stage.loadFile(blob, params);
        stage.viewerControls.orient(old_orientation);
        representations.forEach((repr) => {
            var repr_name = repr.name;
            var repr_params = repr.repr.getParameters();
            // Note: not using repr.repr.type, repr.repr.params
            // since seems to me that repr.repr.params won't return correct "sele"
            if (component) {
                component.addRepresentation(repr_name, repr_params);
            }
        });
        stage.removeComponent(comp);
        this._handleLoadFileFinished();
    }

    superpose(cindex0, cindex1, align, sele0, sele1) {
        // superpose two components with given params
        var component0 = this.stage.compList[cindex0];
        var component1 = this.stage.compList[cindex1];
        (component1 as NGL.StructureComponent).superpose(component0 as NGL.StructureComponent, align, sele0, sele1);
    }

    decode_base64(base64) {
        // lightly adapted from Niklas

        /*
         * base64-arraybuffer
         * https://github.com/niklasvh/base64-arraybuffer
         *
         * Copyright (c) 2012 Niklas von Hertzen
         * Licensed under the MIT license.
         */
        var chars =
            "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
        var bufferLength = base64.length * 0.75,
            len = base64.length,
            i, p = 0,
            encoded1, encoded2, encoded3, encoded4;

        if (base64[base64.length - 1] === "=") {
            bufferLength--;
            if (base64[base64.length - 2] === "=") {
                bufferLength--;
            }
        }

        var arraybuffer = new ArrayBuffer(bufferLength),
            bytes = new Uint8Array(arraybuffer);

        for (i = 0; i < len; i += 4) {
            encoded1 = chars.indexOf(base64[i]);
            encoded2 = chars.indexOf(base64[i + 1]);
            encoded3 = chars.indexOf(base64[i + 2]);
            encoded4 = chars.indexOf(base64[i + 3]);

            bytes[p++] = (encoded1 << 2) | (encoded2 >> 4);
            bytes[p++] = ((encoded2 & 15) << 4) | (encoded3 >> 2);
            bytes[p++] = ((encoded3 & 3) << 6) | (encoded4 & 63);
        }

        return arraybuffer;
    }

    updateCoordinates(coordinates: ArrayBuffer, model: number) {
        // coordinates must be ArrayBuffer (use this.decode_base64)
        var component = this.stage.compList[model];
        if (coordinates && component) {
            var coords = new Float32Array(coordinates);
            (component as NGL.StructureComponent).structure.updatePosition(coords);
            component.updateRepresentations({
                "position": true
            });
        }
    }

    handleResizable() {
        this.$container.resizable({
            resize: function (event, ui) {
                this.setSize(ui.size.width + "px", ui.size.height + "px");
            }.bind(this)
        })
    }

    handleResize() {
        var width = this.$el.width()
        var height = this.$el.height() + "px"
        if (this.stage_widget) {
            width = width - $(this.stage_widget.sidebar.dom).width()
        }
        width = width + "px"
        this.setSize(width, height)
    }

    setSize(width, height) {
        this.stage.viewer.container.style.width = width;
        this.stage.viewer.container.style.height = height;
        this.stage.handleResize();
    }

    GUIStyleChanged() {
        var style = this.model.get("gui_style");
        if (style === 'ngl') {
            this.createNglGUI();
        } else {
            if (this.stage_widget) {
                this.stage_widget.dispose()
                this.stage_widget = undefined
                this.$container.resizable("enable")
                var width = this.$el.parent().width() + "px";
                var height = this.$el.parent().height() + "px";
                this.setSize(width, height);
            }
        }
    }

    parametersChanged() {
        var _parameters = this.model.get("_parameters");
        this.setParameters(_parameters);
    }

    setParameters(parameters) {
        this.stage.setParameters(parameters);

        // do not set _ngl_full_stage_parameters here
        // or parameters will be never updated (not sure why)
        // use observe in python side
        var updated_params = this.stage.getParameters();
        this.send({
            'type': 'stage_parameters',
            'data': updated_params
        })
    }

    async _downloadImage(filename, params) {
        if (this.ngl_view_id == this.get_last_child_id()) {
            var blob = await this.stage.makeImage(params)
            NGL.download(blob, filename);
        }
    }

    async _exportImage(wid, params) {
        if (this.ngl_view_id == this.get_last_child_id()) {
            var blob = await this.stage.makeImage(params)
            var reader = new FileReader();
            var arr_str;
            reader.onload = function () {
                arr_str = (reader.result as string).replace("data:image/png;base64,", "");
                // this.model.set("_image_data", arr_str);
                // this.touch();
                this.send({
                    "data": arr_str,
                    "type": "image_data",
                    "ID": wid,
                });
                this.send({ 'type': 'async_message', 'data': 'ok' });
            }.bind(this);
            reader.readAsDataURL(blob);
        }
    }

    async handleMovieMaking(render_params) {
        if (this.ngl_view_id == this.get_last_child_id()) {
            var blob = await this.stage.makeImage(render_params)
            var reader = new FileReader();
            var arr_str;
            reader.onload = function () {
                arr_str = (reader.result as string).replace("data:image/png;base64,", "");
                this.send({
                    "data": arr_str,
                    "type": "movie_image_data",
                }); // tell backend that image render is finished,
                // backend will send next frame's coordinates.
                this.send({ 'type': 'async_message', 'data': 'ok' });
            }.bind(this);
            reader.readAsDataURL(blob);
        }
    }

    async _handleLoadFileFinished() {
        this.send({ 'type': 'async_message', 'data': 'ok' });
    }

    _getLoadFilePromise(msg) {
        // args = [{'type': ..., 'data': ...}]
        var args0 = msg.args[0];
        if (args0.type == 'blob') {
            var blob;
            if (args0.binary) {
                var decoded_data = this.decode_base64(args0.data);
                blob = new Blob([decoded_data], {
                    type: "application/octet-binary"
                });
            } else {
                blob = new Blob([args0.data], {
                    type: "text/plain"
                });
            }
            return this.stage.loadFile(blob, msg.kwargs)
        } else {
            var file = new File([""], args0.data);
            // FIXME: if not "any", typescipt complains there is no
            // "exists" method.
            var path = "";
            // FIXME: uncomment ("exists" does not exist)
            // if (args0.type == 'path' && ! file.exists){
            // hacky fix for jupyterlab
            // https://github.com/arose/nglview/issues/783
            // https://github.com/jupyterlab/jupyterlab/issues/6218
            // e.g: http://localhost:8888/files/ala3.pdb
            //    path = 'files/' + args0.data;
            // } else {
            // path = args0.data;
            //}
            path = args0.data;
            return this.stage.loadFile(path, msg.kwargs)
        }
    }

    get_last_child_id() {
        var keys = this.model.get('_ngl_view_id')
        return keys[keys.length - 1]
    }

    async _handleStageLoadFile(msg) {
        // args = [{'type': ..., 'data': ...}]
        if (this.ngl_view_id != this.get_last_child_id() && msg.last_child) {
            return
        }
        var o = await this._getLoadFilePromise(msg)
        this._handleLoadFileFinished();
    }

    addColorScheme(args, label) {
        var id = NGL.ColormakerRegistry.addSelectionScheme(args, label);
        var scheme = NGL.ColormakerRegistry.userSchemes[id];
        NGL.ColormakerRegistry.removeScheme(id);
        // hard code the scheme ID
        NGL.ColormakerRegistry.add(label, scheme);
        return label
    }

    async on_msg(msg) {
        this.messageHandler.on_msg(msg);
    }

    async handleCallMethod(msg) {
        this.messageHandler.handleCallMethod(msg);
    }

    async handleStageMethod(msg, new_args) {
        this.messageHandler.handleStageMethod(msg, new_args);
    }

    handleBase64Single(coordinatesDict) {
        this.messageHandler.handleBase64Single(coordinatesDict);
    }

    handleBinarySingle(msg) {
        this.messageHandler.handleBinarySingle(msg);
    }

    handleGetRequest(data) {
        this.messageHandler.handleGetRequest(data);
    }
}

// export all models and views here to make embeding a bit easier
module.exports = {
    'NGLView': NGLView,
    'NGLModel': NGLModel,
    'NGL': NGL,
    'FullscreenModel': FullscreenModel,
    'FullscreenView': FullscreenView,
    'ColormakerRegistryModel': ColormakerRegistryModel,
    'ColormakerRegistryView': ColormakerRegistryView,
    'ThemeManagerModel': ThemeManagerModel,
    'ThemeManagerView': ThemeManagerView,
    'MolstarModel': MolstarModel,
    'MolstarView': MolstarView,
}

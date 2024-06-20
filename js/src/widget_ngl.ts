var Jupyter
var widgets = require("@jupyter-widgets/base")
var NGL = require('ngl')
var BaseView = require('./base').BaseView
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
import { ThemeManagerModel, ThemeManagerView} from "./theme"

NGL.nglview_debug = false

// From NGL
// http://www.broofa.com/Tools/Math.uuid.htm
const chars = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'.split('')
const uuid = new Array(36)

function generateUUID () {
  let rnd = 0
  let r

  for (let i = 0; i < 36; i++) {
    if (i === 8 || i === 13 || i === 18 || i === 23) {
      uuid[ i ] = '-'
    } else if (i === 14) {
      uuid[ i ] = '4'
    } else {
      if (rnd <= 0x02) rnd = 0x2000000 + (Math.random() * 0x1000000) | 0
      r = rnd & 0xf
      rnd = rnd >> 4
      uuid[ i ] = chars[ (i === 19) ? (r & 0x3) | 0x8 : r ]
    }
  }

  return uuid.join('')
}

enum Targets {
    Stage = 'Stage',
    Viewer = 'Viewer',
    ViewerControls = 'viewerControls',
    CompList = 'compList',
    Widget = 'Widget',
    Representation = 'Representation'
}

async function createView(that, trait_name){
    // Create a view for the model with given `trait_name`
    // e.g: in backend, 'view.<trait_name>`
    console.log("Creating view for model " + trait_name);
    var manager = that.model.widget_manager
    var model_id = that.model.get(trait_name).replace("IPY_MODEL_", "");
    return await manager.create_view(await manager.get_model(model_id))
}


export
class NGLModel extends widgets.DOMWidgetModel{
    defaults(){
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

export
class NGLView extends widgets.DOMWidgetView{
    render(){
        this.beforeDisplay()
        this.displayed.then(function() {
            // move all below code inside 'displayed'
            // to make sure the NGLView and NGLModel are created
            this.createStage()
            this.handlePicking()
            this.handleSignals()
            this.handleMessage()
            this.finalizeDisplay()
        }.bind(this));

    }

    beforeDisplay(){
        this.model.on("change:_parameters", this.parametersChanged, this);
        this.model.on("change:gui_style", this.GUIStyleChanged, this);
        this.model.set('_ngl_version', NGL.Version);
        this._ngl_focused = 0
        this.uuid = generateUUID()
        this.stage_widget = undefined
        this.comp_uuids = []
        this._synced_model_ids = this.model.get("_synced_model_ids");
        this._synced_repr_model_ids = this.model.get("_synced_repr_model_ids")

        if (this.isEmbeded()){
            // embed mode
            this._handleEmbedBeforeStage()
        }
    }

    createStage(){
        // init NGL stage
        var stage_params = {
            // Shallow copy so that _ngl_full_stage_parameters is not updated yet
            ...this.model.get("_ngl_full_stage_parameters")
        };
        if (!("backgroundColor" in stage_params)){
            stage_params["backgroundColor"] = "white"
        }
        NGL.useWorker = false;
        var view_parent = this.options.parent
        this.stage = new NGL.Stage(undefined)
        this.$container = $(this.stage.viewer.container);
        this.$el.append(this.$container)
        this.stage.setParameters(stage_params);
        this.$container = $(this.stage.viewer.container);
        this.handleResizable()
        this.ngl_view_id = this.uuid
        this.touch();
        var that = this;
        var width = this.model.get("_view_width") || this.$el.parent().width() + "px";
        var height = this.model.get("_view_height") || "300px";
        this.setSize(width, height);
        this.createFullscreenBtn(); // FIXME: move up?
        this.createIPlayer(); // FIXME: move up?
        this.GUIStyleChanged(); // must be called after displaying to get correct width and height

        this.$container.resizable(
            "option", "maxWidth", this.$el.parent().width()
        );
        if (this.isEmbeded()){
            console.log("Embed mode for NGLView")
            that.handleEmbed();
        }else{
            this.requestUpdateStageParameters();
            if (this.model.views.length == 1){
                this.serialize_camera_orientation();
            }else{
                this.set_camera_orientation(that.model.get("_camera_orientation"));
            }
        }
    }

    isEmbeded(){
        return (this.model.get("_ngl_serialize") || (this.model.comm == undefined))
    }

    handleMessage(){
        this.model.on("msg:custom", function(msg){
           this.on_msg(msg);
        }, this);

        if (this.model.comm) {
            this.model.comm.on_msg(function(msg) {
                var buffers = msg.buffers;
                var content = msg.content.data.content;
                if (buffers.length && content) {
                    content.buffers = buffers;
                }
                this.model._handle_comm_msg.call(this.model, msg);
            }.bind(this));
        }
    }


    finalizeDisplay(){
      // for callbacks from Python
      // must be after initializing NGL.Stage
      this.send({
          'type': 'request_loaded',
          'data': true
      })
      var state_params = this.stage.getParameters();
      this.model.set('_ngl_original_stage_parameters', state_params);
      this.touch();
      if (!this.isEmbeded() && this.stage.compList.length < this.model.get("n_components")) {
          // only call this in notebook to avoid calling handleEmbed twice in embeded mode.
          this.handleEmbed()
      }
      var ngl_view_ids = this.model.get("_ngl_view_id")
      ngl_view_ids.push(this.ngl_view_id)
      this.send({"type": "updateIDs", "data": ngl_view_ids})

      // FIXME: Why below doesn't update _ngl_view_id in backend?
      // this.model.set("_ngl_view_id", ngl_view_ids)
      // this.touch()
    }

    handleSignals(){
      var container = this.stage.viewer.container;
      var that = this;
      container.addEventListener('mouseover', function(e) {
          that._ngl_focused = 1;
          e; // linter
          that.mouseOverDisplay('block')
      }, false);

      container.addEventListener('mouseout', function(e) {
          that._ngl_focused = 0;
          e; // linter
          that.mouseOverDisplay('none')
      }, false);

      container.addEventListener('contextmenu', function(e) {
          e.stopPropagation();
          e.preventDefault();
      }, true);

      this.stage.signals.componentAdded.add(function(component) {
          this.comp_uuids.push(component.uuid)
          var len = this.stage.compList.length;
          this.model.set("n_components", len);
          this.touch();
          var comp = this.stage.compList[len - 1];
          comp.signals.representationRemoved.add(function() {
              that.request_repr_dict();
          });
          comp.signals.representationAdded.add(function(repr) {
              that.request_repr_dict();
              repr.signals.parametersChanged.add(function(){
                  console.log("repr.parametersChanged")
                  that.request_repr_dict();
              })
          });
      }, this);

      this.stage.signals.componentRemoved.add(async function(component) {
          var that = this
          var cindex = this.comp_uuids.indexOf(component.uuid)
          this.comp_uuids.splice(cindex, 1)
          var n_components = this.stage.compList.length
          this.model.set("n_components", n_components)
          this.touch()
          console.log('componentRemoved', component, component.uuid)

          var pviews = []
          for (var k in this.model.views){
              pviews.push(this.model.views[k])
          }

          var views = await Promise.all(pviews)
          console.log(views)
          var update_backend = false
          for (var k in views){
              var view = views[k]
              if ((view.uuid != that.uuid) && (view.stage.compList.length > n_components)){
                  // remove component from NGL's GUI
                  // pass
                  view.stage.removeComponent(view.stage.compList[cindex])
                  update_backend = true
              }
          }
          if (update_backend){
              console.log("should update backend")
              that.send({"type": "removeComponent", "data": cindex})
          }
      }, this);

      this.stage.signals.parametersChanged.add(function(){
          this.requestUpdateStageParameters();
      }, this);

      this.stage.viewerControls.signals.changed.add(function() {
          this.serialize_camera_orientation();
          var m = this.stage.viewerControls.getOrientation();
          if (that._synced_model_ids.length > 0 && that._ngl_focused == 1){
              that._synced_model_ids.forEach(async function(mid){
                  var model = await that.model.widget_manager.get_model(mid)
                  for (var k in model.views){
                      var pview = model.views[k];
                      var view = await model.views[k]
                      if (view.uuid != that.uuid){
                          view.stage.viewerControls.orient(m);
                      }
                  }
              })
          }
      }.bind(this));

    }

    handlePicking(){
      this.$pickingInfo = $("<div></div>")
          .css("position", "absolute")
          .css("top", "5%")
          .css("left", "3%")
          .css("background-color", "white")
          .css("padding", "2px 5px 2px 5px")
          .css("opacity", "0.7")
          .appendTo(this.$container);

      var that = this;
      this.stage.signals.clicked.add(function (pd) {
          if (pd) {
              this.model.set('picked', {}); //refresh signal
              this.touch();

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

              var n_components = this.stage.compList.length;
              for (var i = 0; i < n_components; i++) {
                  var comp = this.stage.compList[i];
                  if (comp.uuid == pd.component.uuid) {
                      pd2.component = i;
                  }
              }

              this.model.set('picked', pd2);
              this.touch();

              this.$pickingInfo.text(pickingText);
          }
      }, this);
    }

    async mouseOverDisplay(type){
        var that = this;
        if (this.btn_pview_fullscreen){
            var btn = await this.btn_pview_fullscreen
            btn.el.style.display = type
            if (that.stage_widget){
                // If NGL's GUI exists, use its fullscreen button.
                btn.el.style.display = 'none'
            }
        }

        var that = this;
        if (this.player_pview){
            var v = await this.player_pview
                v.el.style.display = type
                // Need to check if max_frame is available (otherwise NaN)
                // https://github.com/jupyter-widgets/ipywidgets/issues/2485
                if (!that.model.get("max_frame") || (that.model.get("max_frame") == 0)){
                    // always hide if there's no trajectory.
                    v.el.style.display = 'none'
            }
        }
    }

    updateNGLTheme(css_content){
        var ele = document.getElementById("nglview_style")
        if (ele != undefined){
            document.head.removeChild(ele)
        }
        var style = document.createElement("style")
        style.id = 'nglview_style'
        style.type = 'text/css'
        style.appendChild(document.createTextNode(css_content))
        document.head.appendChild(style)
    }

    serialize_camera_orientation(){
         var m = this.stage.viewerControls.getOrientation();
         this.model.set('_camera_orientation', m.elements);
         this.touch();
    }

    set_camera_orientation(orientation){
        if (orientation.length > 0){
            this.stage.viewerControls.orient(orientation);
            this.serialize_camera_orientation();
        }
    }

    executeCode(code){
        eval(code);
    }

    _handleEmbedBeforeStage(){
        // Only need to reconstruct colors in embeding mode (outside notebook)
        // FIXME: remove this function
        var that = this
        var ngl_color_dict = that.model.get("_ngl_color_dict");
        var label
            // Old API (_ColorScheme)
        for (label in ngl_color_dict){
            if (!NGL.ColormakerRegistry.hasScheme(label)){
                that.addColorScheme(ngl_color_dict[label], label);
            }
        }
    }

    async handleEmbed() {
        const { _ngl_msg_archive, _ngl_full_stage_parameters, _camera_orientation, _ngl_coordinate_resource, _ngl_repr_dict, comm } = this.model.attributes;

        if (Object.keys(_ngl_full_stage_parameters).length === 0 && _camera_orientation.length === 0) {
            console.log("No state stored; initializing embedded widget for the first time.");
            for (const msg of _ngl_msg_archive) {
                await this.on_msg(msg);
            }
            return;
        }

        const loadFilePromises = _ngl_msg_archive
            .filter(msg => msg.methodName === 'loadFile')
            .map(msg => {
                if (msg.kwargs && msg.kwargs.defaultRepresentation) {
                    // no need to add default representation as all representations
                    // are serialized separately, also it unwantedly sets the orientation
                    msg.kwargs.defaultRepresentation = false;
                }
                return this._getLoadFilePromise(msg);
            });

        const compList = await Promise.all(loadFilePromises);
        this.stage.setParameters(_ngl_full_stage_parameters);
        this.set_camera_orientation(_camera_orientation);
        this.touch();

        // Outside notebook
        if (comm === undefined) {
            const n_frames = _ngl_coordinate_resource['n_frames'] || 1;
            this.model.set("max_frame", n_frames - 1);  // trigger updating slider and player's max
            this.touch();
            const model = await this.getPlayerModel();
            const pmodel = model.get("children")[0];
            this.listenTo(pmodel, "change:value", () => this.updateCoordinatesFromDict(_ngl_coordinate_resource, pmodel.get("value")));
        }

        // fire any msg with "fire_embed"
        for (const msg of _ngl_msg_archive) {
            if (msg.fire_embed) {
                await this.on_msg(msg);
            }
        }

        // Must call _set_representation_from_repr_dict after "fire_embed"
        // User might add Shape (buffer component) to the view and the buffer component
        // is not created yet via loadFile
        // https://github.com/nglviewer/nglview/issues/1003
        this._set_representation_from_repr_dict(_ngl_repr_dict);
        this.handleResize(); // FIXME: really need this?
    }

    updateCoordinatesFromDict(cdict, frame_index){
        // update coordinates for given "index"
        // cdict = Dict[int, List[base64]]
        var keys = Object.keys(cdict).filter(k => (k !== 'n_frames'));

        for (var i = 0; i < keys.length; i++) {
            var traj_index = keys[i];
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

    requestUpdateStageParameters(){
        var updated_params = this.stage.getParameters();
        this.model.set('_ngl_full_stage_parameters', updated_params);
        this.touch();
    }

    requestReprParameters(component_index, repr_index) {
        var comp = this.stage.compList[component_index];
        var repr = comp.reprList[repr_index];
        var msg = repr.repr.getParameters();

        if (msg) {
            msg['name'] = repr.name;
            this.send({
                'type': 'repr_parameters',
                'data': msg
            });
        }
    }

    request_repr_dict(){
        var repr_dict = this.getReprDictFrontEnd()
        this.send({
            // make sure we are using "request_repr_dict" name
            // in backend too.
            'type': 'request_repr_dict',
            'data': repr_dict,
        });
        var that = this
        if (that._synced_repr_model_ids.length > 0){
            that._synced_repr_model_ids.forEach(async function(mid){
                var model = await that.model.widget_manager.get_model(mid)
                    for (var k in model.views){
                        var view = await model.views[k];
                        // not sync with itself
                        if (view.uuid != that.uuid){
                            view._set_representation_from_repr_dict(repr_dict)
                        }
                    }
            })
        }
    }

    getReprDictFrontEnd(){
        var repr_dict = {};
        var n_components = this.stage.compList.length;
        for (var i = 0; i < n_components; i++) {
            var comp = this.stage.compList[i];
            repr_dict[i] = {};
            var msgi = repr_dict[i];
            for (var j = 0; j < comp.reprList.length; j++) {
                var repr = comp.reprList[j];
                msgi[j] = {};
                msgi[j]['type'] = repr.name;
                msgi[j]['params'] = repr.repr.getParameters();
            }
        }
        return repr_dict
    }

    syncReprForAllViews(){
        var repr_dict_backend = this.model.get("_ngl_repr_dict")
        var repr_dict_frontend = this.getReprDictFrontEnd()
        if (JSON.stringify(repr_dict_frontend) !== JSON.stringify(repr_dict_backend)){
            this._set_representation_from_repr_dict(repr_dict_backend)
        }
    }

    async syncReprWithMe(){
        // Make sure views of the same model has the same representations
        // Only needed if we use Sidebar that connects to specific view.
        var that = this
        var repr_dict = this.getReprDictFrontEnd()
        for (var k in this.model.views){
            var v = await this.model.views[k]
            if (v.uuid != that.uuid){
                v._set_representation_from_repr_dict(repr_dict)
            }
        }
        this.request_repr_dict()
    }

    setSyncRepr(model_ids){
        this._synced_repr_model_ids = model_ids
    }

    setSyncCamera(model_ids){
        this._synced_model_ids = model_ids
    }

    viewXZPlane(){
        var m = new NGL.Matrix4().makeRotationX( Math.PI / 2 );
        var q = new NGL.Quaternion().setFromRotationMatrix( m );
        this.stage.viewerControls.rotate( q );
    }

    set_representation_from_backend(){
        var repr_dict = this.model.get('_ngl_repr_dict')
        this._set_representation_from_repr_dict(repr_dict)
    }

    _set_representation_from_repr_dict(repr_dict){
        var compList = this.stage.compList
        if (compList.length > 0){
            for (var index in repr_dict){
                var comp = compList[index];
                comp.removeAllRepresentations();
                var reprlist = repr_dict[index];
                for (var j in reprlist){
                    var repr = reprlist[j];
                    if (repr){
                        comp.addRepresentation(repr.type, repr.params);
                    }
                }
            }
        }
    }

    async createView(trait_name){
        // Create a view for the model with given `trait_name`
        // e.g: in backend, 'view.<trait_name>`
        console.log("Creating view for model " + trait_name);
        var manager = this.model.widget_manager;
        var model_id = this.model.get(trait_name).replace("IPY_MODEL_", "");
        return await manager.create_view(await manager.get_model(model_id))
    }

    getPlayerModel(){
        // return a Promise
        var model_id = this.model.get("_iplayer").replace("IPY_MODEL_", "");
        return this.model.widget_manager.get_model(model_id)
    }

    async createIPlayer() {
        this.player_pview = this.createView("_iplayer");
        const view = await this.player_pview;
        const pe = view.el;
        pe.style.position = 'absolute';
        pe.style.zIndex = '100';
        pe.style.bottom = '5%';
        pe.style.left = '10%';
        pe.style.opacity = '0.7';
        this.stage.viewer.container.append(view.el);
        pe.style.display = 'none';
    }

    async createImageBtn() {
        this.image_btn_pview = this.createView("_ibtn_image");
        const view = await this.image_btn_pview;
        const pe = view.el;
        pe.style.position = 'absolute';
        pe.style.zIndex = '100';
        pe.style.top = '5%';
        pe.style.right = '10%';
        pe.style.opacity = '0.7';
        pe.style.width = '35px';
        this.stage.viewer.container.append(view.el);
    }

    async createFullscreenBtn() {
        this.btn_pview_fullscreen = this.createView("_ibtn_fullscreen");
        const view = await this.btn_pview_fullscreen;
        const pe = view.el;
        pe.style.position = 'absolute';
        pe.style.zIndex = '100';
        pe.style.top = '5%';
        pe.style.right = '5%';
        pe.style.opacity = '0.7';
        pe.style.width = '35px';
        pe.style.background = 'white';
        pe.style.opacity = '0.3';
        pe.style.display = 'none';
        pe.onclick = () => {
            this.stage.toggleFullscreen();
        };
        this.stage.viewer.container.append(view.el);
        this.stage.signals.fullscreenChanged.add((isFullscreen) => {
            if (isFullscreen) {
                view.model.set("icon", "compress");
            } else {
                view.model.set("icon", "expand");
            }
        });
    }

    async createGUI() {
        this.pgui_view = this.createView("_igui");
        const view = await this.pgui_view;
        const pe = view.el;
        pe.style.position = 'absolute';
        pe.style.zIndex = '100';
        pe.style.top = '5%';
        pe.style.right = '10%';
        pe.style.width = '300px';
        this.stage.viewer.container.append(view.el);
    }


    createNglGUI(){
      this.stage_widget = new StageWidget(this)
    }


        setVisibilityForRepr(component_index: number, repr_index: number, value: boolean) {
        const component = this.stage.compList[component_index];
        const repr = component?.reprList[repr_index];

        repr?.setVisibility(value);
    }

    removeRepresentation(component_index: number, repr_index: number) {
        const component = this.stage.compList[component_index];
        const repr = component?.reprList[repr_index];

        if (repr) {
            component.removeRepresentation(repr);
        }
    }

    removeRepresentationsByName(repr_name: string, component_index: number) {
        const component = this.stage.compList[component_index];

        component?.reprList.forEach((repr) => {
            if (repr.name === repr_name) {
                component.removeRepresentation(repr);
            }
        });
    }

    updateRepresentationForComponent(repr_index: number, component_index: number, params: any) {
        const component = this.stage.compList[component_index];
        const repr = component?.reprList[repr_index];

        repr?.setParameters(params);
    }

    updateRepresentationsByName(repr_name: string, component_index: number, params: any) {
        const component = this.stage.compList[component_index];

        component?.reprList.forEach((repr) => {
            if (repr.name === repr_name) {
                repr.setParameters(params);
                this.request_repr_dict();
            }
        });
    }

    setRepresentation(name: string, params: any, component_index: number, repr_index: number) {
        const component = this.stage.compList[component_index];
        const repr = component?.reprList[repr_index];

        if (repr) {
            params['useWorker'] = false;
            const new_repr = NGL.makeRepresentation(name, component.structure, this.stage.viewer, params);
            if (new_repr) {
                repr.setRepresentation(new_repr);
                repr.name = name;
                component.reprList[repr_index] = repr;
                this.request_repr_dict();
            }
        }
    }
    setColorByResidue(colors, component_index, repr_index){
        var repr = this.stage.compList[component_index].reprList[repr_index];
        var schemeId = NGL.ColormakerRegistry.addScheme(function(params){
            this.atomColor = function(atom){
                 var color = colors[atom.residueIndex];
                 return color
            };
            params; // to pass eslint; ack;
        });
        repr.setColor(schemeId);
    }

    addShape(name: string, shapes: any[]) {
        const shape = new NGL.Shape(name);
        const shapeDict = {
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

        shapes.forEach(([shapeType, ...params]) => {
            const func = shapeDict[shapeType];
            func.apply(shape, params);
        });

        const shapeComp = this.stage.addComponentFromObject(shape);
    }

    addBuffer(name, kwargs){
        var class_dict = {
            "arrow": NGL.ArrowBuffer,
            "box": NGL.BoXbuffer,
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
        for (var key in kwargs){
            params[key] = new Float32Array(kwargs[key]);
        }
        var shape = new NGL.Shape("shape");
        var buffer_class = class_dict[name];
        var buffer = new buffer_class(params);
        shape.addBuffer(buffer);
        var shapeComp = this.stage.addComponentFromObject(shape);
        shapeComp.addRepresentation("buffer");
    }

    async g(structure) {
        const { data, params = {}, ext } = structure;
        const blob = new Blob([data], { type: "text/plain" });
        const { compList, viewerControls } = this.stage;
        const [comp] = compList;
        const representations = comp.reprList.slice();
        const oldOrientation = viewerControls.getOrientation();

        const component = await this.stage.loadFile(blob, { ...params, ext, defaultRepresentation: false });

        viewerControls.orient(oldOrientation);

        representations.forEach(({ name, repr }) =>
            component.addRepresentation(name, repr.getParameters())
        );

        this.stage.removeComponent(comp);
        this._handleLoadFileFinished();
    }
    superpose(cindex0, cindex1, align, sele0, sele1) {
        // superpose two components with given params
        var component0 = this.stage.compList[cindex0];
        var component1 = this.stage.compList[cindex1];
        component1.superpose(component0, align, sele0, sele1);
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

    updateCoordinates(coordinates, model) {
        // coordinates must be ArrayBuffer (use this.decode_base64)
        var component = this.stage.compList[model];
        if (coordinates && component) {
            var coords = new Float32Array(coordinates);
            component.structure.updatePosition(coords);
            component.updateRepresentations({
                "position": true
            });
        }
    }

    handleResizable() {
        this.$container.resizable({
            resize: function(event, ui) {
                this.setSize(ui.size.width + "px", ui.size.height + "px");
            }.bind(this)
        })
    }

    handleResize(){
        var width = this.$el.width()
        var height = this.$el.height() + "px"
        if (this.stage_widget){
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

    GUIStyleChanged(){
        var style = this.model.get("gui_style");
        if (style === 'ngl'){
            this.createNglGUI();
        }else{
            if (this.stage_widget){
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
        if (this.ngl_view_id == this._getLastChildId()){
            const blob = await this.stage.makeImage(params);
            NGL.download(blob, filename);
        }
    }

    async _exportImage(wid, params) {
        if (this.ngl_view_id == this._getLastChildId()){
            const blob = await this.stage.makeImage(params);
            this.processBlob(blob, "image_data", wid);
        }
    }

    async handleMovieMaking(render_params) {
        console.log('handleMovieMaking: render_params', render_params);
        if (this.ngl_view_id == this._getLastChildId()){
            const blob = await this.stage.makeImage(render_params);
            this.processBlob(blob, "movie_image_data");
        }
    }

    async processBlob(blob, type, wid = null) {
        const reader = new FileReader();
        reader.onload = function() {
            const arr_str = (reader.result as string).replace("data:image/png;base64,", "");
            this.send({
                "data": arr_str,
                "type": type,
                "ID": wid,
            });
            this.send({'type': 'async_message', 'data': 'ok'});
        }.bind(this);
        reader.readAsDataURL(blob);
    }


    _handleLoadFileFinished() {
        this.send({'type': 'async_message', 'data': 'ok'});
    }

    _getLoadFilePromise(msg){
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

    _getLastChildId(){
        var keys = this.model.get('_ngl_view_id')
        return keys[keys.length-1]
    }

    async _handleStageLoadFile(msg){
        // args = [{'type': ..., 'data': ...}]
        if (this.ngl_view_id != this._getLastChildId() && msg.last_child){
            return
        }
        var o = await this._getLoadFilePromise(msg)
        this._handleLoadFileFinished();
    }

	addColorScheme(args, label){
        var id = NGL.ColormakerRegistry.addSelectionScheme(args, label);
        var scheme = NGL.ColormakerRegistry.userSchemes[id];
        NGL.ColormakerRegistry.removeScheme(id);
        // hard code the scheme ID
        NGL.ColormakerRegistry.add(label, scheme);
        return label
	}

    async on_msg(msg) {
        if (('ngl_view_id' in msg) && (msg.ngl_view_id !== this.ngl_view_id)) {
            return;
        }

        if (msg.type == 'call_method') {
            this.handleCallMethod(msg);
        } else if (msg.type == 'base64_single') {
            this.handleBase64Single(msg);
        } else if (msg.type == 'binary_single') {
            this.handleBinarySingle(msg);
        } else if (msg.type == 'get') {
            this.handleGet(msg);
        }
    }

    handleCallMethod(msg) {
        var new_args = msg.args.slice();
        new_args.push(msg.kwargs);

        this.handleColor(msg, new_args);

        switch (msg.target) {
            case Targets.Stage:
                this.handleStageTarget(msg, new_args);
                break;
            case Targets.Viewer:
                this.handleViewerTarget(msg, new_args);
                break;
            case Targets.ViewerControls:
                this.handleViewerControlsTarget(msg, new_args);
                break;
            case Targets.CompList:
                this.handleCompListTarget(msg, new_args);
                break;
            case Targets.Widget:
                this.handleWidgetTarget(msg, new_args);
                break;
            case Targets.Representation:
                this.handleRepresentationTarget(msg, new_args);
                break;
            default:
                console.log('there is no method for ' + msg.target);
                break;
        }
    }

    handleBase64Single(msg) {
        var coordinatesDict = msg.data;
        var keys = Object.keys(coordinatesDict);

        for (var i = 0; i < keys.length; i++) {
            var traj_index = keys[i];
            var coordinates = this.decode_base64(coordinatesDict[traj_index]);
            if (coordinates && coordinates.byteLength > 0) {
                this.updateCoordinates(coordinates, traj_index);
            }
        }
    }

    handleBinarySingle(msg) {
        var coordinateMeta = msg.data;
        var keys = Object.keys(coordinateMeta);

        for (var i = 0; i < keys.length; i++) {
            var traj_index = keys[i];
            var coordinates = new Float32Array(msg.buffers[i].buffer);
            if (coordinates.byteLength > 0) {
                this.updateCoordinates(coordinates, traj_index);
            }
        }
        if (msg.movie_making){
            this.handleMovieMaking(msg.render_params)
        }
    }

    handleGet(msg) {
        if (msg.data == 'camera') {
            this.send(JSON.stringify(this.stage.viewer.camera));
        } else if (msg.data == 'parameters') {
            this.send(JSON.stringify(this.stage.parameters));
        } else {
            console.log("Number of components", this.stage.compList.length);
            console.log("ngl_view_id", this.ngl_view_id);
        }
    }

    handleColor(msg, new_args) {
        if (msg.methodName == 'addRepresentation' && msg.reconstruc_color_scheme){
            msg.kwargs.color = this.addColorScheme(msg.kwargs.color, msg.kwargs.color_label);
        }
        if ("colorVolume" in msg.kwargs){
            // backend only send component index
            // so we need to convert to 'volume' data
            var index = msg.kwargs["colorVolume"];
            msg.kwargs["colorVolume"] = this.stage.compList[index].volume;
        }
    }

    async handleStageTarget(msg, new_args) {
        var stage_func = this.stage[msg.methodName];
        var stage = this.stage;
        if (msg.methodName == 'screenshot') {
            NGL.screenshot(this.stage.viewer, msg.kwargs);
        } else if (msg.methodName == 'removeComponent') {
            var index = msg.args[0];
            var component = this.stage.compList[index];
            this.stage.removeComponent(component);
        } else if (msg.methodName == 'loadFile') {
            if (this.model.views.length > 1 && msg.kwargs && msg.kwargs.defaultRepresentation) {
                // no need to add default representation as all representations
                // are serialized separately, also it unwantedly sets the orientation
                msg.kwargs.defaultRepresentation = false
            }
            await this._handleStageLoadFile(msg);
        } else {
            stage_func.apply(stage, new_args);
        }
    }

    handleViewerTarget(msg, new_args) {
        var viewer = this.stage.viewer;
        var func = this.stage.viewer[msg.methodName];
        func.apply(viewer, new_args);
    }

    handleViewerControlsTarget(msg, new_args) {
        var controls = this.stage.viewerControls;
        var func = controls[msg.methodName];
        func.apply(controls, new_args);
    }

    handleCompListTarget(msg, new_args) {
        var index = msg['component_index'];
        var component = this.stage.compList[index];
        var func = component[msg.methodName];
        func.apply(component, new_args);
    }

    handleWidgetTarget(msg, new_args) {
        var func = this[msg.methodName];
        if (func) {
            func.apply(this, new_args);
        } else {
            // send error message to Python?
            console.log('can not create func for ' + msg.methodName);
        }
    }

    handleRepresentationTarget(msg, new_args) {
        var component_index = msg['component_index'];
        var repr_index = msg['repr_index'];
        var component = this.stage.compList[component_index];
        var repr = component.reprList[repr_index];
        var func = repr[msg.methodName];
        if (repr && func) {
            func.apply(repr, new_args);
        }
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
}

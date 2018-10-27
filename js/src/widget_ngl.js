var widgets = require("@jupyter-widgets/base");
var NGL = require('ngl');
var $ = require('jquery');
var _ = require('underscore');
require('jquery-ui/ui/widgets/draggable');
require('jquery-ui/ui/widgets/slider');
require('jquery-ui/ui/widgets/dialog');
require('jquery-ui/themes/base/all.css');

var Jupyter;
if (typeof window !== 'undefined') {
  Jupyter = window['Jupyter'] = window['Jupyter'] || {};
} else {
  Jupyter = Jupyter || {};
}

var NGLModel = widgets.DOMWidgetModel.extend({
    defaults: function(){
        return _.extend(widgets.DOMWidgetModel.prototype.defaults(), {
            _model_name: 'NGLModel',
            _model_module: 'nglview-js-widgets',
            _model_module_version: require("../package.json").version,
            _view_name: "NGLView",
            _view_module: "nglview-js-widgets",
            _view_module_version: require("../package.json").version,
        });
    }
})

var NGLView = widgets.DOMWidgetView.extend({
    render: function() {
        // init setting of frame
        this.model.on("change:frame", this.frameChanged, this);

        // init setting of frame
        this.model.on("change:count", this.countChanged, this);

        // init _parameters handling
        this.model.on("change:_parameters", this.parametersChanged, this);

        this.model.set('_ngl_version', NGL.Version);

        // for player
        this.delay = 100;
        this.sync_frame = false;
        this.sync_camera = false;

        // get message from Python
        this.model.on("msg:custom", function(msg) {
            this.on_msg(msg);
        }, this);

        if (this.model.comm) {
            // for embeding in website
            this.model.comm.on_msg(function(msg) {
                var buffers = msg.buffers;
                var content = msg.content.data.content;
                if (buffers.length && content) {
                    content.buffers = buffers;
                }
                this.model._handle_comm_msg.call(this.model, msg);
            }.bind(this));
        }

        // init NGL stage
        NGL.useWorker = false;
        this.stage = new NGL.Stage(undefined, {
            backgroundColor: "white"
        });
        this.structureComponent = undefined;
        this.$container = $(this.stage.viewer.container);
        this.$el.append(this.$container);
        this.$container.resizable({
            resize: function(event, ui) {
                this.setSize(ui.size.width + "px", ui.size.height + "px");
            }.bind(this)
        });
        this.displayed.then(function() {
            this.ngl_view_id = this.get_last_child_id();
            var that = this;
            var width = this.$el.parent().width() + "px";
            var height = "300px";

            this.setSize(width, height);
            this.$container.resizable(
                "option", "maxWidth", this.$el.parent().width()
            );
            if (this.model.get("_ngl_serialize")){
                that.handle_embed();
            }else{
                this.requestUpdateStageParameters();
                if (this.model.views.length == 1){
                    this.serialize_camera_orientation();
                }else{
                    this.set_camera_orientation(that.model.get("_camera_orientation"));
                }
            }
        }.bind(this));

        this.stage.viewerControls.signals.changed.add(function() {
            this.serialize_camera_orientation();
        }.bind(this));

        // init toggle fullscreen
        $(this.stage.viewer.container).dblclick(function(e) {
            if (!e.altKey || !e.ctrlKey){
                this.stage.toggleFullscreen();
            }
        }.bind(this));

        // init picking handling
        this.$pickingInfo = $("<div></div>")
            .css("position", "absolute")
            .css("top", "5%")
            .css("left", "3%")
            .css("background-color", "white")
            .css("padding", "2px 5px 2px 5px")
            .css("opacity", "0.7")
            .appendTo(this.$container);

        var $inputNotebookCommand = $('<input id="input_notebook_command" type="text" style="border:1px solid skyblue" size="50"></input>');
        var that = this;

        $inputNotebookCommand.keypress(function(e) {
            var command = $("#input_notebook_command").val();
            if (e.which == 13) {
                $("#input_notebook_command").val("")
                Jupyter.notebook.kernel.execute(command);
            }
        });

        this.$notebook_text = $("<div></div>")
            .css("position", "absolute")
            .css("bottom", "5%")
            .css("left", "3%")
            .css("padding", "2px 5px 2px 5px")
            .css("opacity", "0.7")
            .append($inputNotebookCommand)
            .appendTo(this.$container);
        this.$notebook_text.hide();

        this.stage.signals.clicked.add(function(pd) {
            if (pd) {
                var pd2 = {};
                if (pd.atom) pd2.atom = pd.atom.toObject();
                if (pd.bond) pd2.bond = pd.bond.toObject();
                if (pd.instance) pd2.instance = pd.instance;
                this.model.set("picked", pd2);
                this.touch();
                var pickingText = "";
                if (pd.atom) {
                    pickingText = "Atom: " + pd.atom.qualifiedName();
                } else if (pd.bond) {
                    pickingText = "Bond: " + pd.bond.atom1.qualifiedName() + " - " + pd.bond.atom2.qualifiedName();
                }
                this.$pickingInfo.text(pickingText);
            }
        }, this);

        this.initPlayer();

        var container = this.stage.viewer.container;
        that = this;
        container.addEventListener('dragover', function(e) {
            e.stopPropagation();
            e.preventDefault();
            e.dataTransfer.dropEffect = 'copy';
        }, false);

        container.addEventListener('drop', function(e) {
            e.stopPropagation();
            e.preventDefault();
            var file = e.dataTransfer.files[0];

            that.stage.loadFile(file).then(function(o){
                that._handle_loading_file_finished();
                o;
            });
            var numDroppedFiles = that.model.get("_n_dragged_files");
            that.model.set("_n_dragged_files", numDroppedFiles + 1);
            that.touch();
        }, false);

        that = this;
        this.stage.signals.componentAdded.add(function() {
            var len = this.stage.compList.length;
            this.model.set("n_components", len);
            this.touch();
            var comp = this.stage.compList[len - 1];
            comp.signals.representationRemoved.add(function() {
                that.request_repr_dict();
            });
            comp.signals.representationAdded.add(function() {
                that.request_repr_dict();
            });
        }, this);

        this.stage.signals.componentRemoved.add(function() {
            this.model.set("n_components", this.stage.compList.length);
            this.touch();
        }, this);

        // for callbacks from Python
        // must be after initializing NGL.Stage
        this.send({
            'type': 'request_loaded',
            'data': true
        })
        var state_params = this.stage.getParameters();
        this.model.set('_ngl_original_stage_parameters', state_params);
        this.touch();
    },

    serialize_camera_orientation: function(){
         var m = this.stage.viewerControls.getOrientation();
         this.model.set('_camera_orientation', m.elements);
         this.touch();
    },

    set_camera_orientation: function(orientation){
        if (orientation.length > 0){
            this.stage.viewerControls.orient(orientation);
            this.serialize_camera_orientation();
        }
    },

    handle_embed: function(){
        var that = this;
        var ngl_coordinate_resource = that.model.get("_ngl_coordinate_resource");
        var ngl_msg_archive = that.model.get("_ngl_msg_archive");
        var ngl_stage_params = that.model.get('_ngl_full_stage_parameters_embed');
        var loadfile_list = [];

        _.each(ngl_msg_archive, function(msg){
            if (msg.methodName == 'loadFile'){
                if (msg.kwargs && msg.kwargs.defaultRepresentation) {
                    // no need to add default representation as all representations
                    // are serialized separately, also it unwantedly sets the orientation
                    msg.kwargs.defaultRepresentation = false
                 }
                loadfile_list.push(that._get_loadFile_promise(msg));
            }
        });


        Promise.all(loadfile_list).then(function(compList){
            that._set_representation_from_backend(compList);
            that.stage.setParameters(ngl_stage_params);
            that.set_camera_orientation(that.model.get("_camera_orientation"));

            var frame = 0;
            var count = ngl_coordinate_resource['n_frames'];
            delete ngl_coordinate_resource['n_frames'];

            var play = function(){
                that.$playerButton.text("pause");
                that.playerInterval = setInterval(function(){
                    frame = frame + 1;
                    if (frame > count - 1){
                        frame = 0;
                    }
                    that.$playerSlider.slider("option", "value", frame);
                    that.updateCoordinatesFromDict(ngl_coordinate_resource, frame);
                }, that.delay)
            }

            var pause = function() {
                that.$playerButton.text("play");
                if (that.playerInterval !== undefined) {
                    clearInterval(that.playerInterval);
                }
            }.bind(that);

            if (that.$playerButton){
                that.$playerButton
                    .off('click')
                    .click(function(event) {
                        if (that.$playerButton.text() === "play") {
                            play();
                        } else if (that.$playerButton.text() === "pause") {
                            pause();
                        }
                        event; // to pass eslint
                    }.bind(that));
                that.$playerSlider.slider({
                    max : count-1,
                    slide: function(event, ui) {
                        pause();
                        that.updateCoordinatesFromDict(ngl_coordinate_resource, ui.value);
                        frame = ui.value;
                    }.bind(that)
                })
           }
        });
    },

    updateCoordinatesFromDict: function(cdict, frame_index){
        // update coordinates for given "index"
        // cdict = Dict[int, List[base64]]
        var keys = Object.keys(cdict);

        for (var i = 0; i < keys.length; i++) {
            var traj_index = keys[i];
            var coordinates = this.decode_base64(cdict[traj_index][frame_index]);
            if (coordinates && coordinates.byteLength > 0) {
                this.updateCoordinates(coordinates, traj_index);
            }
        }
    },

    setSelector: function(selector_id) {
        // id is uuid that will be set from Python
        var selector = "<div class='" + selector_id + "'></div>";
        this.$ngl_selector = $(selector)
            .css("position", "absolute")
            .css("bottom", "5%")
            .css("left", "3%")
            .css("padding", "2px 5px 2px 5px")
            .css("opacity", "0.7")
            .appendTo(this.$container);
    },

    setIPythonLikeCell: function() {
        var cell = Jupyter.notebook.insert_cell_at_bottom();

        var handler = function(event) {
            var selected_cell = Jupyter.notebook.get_selected_cell();
            if (selected_cell.cell_id === cell.cell_id) {
                selected_cell.execute();
                selected_cell.set_text('');
            }
            event; // to pass eslint
            return false;
        };

        var action = {
            help: 'run cell',
            help_index: 'zz',
            handler: handler
        };

        Jupyter.keyboard_manager.edit_shortcuts.add_shortcut('enter', action);
    },

    hideNotebookCommandBox: function() {
        this.$notebook_text.hide();
    },

    showNotebookCommandBox: function() {
        this.$notebook_text.show();
    },

    requestFrame: function() {
        this.send({
            'type': 'request_frame',
            'data': 'frame'
        });
    },

    requestUpdateStageParameters: function() {
        var updated_params = this.stage.getParameters();
        this.model.set('_ngl_full_stage_parameters', updated_params);
        this.touch();
    },

    requestReprParameters: function(component_index, repr_index) {
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
    },

    request_repr_dict: function() {
        var n_components = this.stage.compList.length;
        var msg = {};

        for (var i = 0; i < n_components; i++) {
            var comp = this.stage.compList[i];
            msg[i] = {};
            var msgi = msg[i];
            for (var j = 0; j < comp.reprList.length; j++) {
                var repr = comp.reprList[j];
                msgi[j] = {};
                msgi[j]['type'] = repr.name;
                msgi[j]['params'] = repr.repr.getParameters();
            }
        }
        this.send({
            // make sure we are using "request_repr_dict" name
            // in backend too.
            'type': 'request_repr_dict',
            'data': msg
        });
    },

    // setDraggable: function(params) {
    //     if (params) {
    //         this.$container.draggable(params);
    //     } else {
    //         this.$container.draggable();
    //     }
    // },
    setDelay: function(delay) {
        this.delay = delay;
    },

    setSyncFrame: function() {
        this.sync_frame = true;
    },

    setUnSyncFrame: function() {
        this.sync_frame = false;
    },

    setSyncCamera: function() {
        this.sync_camera = true;
        this.serialize_camera_orientation();
    },

    setUnSyncCamera: function() {
        this.sync_camera = false;
    },

    viewXZPlane: function() {
        var m = new NGL.Matrix4().makeRotationX( Math.PI / 2 );
        var q = new NGL.Quaternion().setFromRotationMatrix( m );
        this.stage.viewerControls.rotate( q );
    },

    set_representation_from_backend: function(){
        this._set_representation_from_backend(this.stage.compList);
    },

    _set_representation_from_backend: function(compList){
        if (compList.length > 0){
            var ngl_repr_dict = this.model.get('_ngl_repr_dict');
            for (var index in ngl_repr_dict){
                var comp = compList[index];
                comp.removeAllRepresentations();
                var reprlist = ngl_repr_dict[index];
                for (var j in reprlist){
                    var repr = reprlist[j];
                    if (repr){
                        comp.addRepresentation(repr.type, repr.params);
                    }
                }
            }
        }
    },

    initPlayer: function() {
        // init player
        if (this.model.get("count")) {
            var frame = this.model.get("frame");
            var play = function() {
                this.$playerButton.text("pause");
                this.playerInterval = setInterval(function() {
                    var frame = this.model.get("frame") + 1;
                    var count = this.model.get("count");
                    if (frame >= count) frame = 0;

                    if (this.sync_frame) {
                        this.model.set("frame", frame);
                        this.touch();
                    } else {
                        this.requestFrame();
                    }
                }.bind(this), this.delay);
            }.bind(this);
            var pause = function() {
                this.$playerButton.text("play");
                if (this.playerInterval !== undefined) {
                    clearInterval(this.playerInterval);
                }
            }.bind(this);
            this.$playerButton = $("<button>play</button>")
                .css("float", "left")
                .css("width", "55px")
                .css("opacity", "0.7")
                .click(function(event) {
                    if (this.$playerButton.text() === "play") {
                        play();
                    } else if (this.$playerButton.text() === "pause") {
                        pause();
                    }
                    event; // to pass eslint
                }.bind(this));
            this.$playerSlider = $("<div></div>")
                .css("margin-left", "70px")
                .css("position", "relative")
                .css("bottom", "-7px")
                .slider({
                    min: 0,
                    max: this.model.get("count") - 1,
                    value: frame,
                    slide: function(event, ui) {
                        pause();
                        this.model.set("frame", ui.value);
                        this.touch();
                    }.bind(this)
                });
            this.$player = $("<div></div>")
                .css("position", "absolute")
                .css("bottom", "5%")
                .css("width", "94%")
                .css("margin-left", "3%")
                .css("opacity", "0.7")
                .append(this.$playerButton)
                .append(this.$playerSlider)
                .appendTo(this.$container);
            this.model.on("change:frame", function() {
                this.$playerSlider.slider("value", this.model.get("frame"));
            }, this);

            if (this.model.get("count") < 2) {
                this.$player.hide()
            }
        }
    },

    countChanged: function() {
        var count = this.model.get("count");
        this.$playerSlider.slider({
            max: count - 1
        });
        if (this.model.get("count") > 1) {
            this.$player.show()
        }
    },

    setVisibilityForRepr: function(component_index, repr_index, value) {
        // value = True/False
        var component = this.stage.compList[component_index];
        var repr = component.reprList[repr_index];

        if (repr) {
            repr.setVisibility(value);
        }
    },

    removeRepresentation: function(component_index, repr_index) {
        var component = this.stage.compList[component_index];
        var repr = component.reprList[repr_index]

        if (repr) {
            component.removeRepresentation(repr);
        }
    },

    removeRepresentationsByName: function(repr_name, component_index) {
        var component = this.stage.compList[component_index];

        if (component) {
            component.reprList.forEach(function(repr) {
                if (repr.name == repr_name) {
                    component.removeRepresentation(repr);
                }
            })
        }
    },

    updateRepresentationForComponent: function(repr_index, component_index, params) {
        var component = this.stage.compList[component_index];
        var that = this;
        var repr = component.reprList[repr_index];
        if (repr) {
            repr.setParameters(params);
            that.request_repr_dict();
        }
    },

    updateRepresentationsByName: function(repr_name, component_index, params) {
        var component = this.stage.compList[component_index];
        var that = this;

        if (component) {
            component.reprList.forEach(function(repr) {
                if (repr.name == repr_name) {
                    repr.setParameters(params);
                    that.request_repr_dict();
                }
            })
        }
    },

    setRepresentation: function(name, params, component_index, repr_index) {
        var component = this.stage.compList[component_index];
        var repr = component.reprList[repr_index];
        var that = this;

        if (repr) {
            params['useWorker'] = false;
            var new_repr = NGL.makeRepresentation(name, component.structure,
                this.stage.viewer, params);
            if (new_repr) {
                repr.setRepresentation(new_repr);
                repr.name = name;
                component.reprList[repr_index] = repr;
                that.request_repr_dict();
            }
        }
    },

    setColorByResidue: function(colors, component_index, repr_index){
        console.log('NGL.ColormakerRegistry', NGL.ColormakerRegistry);
        var repr = this.stage.compList[component_index].reprList[repr_index];
        var schemeId = NGL.ColormakerRegistry.addScheme(function(params){
            this.atomColor = function(atom){
                 var color = colors[atom.residueIndex];
                 return color
            };
            params; // to pass eslint; ack;
        });
        repr.setColor(schemeId);
    },

    addShape: function(name, shapes) {
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
        shapeComp.addRepresentation("buffer");
    },

    addBuffer: function(name, kwargs){
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
    },

    replaceStructure: function(structure){
         var blob = new Blob([structure.data], {type: "text/plain"});
         var params = structure.params || {};
         params.ext = structure.ext;
         params.defaultRepresentation = false;
         var comp = this.stage.compList[0];
         var representations = comp.reprList.slice();
         var old_orientation = this.stage.viewerControls.getOrientation();
         var that = this;
         this.stage.loadFile(blob, params).then(function(component) {
             that.stage.viewerControls.orient(old_orientation);
             representations.forEach(function(repr) {
                 var repr_name = repr.name;
                 var repr_params = repr.repr.getParameters();
                 // Note: not using repr.repr.type, repr.repr.params
                 // since seems to me that repr.repr.params won't return correct "sele"
                 component.addRepresentation(repr_name, repr_params);
             });
             that.stage.removeComponent(comp);
             that._handle_loading_file_finished();
         });
    },

    superpose: function(cindex0, cindex1, align, sele0, sele1) {
        // superpose two components with given params
        var component0 = this.stage.compList[cindex0];
        var component1 = this.stage.compList[cindex1];
        component1.superpose(component0, align, sele0, sele1);
    },

    decode_base64: function(base64) {
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
    },

    updateCoordinates: function(coordinates, model) {
        // coordinates must be ArrayBuffer (use this.decode_base64)
        var component = this.stage.compList[model];
        if (coordinates && component) {
            var coords = new Float32Array(coordinates);
            component.structure.updatePosition(coords);
            component.updateRepresentations({
                "position": true
            });
        }
    },

    handleResize: function() {
        this.$container.resizable({
            resize: function(event, ui) {
                this.setSize(ui.size.width + "px", ui.size.height + "px");
            }.bind(this)
        })
    },

    setSize: function(width, height) {
        this.stage.viewer.container.style.width = width;
        this.stage.viewer.container.style.height = height;
        this.stage.handleResize();
    },

    openNotebookCommandDialog: function() {
        var that = this;
        var dialog = this.$notebook_text.dialog({
            draggable: true,
            resizable: true,
            modal: false,
            show: {
                effect: "blind",
                duration: 150
            },
            close: function(event, ui) {
                that.$container.append(that.$notebook_text);
                that.$notebook_text.dialog('destroy');
                event; ui; // to pass eslint; ack;
            },
        });
        dialog.css({
            overflow: 'hidden'
        });
        dialog.prev('.ui-dialog-titlebar')
            .css({
                'background': 'transparent',
                'border': 'none'
            });
        Jupyter.keyboard_manager.register_events(dialog);
    },

    setDialog: function() {
        var $nb_container = Jupyter.notebook.container;
        var that = this;
        var dialog = this.$container.dialog({
            title: "NGLView",
            draggable: true,
            resizable: true,
            modal: false,
            width: window.innerWidth - $nb_container.width() - $nb_container.offset().left - 50,
            height: 'auto',
            position: {
                my: 'right',
                at: 'right',
                of: window
            },
            show: {
                effect: "blind",
                duration: 150
            },
            close: function(event, ui) {
                that.$el.append(that.$container);
                that.$container.dialog('destroy');
                that.handleResize();
                event; ui; // to pass eslint; ack;
            },
            resize: function(event, ui) {
                that.stage.handleResize();
                that.setSize(ui.size.width + "px", ui.size.height + "px");
            }.bind(that),
        });
        dialog.css({
            overflow: 'hidden'
        });
        dialog.prev('.ui-dialog-titlebar')
            .css({
                'background': 'transparent',
                'border': 'none'
            });
    },

    resizeNotebook: function(width) {
        var $nb_container = Jupyter.notebook.container;
        $nb_container.width(width);

        if (this.$container.dialog) {
            this.$container.dialog({
                width: $nb_container.offset().left
            });
        }
    },

    parametersChanged: function() {
        var _parameters = this.model.get("_parameters");
        this.setParameters(_parameters);
    },

    setParameters: function(parameters) {
        this.stage.setParameters(parameters);

        // do not set _ngl_full_stage_parameters here
        // or parameters will be never updated (not sure why)
        // use observe in python side
        var updated_params = this.stage.getParameters();
        this.send({
            'type': 'stage_parameters',
            'data': updated_params
        })
    },

    _downloadImage: function(filename, params) {
        if (this.ngl_view_id == this.get_last_child_id()){
            this.stage.makeImage(params).then(function(blob) {
                NGL.download(blob, filename);
            })
        }
    },

    _exportImage: function(params) {
        if (this.ngl_view_id == this.get_last_child_id()){
            this.stage.makeImage(params).then(function(blob) {
                var reader = new FileReader();
                var arr_str;
                reader.onload = function() {
                    arr_str = reader.result.replace("data:image/png;base64,", "");
                    this.model.set("_image_data", arr_str);
                    this.touch();
                }.bind(this);
                reader.readAsDataURL(blob);
            }.bind(this));
        }
    },

    cleanOutput: function() {

        var cells = Jupyter.notebook.get_cells();

        for (var i = 0; i < cells.length; i++) {
            var cell = cells[i];
            if (cell.output_area.outputs.length > 0) {
                var out = cell.output_area.outputs[0];
                if (out.output_type == 'display_data') {
                    cell.clear_output();
                }
            }
        }
    },

    _handle_loading_file_finished: function() {
        this.send({'type': 'async_message', 'data': 'ok'});
    },

    _get_loadFile_promise: function(msg){
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
             return this.stage.loadFile(msg.args[0].data, msg.kwargs)
         }
    },

    get_last_child_id: function(){
        var keys = Object.keys(this.model.views);
        return keys[keys.length-1]
    },

    _handle_stage_loadFile: function(msg){
        // args = [{'type': ..., 'data': ...}]
        if (this.ngl_view_id != this.get_last_child_id() && msg.last_child){
            return
        }
        var that = this;
        this._get_loadFile_promise(msg).then(function(o){
            that._handle_loading_file_finished();
            o;
        });
    },

    _make_color_scheme: function(args){
        console.log("making color scheme", args);
        return NGL.ColormakerRegistry.addSelectionScheme(args)
    },

    on_msg: function(msg) {
        // TODO: re-organize
        if (msg.type == 'call_method') {
            var index, component, func, stage;
            var new_args = msg.args.slice();
            new_args.push(msg.kwargs);

            if (msg.methodName == 'addRepresentation' && 
                msg.reconstruc_color_scheme){
                msg.kwargs.color = this._make_color_scheme(msg.kwargs.color);
            };

            switch (msg.target) {
                case 'Stage':
                    var stage_func = this.stage[msg.methodName];
                    stage = this.stage;
                    if (msg.methodName == 'screenshot') {
                        NGL.screenshot(this.stage.viewer, msg.kwargs);
                    } else if (msg.methodName == 'removeComponent') {
                        index = msg.args[0];
                        component = this.stage.compList[index];
                        this.stage.removeComponent(component);
                    } else if (msg.methodName == 'loadFile') {
                        if (this.model.views.length > 1 && msg.kwargs &&
                            msg.kwargs.defaultRepresentation) {
                            // no need to add default representation as all representations
                            // are serialized separately, also it unwantedly sets the orientation
                            msg.kwargs.defaultRepresentation = false
                        }
                        this._handle_stage_loadFile(msg);
                    } else {
                            stage_func.apply(stage, new_args);
                    }
                    break;
                case 'Viewer':
                    var viewer = this.stage.viewer;
                    func = this.stage.viewer[msg.methodName];
                    func.apply(viewer, new_args);
                    break;
                case 'viewerControls':
                    var controls = this.stage.viewerControls;
                    func = controls[msg.methodName];
                    func.apply(controls, new_args);
                    break;
                case 'compList':
                    index = msg['component_index'];
                    component = this.stage.compList[index];
                    func = component[msg.methodName];
                    func.apply(component, new_args);
                    break;
                case 'StructureComponent':
                    component = this.structureComponent;
                    func = component[msg.methodName];
                    func.apply(component, new_args);
                    break;
                case 'Widget':
                    func = this[msg.methodName];
                    if (func) {
                        func.apply(this, new_args);
                    } else {
                        // send error message to Python?
                        console.log('can not create func for ' + msg.methodName);
                    }
                    break;
                case 'player':
                    func = this.$player[msg.methodName];
                    if (func) {
                        func.apply(this.$player, new_args);
                    } else {
                        // send error message to Python?
                        console.log('can not create func for ' + msg.methodName);
                    }
                    break;
                case 'Representation':
                    var component_index = msg['component_index'];
                    var repr_index = msg['repr_index'];
                    component = this.stage.compList[component_index];
                    var repr = component.reprList[repr_index];
                    func = repr[msg.methodName];
                    if (repr && func) {
                        func.apply(repr, new_args);
                    }
                    break;
                default:
                    console.log('there is no method for ' + msg.target);
                    break;
            }
        } else if (msg.type == 'base64_single') {
            var coordinatesDict = msg.data;
            var keys = Object.keys(coordinatesDict);

            for (var i = 0; i < keys.length; i++) {
                var traj_index = keys[i];
                var coordinates = this.decode_base64(coordinatesDict[traj_index]);
                if (coordinates && coordinates.byteLength > 0) {
                    this.updateCoordinates(coordinates, traj_index);
                }
            }
        } else if (msg.type == 'binary_single') {
            var coordinateMeta = msg.data;
            keys = Object.keys(coordinateMeta);

            for (i = 0; i < keys.length; i++) {
                traj_index = keys[i];
                coordinates = new Float32Array(msg.buffers[i].buffer);
                if (coordinates.byteLength > 0) {
                    this.updateCoordinates(coordinates, traj_index);
                }
            }
        } else if (msg.type == 'get') {
            if (msg.data == 'camera') {
                this.send(JSON.stringify(this.stage.viewer.camera));
            } else if (msg.data == 'parameters') {
                this.send(JSON.stringify(this.stage.parameters));
            } else {
                console.log("Number of components", this.stage.compList.length);
                console.log("ngl_view_id", this.ngl_view_id);
            }
        }
    },
});

module.exports = {
    'NGLView': NGLView,
    'NGLModel': NGLModel,
    'NGL': NGL,
};

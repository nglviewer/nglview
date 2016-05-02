
require.config( {
    paths: {
        "THREE": "../nbextensions/nglview/three.custom.min",
        "Detector": "../nbextensions/nglview/Detector",
        "async": "../nbextensions/nglview/async.min",
        "Promise": "../nbextensions/nglview/promise.min",
        "sprintf": "../nbextensions/nglview/sprintf.min",
        "JSZip": "../nbextensions/nglview/jszip.min",
        "pako": "../nbextensions/nglview/pako.min",
        "LZMA": "../nbextensions/nglview/lzma.min",
        "bzip2": "../nbextensions/nglview/bzip2.min",
        "chroma": "../nbextensions/nglview/chroma.min",
        "jsfeat": "../nbextensions/nglview/svd.min",
        "signals": "../nbextensions/nglview/signals.min",
        "NGL": "../nbextensions/nglview/ngl",
        "mdsrv": "../nbextensions/nglview/mdsrv",
    },
    shim: {
        THREE: { exports: "THREE" },
        Detector: { exports: "Detector" },
        Promise: { exports: "Promise" },
        sprintf: { exports: "sprintf" },
        JSZip: { exports: "JSZip" },
        LZMA: { exports: "LZMA" },
        bzip2: { exports: "bzip2" },
        jsfeat: { exports: "jsfeat" },
        NGL: {
            deps: [ "THREE", "Detector", "async", "Promise", "sprintf", "JSZip", "pako", "LZMA", "bzip2", "chroma", "jsfeat", "signals" ],
            exports: "NGL"
        }
    },
} );

define( [
    "nbextensions/widgets/widgets/js/widget",
    "nbextensions/widgets/widgets/js/manager",
    "jqueryui",
    "THREE", "Detector", "async", "Promise", "sprintf", "JSZip", "pako",
    "LZMA", "bzip2", "chroma", "jsfeat", "signals", "NGL", "mdsrv"
], function(
    widget, manager, $, _THREE, _Detector, async, _Promise, _sprintf, _JSZip, pako,
    _LZMA, _bzip2, chroma, _jsfeat, signals, _NGL, _NGL_mdsrv
){

    window.async = async;
    window.pako = pako;
    window.signals = signals;
    window.chroma = chroma;

    var NGLView = widget.DOMWidgetView.extend( {

        render: function(){

            // console.log( "NGLView", this )

            if( !Detector.webgl ){
                Detector.addGetWebGLMessage( { "parent": this.el } );
                return;
            }

            NGL.init( function(){

                // init representations handling
                this.model.on( "change:_init_representations", this.representationsChanged, this );

                // init structure loading
                this.model.on( "change:structure", this.structureChanged, this );

                // init setting of coordinates
                this.model.on( "change:_coordinates_meta", this.coordinatesChanged, this );

                // init setting of coordinates
                this.model.on( "change:coordinates_dict", this.coordsDictChanged, this );

                // init setting of frame
                this.model.on( "change:frame", this.frameChanged, this );

                // init parameters handling
                this.model.on( "change:parameters", this.parametersChanged, this );

                // init parameters handling
                this.model.on( "change:cache", this.cacheChanged, this );

                // get message from Python
                this.coordsDict = {};
                this.model.on( "msg:custom", function (msg) {
                    this.on_msg( msg );
                }, this);

                // init NGL stage
                NGL.useWorker = false;
                this.stage = new NGL.Stage();
                this.stage.setTheme( "light" );
                this.structureComponent = undefined;
                this.$container = $( this.stage.viewer.container );
                this.$el.append( this.$container );
                this.$container.resizable( {
                    resize: function( event, ui ){
                        this.setSize( ui.size.width + "px", ui.size.height + "px" );
                    }.bind( this )
                } );
                this.displayed.then( function(){
                    var width = this.$el.parent().width() + "px";
                    var height = "300px";
                    this.setSize( width, height );
                    this.$container.resizable(
                        "option", "maxWidth", this.$el.parent().width()
                    );
                }.bind( this ) );

                // init toggle fullscreen
                $( this.stage.viewer.container ).dblclick( function(){
                    this.stage.toggleFullscreen();
                }.bind( this ) );

                // init model data
                this.structureChanged();

                // init picking handling
                this.$pickingInfo = $( "<div></div>" )
                    .css( "position", "absolute" )
                    .css( "top", "5%" )
                    .css( "left", "3%" )
                    .css( "background-color", "white" )
                    .css( "padding", "2px 5px 2px 5px" )
                    .css( "opacity", "0.7" )
                    .appendTo( this.$container );
                this.stage.signals.onPicking.add( function( pd ){
                    var pd2 = {};
                    if( pd.atom ) pd2.atom = pd.atom.toJSON();
                    if( pd.bond ) pd2.bond = pd.bond.toJSON();
                    if( pd.instance ) pd2.instance = pd.instance;
                    this.model.set( "picked", pd2 );
                    this.model.save();
                    var pickingText = "";
                    if( pd.atom ){
                        pickingText = "Atom: " + pd.atom.qualifiedName();
                    }else if( pd.bond ){
                        pickingText = "Bond: " + pd.bond.atom1.qualifiedName() + " - " + pd.bond.atom2.qualifiedName();
                    }
                    this.$pickingInfo.text( pickingText );
                }, this );

                // init player
                if( this.model.get( "count" ) ){
                    var play = function(){
                        this.$playerButton.text( "pause" );
                        this.playerInterval = setInterval( function(){
                            var frame = this.model.get( "frame" ) + 1;
                            var count = this.model.get( "count" );
                            if( frame >= count ) frame = 0;
                            this.model.set( "frame", frame );
                            this.model.save();
                        }.bind( this ), 100 );
                    }.bind( this );
                    var pause = function(){
                        this.$playerButton.text( "play" );
                        if( this.playerInterval !== undefined ){
                            clearInterval( this.playerInterval );
                        }
                    }.bind( this );
                    this.$playerButton = $( "<button>play</button>" )
                        .css( "float", "left" )
                        .css( "width", "55px" )
                        .css( "opacity", "0.7" )
                        .click( function( event ){
                            if( this.$playerButton.text() === "play" ){
                                play();
                            }else if( this.$playerButton.text() === "pause" ){
                                pause();
                            }
                        }.bind( this ) );
                    this.$playerSlider = $( "<div></div>" )
                        .css( "margin-left", "70px" )
                        .css( "position", "relative" )
                        .css( "bottom", "-7px" )
                        .slider( {
                            min: 0,
                            max: this.model.get( "count" ) - 1,
                            slide: function( event, ui ){
                                pause();
                                this.model.set( "frame", ui.value );
                                this.model.save();
                            }.bind( this )
                        } );
                    this.$player = $( "<div></div>" )
                        .css( "position", "absolute" )
                        .css( "bottom", "5%" )
                        .css( "width", "94%" )
                        .css( "margin-left", "3%" )
                        .css( "opacity", "0.7" )
                        .append( this.$playerButton )
                        .append( this.$playerSlider )
                        .appendTo( this.$container );
                    this.model.on( "change:frame", function(){
                        this.$playerSlider.slider( "value", this.model.get( "frame" ) );
                    }, this );
                }

            }.bind( this ) );

        },

        representationsChanged: function(){
            var representations = this.model.get( "_init_representations" );
            var component = this.structureComponent;
            if( representations && component ){
                component.clearRepresentations();
                representations.forEach( function( repr ){
                    component.addRepresentation( repr.type, repr.params );
                } );
            }
        },

        structureChanged: function(){
            this.structureComponent = undefined;
            var structure = this.model.get( "structure" );
            if( structure.data && structure.ext ){
                var blob = new Blob( [ structure.data ], { type: "text/plain" } );
                var params = structure.params || {};
                params.ext = structure.ext;
                params.defaultRepresentation = false;
                this.stage.loadFile( blob, params ).then( function( component ){
                    component.centerView();
                    this.structureComponent = component;
                    this.representationsChanged();
                }.bind( this ) );
            }
        },

        mydecode: function(base64) {
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

        frameChanged: function(){
            if( this._cache ){
                var frame = this.model.get( "frame" );
                if( frame in this.coordsDict ) {
                    var coordinates = this.coordsDict[frame];
                    this._update_coords(coordinates);
                } // else: just wait
            }
            // else: listen to coordinatesChanged
        },


        coordinatesChanged: function(){
            if (! this._cache ){
                var coordinates_meta = this.model.get( "_coordinates_meta" );

                // not checking dtype yet
                var coordinates = this.mydecode( coordinates_meta['data'] );
                this._update_coords(coordinates);
            }
        },

        _update_coords: function( coordinates ) {
            // coordinates must be ArrayBuffer (use this.mydecode)
            var component = this.structureComponent;
            if( coordinates && component ){
                var coords = new Float32Array( coordinates );
                component.structure.updatePosition( coords );
                component.updateRepresentations( { "position": true } );
            }
        },

        coordsDictChanged: function(){
            this.coordsDict = this.model.get( "coordinates_dict" );
            var cdict = this.coordsDict
            var clen = Object.keys(cdict).length
            if ( clen != 0 ){
                this._cache = true;
            }else{
                this._cache = false;
            }
            this.model.set( "cache", this._cache);
            
            for (var i = 0; i < Object.keys(coordsDict).length; i++) {
                this.coordsDict[i] = this.mydecode( coordsDict[i]);
            }
        },

        setSize: function( width, height ){
            this.stage.viewer.container.style.width = width;
            this.stage.viewer.container.style.height = height;
            this.stage.handleResize();
        },

        parametersChanged: function(){
            var parameters = this.model.get( "parameters" );
            this.stage.setParameters( parameters );
        },

        cacheChanged: function(){
            this._cache = this.model.get( "cache" );
        },

        on_msg: function(msg){
            // TODO: re-organize
            if( msg.type == 'call_method' ){
                var new_args = msg.args.slice();
                new_args.push( msg.kwargs );

                switch( msg.target ) {
                    case 'Stage':
                        var stage_func = this.stage[msg.methodName];
                        var stage = this.stage;
                        if ( msg.methodName == 'screenshot' ){
                             NGL.screenshot( this.stage.viewer, msg.kwargs );
                        }else if( msg.methodName == 'removeComponent' ) {
                            console.log( "removeComponent" );
                            var index = msg.args[0];
                            var component = this.stage.compList[ index ];
                            this.stage.removeComponent( component );
                        }else{
                            if( msg.methodName == 'loadFile' ) {
                                // args = [{'type': ..., 'data': ...}]
                                if( msg.args[0].type == 'blob' ) {
                                    var blob = new Blob( [ msg.args[0].data ], { type: "text/plain" } );
                                    this.stage.loadFile( blob, msg.kwargs );
                                }else if( msg.args[0].type == 'path' ) {
                                    this.stage.loadFile( msg.args[0].data, msg.kwargs );
                                }
                            }else{
                                stage_func.apply( stage, new_args );
                            }
                        }
                        break;
                    case 'Viewer': 
                        var viewer = this.stage.viewer;
                        var func = this.stage.viewer[msg.methodName];
                        func.apply( viewer, new_args );
                        break;
                    case 'compList':
                        var index = msg['component_index'];
                        var component = this.stage.compList[index];
                        var func = component[msg.methodName];
                        func.apply( component, new_args );
                        break;
                    case 'StructureComponent':
                        var component = this.structureComponent;
                        var func = component[msg.methodName];
                        func.apply( component, new_args );
                        break;
                    default: 
                        console.log( "nothing done for " + msg.target );
                        break;
                }
            }else if( msg.type == 'base64' ){
                console.log( "received base64 dict" );
                var base64Dict = msg.data;
                this.coordsDict = {};
                if ( "cache" in msg ){
                    this._cache = msg.cache;
                    this.model.set( "cache", this._cache );
                }
                for (var i = 0; i < Object.keys(base64Dict).length; i++) {
                     this.coordsDict[i] = this.mydecode( base64Dict[i]);
                }
            }else if( msg.type == 'get') {
                console.log( msg.data );

                if( msg.data == 'camera' ) {
                    this.send( JSON.stringify( this.stage.viewer.camera ) );
                }else if( msg.data == 'parameters' ){
                    this.send( JSON.stringify( this.stage.parameters ));
                }else{
                    console.log( "nothing done");
                }
            }
    },
    } );

    manager.WidgetManager.register_widget_view( 'NGLView', NGLView );

    return {
        'NGLView': NGLView,
    };

} );

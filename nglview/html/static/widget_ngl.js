
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
        "mdsrv": "../nbextensions/nglview/mdsrv"
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
                this.model.on( "change:representations", this.representationsChanged, this );

                // init structure loading
                this.model.on( "change:structure", this.structureChanged, this );

                // init setting of coordinates
                this.model.on( "change:coordinates", this.coordinatesChanged, this );

                // init setting of coordinates
                this.model.on( "change:coordinatesdict", this.coordinatesdictChanged, this );

                // init setting of frame
                this.model.on( "change:frame", this.frameChanged, this );

                // init parameters handling
                this.model.on( "change:parameters", this.parametersChanged, this );

                // init parameters handling
                this.model.on( "change:cache", this.cacheChanged, this );

                // get message from Python
                this.coordinatesdict = undefined;
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
            var representations = this.model.get( "representations" );
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

        coordinatesChanged: function(){
            var coordinates = this.model.get( "coordinates" );
            var component = this.structureComponent;
            if( coordinates && component ){
                var coords = new Float32Array( coordinates );
                component.structure.updatePosition( coords );
                component.updateRepresentations( { "position": true } );
            }

        },

        coordinatesdictChanged: function(){
            this.coordinatesdict = this.model.get( "coordinatesdict" );
        },

        frameChanged: function(){
            var frame = this.model.get( "frame" );
            if( this._cache ){
                var coordinates = this.coordinatesdict[frame];
                var component = this.structureComponent;
                if( coordinates && component ){
                    var coords = new Float32Array( coordinates );
                    component.structure.updatePosition( coords );
                    component.updateRepresentations( { "position": true } );
                }
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
            if ( msg.type == 'call_method' ){
               var new_args = msg.args.slice();
               new_args.push( msg.kwargs );

               if( msg.who == 'stage' ){
                   var stage_func = this.stage[msg.methodName];
                   var stage = this.stage;
                   if ( msg.methodName == 'screenshot' ){
                        NGL.screenshot( this.stage.viewer, msg.kwargs );
                   }else{
                       stage_func.apply( stage, new_args );
                   }
            }else if( msg.who == 'viewer' ){
                    var viewer = this.stage.viewer;
                    var viewer_func = this.stage.viewer[msg.methodName];
                    viewer_func.apply( viewer, new_args );
                }
            }else if( msg.type == 'coordinatesdict'){
                this.coordinatesdict = msg.data;
                console.log ( "received coordinatesdict" );
            }
        }

    } );

    manager.WidgetManager.register_widget_view( 'NGLView', NGLView );

    return {
        'NGLView': NGLView,
    };

} );

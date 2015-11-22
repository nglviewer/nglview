
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
        "NGL": "../nbextensions/nglview/ngl.min",
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


    var WIDTH = 400;
    var HEIGHT = 300;

    var NGLView = widget.DOMWidgetView.extend( {

        render: function(){

            // console.log( "NGLView", this )

            if( !Detector.webgl ) Detector.addGetWebGLMessage();

            NGL.init( function(){

                // init representations handling
                this.model.on( "change:representations", this.representationsChanged, this );

                // init structure loading
                this.model.on( "change:structure", this.structureChanged, this );

                // init setting of coordinates
                this.model.on( "change:coordinates", this.coordinatesChanged, this );

                // init NGL stage
                NGL.useWorker = false;
                this.stage = new NGL.Stage();
                this.stage.setTheme( "light" );
                this.structureComponent = undefined;
                this.$el.append( this.stage.viewer.container );
                this.setSize( WIDTH, HEIGHT );
                this.displayed.then( function(){
                    this.stage.handleResize();
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
                    .css( "top", "20px" )
                    .css( "left", "30px" )
                    .css( "background-color", "white" )
                    .css( "padding", "2px 5px 2px 5px" )
                    .css( "opacity", "0.7" )
                    .appendTo( this.stage.viewer.container );
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
                    this.$pickingInfo.prop( "innerText", pickingText );
                }, this );

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
                this.stage.loadFile(
                    new Blob( [ structure.data ], { type: "text/plain" } ),
                    { ext: structure.ext, defaultRepresentation: false }
                ).then( function( component ){
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

        setSize: function( width, height ){
            this.stage.viewer.container.style.width = width + "px";
            this.stage.viewer.container.style.height = height + "px";
            this.stage.handleResize();
        }

    } );

    manager.WidgetManager.register_widget_view( 'NGLView', NGLView );

    return {
        'NGLView': NGLView,
    };

} );

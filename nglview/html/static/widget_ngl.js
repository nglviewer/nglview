
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
    "THREE", "Detector", "async", "Promise", "sprintf", "JSZip", "pako",
    "LZMA", "bzip2", "chroma", "jsfeat", "signals", "NGL", "mdsrv"
], function(
    widget, manager, _THREE, _Detector, async, _Promise, _sprintf, _JSZip, pako,
    _LZMA, _bzip2, chroma, _jsfeat, signals, _NGL, _NGL_mdsrv
){

    window.async = async;
    window.pako = pako;
    window.signals = signals;
    window.chroma = chroma;

    console.log(NGL)
    console.log(_NGL_mdsrv)

    var WIDTH = 400;
    var HEIGHT = 300;

    var NGLView = widget.DOMWidgetView.extend( {

        render: function(){

            console.log( "NGLView", this )

            if( !Detector.webgl ) Detector.addGetWebGLMessage();

            NGL.init( function(){

                // // init selection input
                // this.$selection = $( "<input />" );
                // this.$selection.val( this.model.get( "selection" ) );
                // this.$selection.change( function(){
                //     this.model.set( "selection", this.$selection.val() );
                //     this.model.save();
                // }.bind( this ) );
                // this.$el.append( this.$selection );
                // this.model.on( "change:selection", this.selectionChanged, this );

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

            }.bind( this ) );

        },

        // selectionChanged: function(){
        //     var selection = this.model.get( "selection" );
        //     this.$selection.val( selection );
        //     this.stage.getRepresentationsByName().setSelection( selection );
        // },

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
                    console.log( "MOIN", component );
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

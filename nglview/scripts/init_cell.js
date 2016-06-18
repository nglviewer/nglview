define([
    'jquery',
    'base/js/namespace',
    'base/js/events'
], function (
    $,
    IPython,
    events
) {
    var run_init_cells = function(){
        var cell = IPython.notebook.get_cell(0);
        cell.execute();
    };

    var load_ipython_extension = function() {
        events.on('kernel_ready.Kernel', run_init_cells);
    };

    return {
        load_ipython_extension : load_ipython_extension
    };
});

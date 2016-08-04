from functools import partial

__alll__ = ['js_clean_error_output', 'js_launch_qtconsole',
            'js_clean_empty_output_area', 'js_open_url_template',
            '_set_ipython_cell', 'ngl_demo', 'init_js_funcs']

def run(command):
    from IPython.display import display, Javascript
    display(Javascript(command))

js_clean_empty_output_area = partial(run, command="""
var output_area = $(".output_area");

for (var i=0; i < output_area.length; i++){
    if (!output_area[i].innerText){
        output_area[i].remove();
    }
}
""")

js_launch_qtconsole = partial(run, command="""
Jupyter.notebook.kernel.execute('%qtconsole')
""")

js_open_url_template = """
window.open({url});
"""

js_clean_error_output = partial(run, command="""
var cells = Jupyter.notebook.get_cells();

for (var i = 0; i < cells.length; i++){
    var cell = cells[i];
    if (cell.output_area.outputs.length > 0) {
        var out = cell.output_area.outputs[0];
        if (out.output_type == 'error') {
            cell.clear_output();
        }
    }
}
""")

def _set_ipython_cell(background='#67a9cf'):
    cm = """
    var setIPythonLikeCell = function(){
        //var cell = Jupyter.notebook.insert_cell_at_bottom();
        var cell = Jupyter.notebook.get_selected_cell();
        var $el = cell.element;
        cell.set_text('');
        $el.css({'background': '%s'});

        var handler = function(event) {
            var selected_cell = Jupyter.notebook.get_selected_cell();
            if (selected_cell.cell_id === cell.cell_id){
                selected_cell.execute();
                selected_cell.set_text('');
            }
            return false;
        };

        var action = {
            help: 'run cell',
            help_index: 'zz',
            handler: handler
        };

        Jupyter.keyboard_manager.edit_shortcuts.add_shortcut('enter', action); 
    };

    setIPythonLikeCell()
    """ % background

    cm2 = """
    var cell = Jupyter.notebook.get_selected_cell();
    cell.element.draggable().resizable();
    """
    from IPython.display import display, Javascript
    display(Javascript(cm))
    display(Javascript(cm2))

def ngl_demo(width=400, height=400):
    """make a viewport, create stage object and populate NGL namespace
    """
    from IPython.display import display, Javascript, HTML

    command = """
    <div id='viewport'></div>
    <script>
        $('#viewport').width(%s).height(%s);
    </script>
    """ % (width, height)

    command2 = """
    <script>
        var NGL = require('nbextensions/nglview/index').NGL;
        var stage = new NGL.Stage('viewport')
    </script>
    """

    display(HTML(command))
    display(HTML(command2))

def init_js_funcs():
    """print
    """
    from IPython.display import display, Javascript, HTML

    command = """
    <script>
        var print = function(x){
            for (var i in x){
                console.log(i)
            }
        }
    </script>
    """
    display(HTML(command))

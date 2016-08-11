from functools import partial
from IPython.display import display, Javascript

__alll__ = ['clean_error_output', 'launch_qtconsole',
            'clean_empty_output_area', 'open_url_template',
            '_set_ipython_cell', 'ngl_demo', 'init_funcs',
            '_move_notebook_to_the_right', '_move_notebook_to_the_left',
            '_reset_notebook',
            '_set_notebook_draggable']

def run(command):
    display(Javascript(command))

def _set_notebook_width(width='20%'):
    script_template = """
    var cb = Jupyter.notebook.container;

    cb.width('{width}');
    cb.offset({{'left': 0}})
    """
    display(Javascript(script_template.format(width=width)))

def _set_notebook_draggable(yes=True):
    script_template = """
    var x = $('#notebook-container');
    x.draggable({args});
    """
    if yes:
        display(Javascript(script_template.replace('{args}', '')))
    else:
        display(Javascript(script_template.replace('{args}', '"destroy"')))

def _move_notebook_to_the_right():
    script_template = """
    var x = $('#notebook-container');
    x.css({position: "relative", left: "20%"});
    """
    display(Javascript(script_template))

def _move_notebook_to_the_left():
    script_template = """
    var cb = Jupyter.notebook.container;

    cb.offset({'left': 0})
    """
    display(Javascript(script_template))

def _reset_notebook():
    script_template = """
    var x = $('#notebook-container');
    x.width('30%');
    x.css({position: "relative", left: "0%"});
    """
    display(Javascript(script_template))


clean_empty_output_area = partial(run, command="""
var output_area = $(".output_area");

for (var i=0; i < output_area.length; i++){
    if (!output_area[i].innerText){
        output_area[i].remove();
    }
}
""")

launch_qtconsole = partial(run, command="""
Jupyter.notebook.kernel.execute('%qtconsole')
""")

open_url_template = """
window.open({url});
"""

clean_error_output = partial(run, command="""
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

def init_funcs():
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

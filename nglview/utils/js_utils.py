from functools import partial

from IPython.display import HTML, Javascript, display

__alll__ = [
    'clean_error_output',
    'launch_qtconsole',
    'clean_empty_output_area',
    'open_url_template',
    '_set_ipython_cell',
    'ngl_demo',
    'init_funcs',
]


def run(command):
    display(Javascript(command))


# cell.clear_output()
def hide_toolbar():
    run("$('#maintoolbar').hide()")
    run("$('#header-container').hide()")


def show_toolbar():
    run("$('#maintoolbar').show()")
    run("$('#header-container').show()")


def execute(command):
    run(f'Jupyter.notebook.kernel.execute("{command}")')


def _set_notebook_width(width='20%', left_padding=0):
    script_template = """
    var cb = Jupyter.notebook.container;

    cb.width('{width}');

    """.format(width=width)

    if left_padding is not None:
        offset_str = f"cb.offset({{'left': {left_padding}}})"
    else:
        offset_str = ''
    command = script_template + offset_str
    run(command)


clean_empty_output_area = partial(run,
                                  command="""
var output_area = $(".output_area");

for (var i=0; i < output_area.length; i++){
    if (!output_area[i].innerText){
        output_area[i].remove();
    }
}
""")

launch_qtconsole = partial(run,
                           command="""
Jupyter.notebook.kernel.execute('%qtconsole')
""")

open_url_template = """
window.open({url});
"""

clean_error_output = partial(run,
                             command="""
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
    run(cm)
    run(cm2)


def ngl_demo(width=400, height=400):
    """make a viewport, create stage object and populate NGL namespace
    """
    command = """
    <div id='viewport'></div>
    <script>
        $('#viewport').width(%s).height(%s);
    </script>
    """ % (width, height)

    command2 = """
    <script>
        var NGL = require('nbextensions/nglview-js-widgets/index').NGL;
        var stage = new NGL.Stage('viewport')
    </script>
    """

    display(HTML(command))
    display(HTML(command2))


def init_funcs():
    """print
    """
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

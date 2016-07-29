from functools import partial

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

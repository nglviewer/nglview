js_clean_empty_output_area = """
var output_area = $(".output_area");

for (var i=0; i < output_area.length; i++){
    if (!output_area[i].innerText){
        output_area[i].remove();
    }
}
"""

js_launch_qtconsole = """
IPython.notebook.kernel.execute('%qtconsole')
"""

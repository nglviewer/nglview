"""highly experimental. Use with your own risk

    >>> from nglview import theme
    >>> theme.oceans16() # oceans16 https://github.com/dunovank/jupyter-themes

Retart your notebook to reset to default Jupyter theme

If you want to set global theme for your notebook, it's better to install jupyter-themes
https://github.com/dunovank/jupyter-themes

"""
import os

style = """
<style id='nglview_style'>
{}
</style>

<script>
$('#nglview_style').appendTo('head');
</script>
"""


def _get_theme(css_file):
    from IPython.display import HTML
    return HTML(_get_css_content(css_file))


def _get_css_content(css_file):
    dirname = os.path.dirname(os.path.abspath(__file__))
    css_file = os.path.join(dirname, css_file)
    css = open(css_file).read()
    return style.format(css)


def oceans16():
    return _get_theme('oceans16.css')


def reset(hide_toolbar=False):
    from IPython.display import Javascript, display
    from nglview import js_utils
    display(Javascript('$("#nglview_style").remove()'))
    js_utils.clean_empty_output_area()

    if hide_toolbar:
        js_utils.hide_toolbar()
    else:
        js_utils.show_toolbar()

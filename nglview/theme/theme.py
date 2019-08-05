import os
from IPython.display import display
from ipywidgets import HTML
from pathlib import Path

from ..base import _singleton


def _get_css_content(css_file):
    p = Path(__file__).resolve().parent / css_file
    return p.read_text()
    
    
class Theme:
    def __init__(self, theme=None):
        self._html = HTML()
        display(self._html)
        if theme == 'light':
            self.light()
        elif theme == 'dark':
            self.dark()
        else:
            raise ValueError("Unsupported theme")


    def oceans16(self):
        self._html.value = (self._html.value + '\n' +
                            '<style>\n' +
                            _get_css_content('oceans16.css') +
                            '</style>')
    
    
    def remove(self):
        self._html.value = ''
    
    
    def dark(self):
        self._html.value = ('<style>\n' + 
                      _get_css_content('dark.css') +
                      _get_css_content('main.css') +
                      '</style>')
    
    def light(self):
        self._html.value = ('<style>\n' + 
                      _get_css_content('light.css') +
                      _get_css_content('main.css') +
                      '</style>')

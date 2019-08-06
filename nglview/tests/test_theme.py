import nglview as nv
from nglview.theme import ThemeManager


def test_theme():
    m  = ThemeManager()
    assert m is ThemeManager() # singleton
    assert m._theme_css == ''
    m.dark()
    assert m._theme == 'dark'
    assert m._theme_css is not None
    m.light()
    assert m._theme == 'light'
    m.remove()
    assert m._theme is None
    assert m._theme_css == ''

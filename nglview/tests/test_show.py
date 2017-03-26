from mock import MagicMock, patch
import nglview
# TODO : add more show_xxx

@patch('schrodinger.structure.Structure')
def test_show_schrodinger_structure(MockStructure):
    s = MockStructure()
    nglview.show_schrodinger_structure(structure)

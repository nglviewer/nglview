import nglview as nv
from nglview.color import ColormakerRegistry as cm
from nglview.color import _ColormakerRegistry
from mock import patch, call, ANY


def test_color_scheme():
    assert _ColormakerRegistry().model_id == _ColormakerRegistry(
    ).model_id  # only one 1 model

    with patch.object(cm, '_call') as mock_call:
        cm.add_scheme('abc', [['blue', '1-10']])
        assert mock_call.call_args_list == [
            call('addSelectionScheme', 'abc', [['blue', '1-10']])
        ]
        mock_call.reset_mock()
        cm.add_scheme(
            'abc', """
              this.colorAtom = function(atom){
              }
        """)
        mock_call.assert_called_once_with("executeCode", ANY)

    # handle before and after display
    cm.add_scheme('abc', [['blue', '1-10']])
    assert cm._msg_ar == cm._msg_q == [{
        'args': ('abc', [['blue', '1-10']]),
        'kwargs': {},
        'methodName': 'addSelectionScheme',
        'type': 'callMethod'
    }]
    cm._ready = True
    assert cm._msg_q == []
    cm._ipython_display_()

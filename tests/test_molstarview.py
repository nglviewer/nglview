import base64
from unittest.mock import patch, MagicMock

import pytest
from nglview.widget_molstar import MolstarView


@pytest.fixture
def molstar_view():
    return MolstarView()


def test_initialization(molstar_view):
    assert isinstance(molstar_view, MolstarView)
    assert molstar_view._view_name == 'MolstarView'
    assert molstar_view._model_name == 'MolstarModel'
    assert molstar_view.frame == 0
    assert not molstar_view.loaded
    assert molstar_view.molstate == {}
    assert molstar_view._state is None


@patch('nglview.widget_molstar.widgets.Widget')
def test_handle_nglview_custom_message_export_image(mock_widget, molstar_view):
    mock_image = MagicMock()
    mock_widget.widgets = {'test_model_id': mock_image}
    msg = {
        "type": "exportImage",
        "model_id": "test_model_id",
        "data": base64.b64encode(b'test_data').decode('utf-8')
    }
    molstar_view._handle_nglview_custom_message(None, msg, None)
    assert mock_image.value == base64.b64decode(msg["data"])


def test_handle_nglview_custom_message_state(molstar_view):
    msg = {"type": "state", "data": "test_state"}
    molstar_view._handle_nglview_custom_message(None, msg, None)
    assert molstar_view._state == "test_state"


def test_handle_nglview_custom_message_request_loaded(molstar_view):
    msg = {"type": "request_loaded", "data": True}
    molstar_view._handle_nglview_custom_message(None, msg, None)
    assert molstar_view.loaded


def test_handle_nglview_custom_message_get_camera(molstar_view):
    msg = {"type": "getCamera", "data": "test_camera"}
    molstar_view._handle_nglview_custom_message(None, msg, None)
    assert molstar_view._molcamera == "test_camera"


@patch('nglview.widget_molstar.MolstarView._fire_callbacks')
def test_on_loaded(mock_fire_callbacks, molstar_view):
    molstar_view._callbacks_before_loaded = ['callback1', 'callback2']
    molstar_view.on_loaded({'new': True})
    mock_fire_callbacks.assert_called_once_with(['callback1', 'callback2'])


@patch('nglview.widget_molstar.MolstarView._remote_call')
def test_load_structure_data(mock_remote_call, molstar_view):
    molstar_view._load_structure_data('test_data', 'pdb', 'default')
    mock_remote_call.assert_called_once_with(
        "loadStructureFromData",
        target="Widget",
        args=['test_data', 'pdb', 'default'])


@patch('nglview.widget_molstar.MolstarView._update_max_frame')
@patch('nglview.widget_molstar.MolstarView._load_structure_data')
def test_add_trajectory(mock_load_structure_data, mock_update_max_frame,
                        molstar_view):
    mock_trajectory = MagicMock()
    mock_trajectory.get_structure_string.return_value = 'test_structure_string'
    mock_trajectory.id = 'test_id'
    molstar_view.add_trajectory(mock_trajectory)
    mock_load_structure_data.assert_called_once_with('test_structure_string',
                                                     'pdb')
    mock_update_max_frame.assert_called_once()
    assert mock_trajectory in molstar_view._trajlist
    assert 'test_id' in molstar_view._molstar_component_ids


@patch('nglview.widget_molstar.MolstarView._load_structure_data')
def test_add_structure(mock_load_structure_data, molstar_view):
    mock_structure = MagicMock()
    mock_structure.get_structure_string.return_value = 'test_structure_string'
    mock_structure.id = 'test_id'
    molstar_view.add_structure(mock_structure)
    mock_load_structure_data.assert_called_once_with('test_structure_string',
                                                     'pdb')
    assert 'test_id' in molstar_view._molstar_component_ids


# FIXME: failed tests. Why?

# @patch('nglview.widget_molstar.MolstarView._set_coordinates'
# def test_on_frame_changed(mock_set_coordinates, molstar_view):
#     molstar_view.frame = 10
#     molstar_view._on_frame_changed({'new': 10})
#     mock_set_coordinates.assert_called_once_with(10)


# @patch('nglview.widget_molstar.MolstarView._remote_call')
# def test_add_representation(mock_remote_call, molstar_view):
#     params = {'component': 'test_component', 'opacity': 0.4}
#     molstar_view.add_representation(**params)
#     expected_params = {'component': 'test_component', 'opacity': 0.4}
#     mock_remote_call.assert_called_once_with('addRepresentation',
#                                              args=[expected_params, 0])

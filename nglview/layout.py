# TODO: reorg
# simplify code
from ipywidgets import Box, Button, Label, Layout, Tab, ToggleButton


def make_form_item_layout():
    # protect nglview if it is not used in notebook (wheen there is not `comm`
    return Layout(
        display='flex',
        flex_flow='row',
        justify_content='space-between',
    )


def _make_box_layout(width='100%'):
    return Layout(display='flex',
                  flex_flow='column',
                  align_items='stretch',
                  width=width)


def _relayout(box, form_item_layout):
    form_items = []
    for kid in box.children:
        if hasattr(kid,
                   'description') and not isinstance(kid,
                                                     (Button, ToggleButton)):
            label_value = kid.description
            kid.description = ''
        else:
            label_value = ''
        if isinstance(kid, Button):
            box2 = Box([
                kid,
            ], layout=form_item_layout)
        else:
            box2 = Box([Label(value=label_value), kid],
                       layout=form_item_layout)
        form_items.append(box2)

    return form_items


def _relayout_master(box, width='20%'):
    """make nicer layout for box.

    This method will take the `description` of each child to make corresponding Label
    The `description` will be cleared.
    """
    old_children = box.children[:]
    form_items = _relayout(box, make_form_item_layout())
    form = Box(form_items, layout=_make_box_layout(width=width))
    form._ngl_children = old_children
    return form


def _make_autofit(box):
    '''

    Parameters
    ----------
    box : ipywidgets.Box
        children is a list of buttons

    Returns
    -------
    relayouted box
    '''
    items_layout = Layout(flex='1 1 auto', width='auto')

    box.layout = items_layout
    return box


def _make_delay_tab(box_factory, selected_index=0):
    """

    Parameters
    ----------
    box_factory : list of (func, tab_name)

    Example of box_factory: [(_make_gen_box, 'General'),
                             (_make_repr_box, 'Representation')]
    """

    tab = Tab([Box() for box, _ in box_factory])
    [tab.set_title(i, title) for i, (_, title) in enumerate(box_factory)]

    # trick
    if not tab.children[selected_index].children:
        tab.selected_index = 1

    def on_update_selected_index(change):
        index = change['new']
        if not tab.children[index].children:
            # make widget on demand
            tab.children[index].children = [
                box_factory[index][0](),
            ]

    tab.observe(on_update_selected_index, names='selected_index')

    # trigger
    tab.selected_index = selected_index

    return tab

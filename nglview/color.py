COLOR_SCHEMES = [
    " ", "picking", "random", "uniform", "atomindex", "residueindex",
    "chainindex", "modelindex", "sstruc", "element", "resname", "bfactor",
    "hydrophobicity", "value", "volume", "occupancy"
]

_USER_COLOR_DICT = {}


class _ColorScheme:
    _color_dict = {}

    def __init__(self, args, label):
        # FIXME: validate `args`
        self._color_scheme = args
        self._label = f'user_{label}'
        _USER_COLOR_DICT[self._label] = self._color_scheme

    @property
    def data(self):
        return {'data': self._color_scheme, 'label': self._label}

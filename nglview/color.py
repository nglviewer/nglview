COLOR_SCHEMES = [
    " ", "picking", "random", "uniform", "atomindex", "residueindex",
    "chainindex", "modelindex", "sstruc", "element", "resname", "bfactor",
    "hydrophobicity", "value", "volume", "occupancy"
]


class _ColorScheme:
    def __init__(self, args):
        # FIXME: validate `args`
        self._color_scheme = args

    @property
    def data(self):
        return self._color_scheme

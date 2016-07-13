REPR_NAME_PAIRS = [
        ('point', 'point'),
        ('line', 'line'),
        ('rope', 'rope'),
        ('tube', 'tube'),
        ('trace', 'trace'),
        ('label', 'label'),
        ('unitcell', 'unitcell'),
        ('cartoon', 'cartoon'),
        ('licorice', 'licorice'),
        ('distance', 'distance'),
        ('ribbon', 'ribbon'),
        ('surface', 'surface'),
        ('backbone', 'backbone'),
        ('contact', 'contact'),
        ('hyperball', 'hyperball'),
        ('rocket', 'rocket'),
        ('helixorient', 'helixorient'),
        ('simplified_base', 'base'),
        ('ball_and_stick', 'ball+stick'),
        ]

REPR_NAMES = list(sorted(set([name for pairs in REPR_NAME_PAIRS for name in pairs])))

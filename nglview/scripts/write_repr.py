#!/usr/bin/env python

reprs = [
        ('point', 'point'),
        ('line', 'line'),
        ('cartoon', 'cartoon'),
        ('licorice', 'licorice'),
        ('ribbon', 'ribbon'),
        ('surface', 'surface'),
        ('trace', 'trace'),
        ('tube', 'tube'),
        ('label', 'label'),
        ('backbone', 'backbone'),
        ('ball_and_stick', 'ball+stick'),
        ('contact', 'contact'),
        ('crossing', 'crossing'),
        ('helixorient', 'helixorient'),
        ('hyperball', 'hyperball'),
        ('rocket', 'rocket'),
        ('rope', 'rope'),
        ('simplified_base', 'base'),
        ]

for repr in reprs:
    print("")
    print("    def add_{method}(self, selection='all', **kwd):".format(method=repr[0]))
    print("        self.add_representation('{representation}', selection, **kwd)".format(representation=repr[1]))

#!/usr/bin/env python

reprs = [
        ('cartoon', 'cartoon'),
        ('line', 'line'),
        ('label', 'label'),
        ]

for repr in reprs:
    print("")
    print("    def {method}(self, selection, **kwd):".format(method=repr[0]))
    print("        self.add_representation('{representation}', selection, **kwd)".format(representation=repr[1]))

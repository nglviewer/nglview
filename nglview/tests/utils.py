import os

def get_fn(fn):
    this_path = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(this_path, 'data', fn)

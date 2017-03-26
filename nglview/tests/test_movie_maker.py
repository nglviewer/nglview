import nglview
import pytraj

# local
from make_dummy_comm import *

def test_movie_maker():
    from nglview.contrib.movie import MovieMaker
    traj = pytraj.datafiles.load_tz2()
    view = nglview.show_pytraj(traj)

    movie = MovieMaker(view)
    movie.make()

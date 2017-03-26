import nglview
import pytraj

# local
from make_dummy_comm import *

def test_movie_maker():
    from nglview.contrib.movie import MovieMaker
    traj = pytraj.datafiles.load_tz2()
    view = nglview.show_pytraj(traj)

    movie = MovieMaker(view, in_memory=False)
    movie.make()

    movie = MovieMaker(view, in_memory=True)
    movie.make()

    movie.interupt()

    movie = MovieMaker(view, in_memory=True, skip_render=True)
    movie.make()

    movie = MovieMaker(view, in_memory=True, skip_render=False,
            start=0, stop=1, step=1)
    movie.make()

    print(movie.thread)

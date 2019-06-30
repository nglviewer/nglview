import os

from mock import MagicMock, patch

import nglview
import pytraj
from make_dummy_comm import *

# local


class FakeEvent:
    def is_set(self):
        return self._event_set


@patch('moviepy.editor.ImageSequenceClip')
def test_movie_maker(ImageSequenceClip):
    from nglview.contrib.movie import MovieMaker
    ImageSequenceClip.write_gif = MagicMock()
    ImageSequenceClip.write_videofile = MagicMock()
    traj = pytraj.datafiles.load_tz2()
    view = nglview.show_pytraj(traj)

    movie = MovieMaker(view, in_memory=False)
    movie.download_folder = os.path.join(os.path.dirname(__file__), 'data')

    # fake _event
    movie._event = FakeEvent()
    movie._event._event_set = True
    movie.make()

    movie._event._event_set = False
    movie.make()

    movie = MovieMaker(view, in_memory=False)
    movie.skip_render = True
    movie.make()

    movie.output = 'hello.mp4'
    movie.make()

    movie.interupt()
    movie._event = None
    movie.interupt()

    movie.in_memory = True
    movie.make()

    movie._event._event_set = False
    movie.make()

    movie = MovieMaker(view,
                       download_folder='here',
                       render_params=dict(factor=4),
                       moviepy_params={},
                       stop=2)
    movie.make()


def test_movie_maker_base64_to_ndarray():
    from nglview.contrib.movie import MovieMaker
    s = 'iVBORw0KGgoAAAANSUhEUgAAAAUAAAAFCAYAAACNbyblAAAAHElEQVQI12P4//8/w38GIAXDIBKE0DHxgljNBAAO9TXL0Y4OHwAAAABJRU5ErkJggg=='
    MovieMaker._base64_to_ndarray(s)

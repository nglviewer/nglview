import time

from ..parameters import REPRESENTATION_NAME_PAIRS

import cv2

def wait(widget, attribute='value', timeout=5):
    """ EXPERIMENTAL. Require `ipython_blocking` package.

    Block further code execution until `attribute` of `widget` is updated.
    """
    from ipython_blocking import CaptureExecution
    c = CaptureExecution()
    with c:
        t0 = time.time()
        while True or time.time() - t0 > timeout:
            attr = getattr(widget, attribute)
            if attr:
                break
            c.step()
    return widget


def get_widget_by_name(box, widget_name):

    if hasattr(box, '_ngl_children'):
        children = box._ngl_children
    elif hasattr(box, 'children'):
        children = box.children
    else:
        children = None

    if children is not None:
        for widget in children:
            if hasattr(widget,
                       '_ngl_name') and widget._ngl_name == widget_name:
                return widget
    return None


def _add_repr_method_shortcut(self, other):
    from types import MethodType

    def make_func_add(rep):
        """return a new function object
        """

        def func(this, selection='all', **kwargs):
            """
            """
            self.add_representation(repr_type=rep[1],
                                    selection=selection,
                                    **kwargs)

        func.__doc__ = """Shortcut for `add_representation` method

        Examples
        --------
        >>> view.add_{name}()
        >>> # is equal to
        >>> view.add_representation('{name}')
        """.format(name=rep[0])
        return func

    def make_func_remove(rep):
        """return a new function object
        """

        def func(this, **kwargs):
            """
            """
            self._remove_representations_by_name(repr_name=rep[1], **kwargs)

        return func

    def make_func_update(rep):
        """return a new function object
        """

        def func(this, **kwargs):
            """
            """
            self._update_representations_by_name(repr_name=rep[1], **kwargs)

        return func

    for rep in REPRESENTATION_NAME_PAIRS:
        for make_func, root_fn in [(make_func_add, 'add'),
                                   (make_func_update, 'update'),
                                   (make_func_remove, 'remove')]:
            func = make_func(rep)
            fn = '_'.join((root_fn, rep[0]))
            setattr(self, fn, MethodType(func, other))

def compare_two_images(image_path_1,
                       image_path_2,
                       standard_width=200,
                       standard_height=200,
                       convert_to_gray=False,
                       threshold=0.97):
    """
    Comparing two images. Two images will be convert to new size
    (standard_width, standard_height) and then compute structural similarity
    to compare them.
    Input:
        image_path_1: string, path for the first image
        image_path_2: string, path for the second image
        standard_width: float (number of pixels) for a standard width
        standard_height: float (number of pixels) for a standard height
        convert_to_gray: boolean; if user want to convert a gray scale after
            resize the images.
        threshold: float, threshold for the structural similarity index
            between two images. The default value is 0.97
    Output:
        Boolean value
    """
    from skimage.measure import compare_ssim as ssim

    # load image and resize
    img_1 = resize_and_load_image(image_path=image_path_1,
                                 new_width=standard_width,
                                 new_height=standard_height,
                                 convert_to_gray=convert_to_gray)
    img_2 = resize_and_load_image(image_path=image_path_2,
                                 new_width=standard_width,
                                 new_height=standard_height,
                                 convert_to_gray=convert_to_gray)
    b_multichanel = not convert_to_gray
    s = ssim(img_1, img_2, multichannel=b_multichanel)

    return s > threshold

def resize_and_load_image(image_path,
                          new_width=200,
                          new_height=200,
                          convert_to_gray=True):
    """
    Load image and resize it to new_width and new_height
    Input:
        image_path: string, path for the image
        new_width, new_height: float (number of pixels); the new width and height
        convert_to_gray: boolean: and option to convert it to gray
    Output:
        pixel values as a nd np.array
    """
    img = cv2.imread(image_path)
    # resize
    img_resize = cv2.resize(img, (new_width, new_height),interpolation=cv2.INTER_AREA)
    # convert images to gray scale
    if convert_to_gray:
        img_resize = cv2.cvtColor(img_resize, cv2.COLOR_BGR2GRAY)

    return img_resize


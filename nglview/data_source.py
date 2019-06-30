from .utils import js_utils


class DatasourceRegistry:
    sources = {}

    @classmethod
    def add(cls, name, url):
        """

        Examples
        --------
        >>> import nglview as nv
        >>> from nglview.data_source import DatasourceRegistry # doctest: +SKIP
        ... DatasourceRegistry.add("data", "//cdn.rawgit.com/arose/ngl/v2.0.0-dev.32/data/")
        ... view = nv.NGLWidget()
        ... v.add_component('data://1CRN.cif') # 1CRN.cif is from above url
        ... v
        """
        cls.sources[name] = url

        js_utils.run(f"""
        var NGL = require("nglview-js-widgets").NGL;

        NGL.DatasourceRegistry.add("{name}",
            new NGL.StaticDatasource("{url}"))
        """)

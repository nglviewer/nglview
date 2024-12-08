from pydantic import BaseModel
from typing import Any, List, Optional

class IMessage(BaseModel):
    type: str
    data: Optional[Any] = None
    ID: Optional[str] = None
    render_params: Optional[Any] = None
    methodName: Optional[str] = None
    target: Optional[str] = None
    args: Optional[List[Any]] = None
    kwargs: Optional[dict] = None
    component_index: Optional[int] = None
    repr_index: Optional[int] = None
    last_child: Optional[bool] = None
    reconstruc_color_scheme: Optional[bool] = None
    movie_making: Optional[bool] = None
    buffers: Optional[List[Any]] = None
    content: Optional[Any] = None
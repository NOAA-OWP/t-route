"""Author: Tadd Bindas"""

from pydantic import BaseModel


class TRouteOuput(BaseModel):
    """
    A schema to define successful t-route output

    Attributes
    ----------
    message: str
        The 
    lid : str
        The location ID belonging to the point being routed
    feature_id : str
        The COMID, or hf_id, belonging to the point being routed
    """

    message: str
    lid: str
    feature_id: str

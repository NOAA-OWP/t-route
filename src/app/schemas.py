from datetime import datetime
from typing import List, Optional, Union

from pydantic import BaseModel, ConfigDict

class TRouteOuput(BaseModel):
    """
    River Forecast Center information.

    Attributes
    ----------
    abbreviation : str
        The abbreviated name of the RFC.
    name : str
        The full name of the RFC.
    """
    message: str
    lid: str
    feature_id: str
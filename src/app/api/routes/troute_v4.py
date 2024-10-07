"""Author: Tadd Bindas"""
from datetime import datetime
from typing import Annotated

from fastapi import APIRouter, Depends
from fastapi.responses import JSONResponse
from nwm_routing import main_v04 as t_route
from pydantic import conint

from app.api.services.initialization import (create_initial_start_file,
                                             create_params, edit_yaml)
from app.core.settings import Settings
from app.core import get_settings
from app.schemas import TRouteStatus

router = APIRouter()


@router.get("/", response_model=TRouteStatus)
async def get_gauge_data(
    lid: str,
    feature_id: str,
    hy_id: str,
    initial_start: float,
    start_time: datetime,
    num_forecast_days: conint(ge=1, le=30),
    settings: Annotated[Settings, Depends(get_settings)],
) -> TRouteStatus:
    """An API call for running T-Route within the context of replace and route

    Parameters:
    ----------
    lid: str
        The Location of the RFC Point
    feature_id: str
        The COMID associated with the LID
    hy_id: str
        The HY_ID associated with the LID
    initial_start: float
        The initial start for T-Route
    start_time: str
        The start time for the forecast
    num_forecast_days: int
        The number of days in the forecast

    Returns:
    -------
    TRouteOutput
        A successful T-Route run
    """
    base_config = settings.base_config
    params = create_params(
        lid, feature_id, hy_id, initial_start, start_time, num_forecast_days, settings
    )
    restart_file = create_initial_start_file(params, settings)
    yaml_file_path = edit_yaml(base_config, params, restart_file)
    try:
        t_route(["-f", yaml_file_path.__str__()])
    except Exception as e:
        JSONResponse(
            status_code=500,
            content={"message": e},
        )

    yaml_file_path.unlink()

    return TRouteStatus(
        message="T-Route run successfully",
        lid=lid,
        feature_id=feature_id,
    )

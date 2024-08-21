from pathlib import Path

from pydantic import BaseSettings

class Settings(BaseSettings):
    """
    Configuration settings for the application.

    This class uses Pydantic's BaseSettings to manage configuration,
    allowing for easy integration with environment variables and
    configuration files.

    Parameters
    ----------
    **data : dict
        Additional keyword arguments to be passed to the parent class.

    Attributes
    ----------
    api_v1_str : str
        The base API string.
    project_name : str
        The project's name for the OPEN API spec

    Notes
    -----
    The configuration is initially read from a 'config.ini' file and can be
    overridden by environment variables.
    """
    api_v1_str: str = "/api/v1"
    base_config: Path = "/t-route/src/app/core/base_config.yaml"
    project_name: str = "T-Route"
    qlat_input_path: str = "/t-route/data/rfc_channel_forcings/{}/"
    geofile_path: str = "/t-route/data/rfc_geopackage_data/{}/subset.gpkg"

"""Author: Tadd Bindas"""

from pathlib import Path

from pydantic import BaseSettings


class Settings(BaseSettings):
    """
    Configuration settings for the application.

    This class uses Pydantic's BaseSettings to manage configuration,
    allowing for easy integration with environment variables and
    configuration files.

    Attributes
    ----------
    api_v1_str : str
        The base API string.
    base_config: str
        The path to the base_config file that we build custom config files from
    project_name : str
        The project's name for the OPENAPI spec
    qlat_input_path: str
        The path to the docker folders where rfc_channel_forcings live in the shared volume
    restart_path: str
        The path to where restart files will be saved to in the shared volume
    restart_file: str
        The regex string for finding restart files
    geofile_path: str
        The path to the docker folders where the geopackage is located in the shared volume

    Notes
    -----
    The configuration is initially read from a 'config.ini' file and can be
    overridden by environment variables.
    """

    api_v1_str: str = "/api/v1"
    base_config: Path = "/t-route/src/app/core/base_config.yaml"
    project_name: str = "T-Route"
    qlat_input_path: str = "/t-route/data/rfc_channel_forcings/{}/"
    restart_path: str = "/t-route/data/troute_restart/{}/"
    restart_file: str = "HYDRO_RST_{}_DOMAIN1"
    geofile_path: str = "/t-route/data/rfc_geopackage_data/{}/downstream.gpkg"

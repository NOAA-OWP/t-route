"""
=================================================
USGS NWIS Instantaneous Values REST Client
=================================================
This module provides an IVDataService class that provides a convenient interface to
the USGS NWIS Instantaneous Values (IV) REST Service API.
Classes
-------
   IVDataService

"""

import datetime
from multiprocessing import Pool
from functools import partial
from typing import List, Tuple, Union, Iterable
import re
import requests

import numpy as np
import pandas as pd
#from evaluation_tools._restclient import RestClient
from _restclient import RestClient


class IVDataService:
    """A REST client class.
    The IVDataService class provides various methods for constructing
    requests, retrieving data, and parsing responses from the NWIS IV
    Service.

    Parameters
    ----------
    processes : int
        Max multiprocessing processes, default 3
    retry : int
        Max number of, per request retries, default 3
    """

    # Class level variables
    _datetime_format = "%Y-%m-%dT%H:%M%z"
    _base_url = "https://waterservices.usgs.gov/nwis/iv/"
    _requests_cache_filename = "nwisiv_cache"
    _headers = {"Accept-Encoding": "gzip, compress"}

    def __init__(self, processes: int = 3, retries: int = 3):
        self._procs = processes

        self._restclient = RestClient(
            base_url=self._base_url,
            headers=self._headers,
            requests_cache_filename=self._requests_cache_filename,
            retries=retries,
        )

    @classmethod
    def get(
        cls,
        sites: Union[str, List[str]],
        parameterCd: str = "00060",
        startDT: Union[
            str, datetime.datetime, np.datetime64, pd.Timestamp, None,
        ] = None,
        endDT: Union[str, datetime.datetime, np.datetime64, pd.Timestamp, None,] = None,
        period: Union[str, None] = None,
        siteStatus: str = "active",
        **params,
    ):
        """Return Pandas DataFrame of NWIS IV data.

        Parameters
        ----------
        sites: str, list, pandas.Series, or numpy.Array, required
            Comma separated list of sites in string format or iterable.
        parameterCd: str, optional, default '00060' (Discharge)
            Comma separated list of parameter codes in string format.
        startDT: str, datetime.datetime, np.datetime64, pd.Timestamp, or None
                 optional, default None
            Observation record start time. If timezone information not provided,
            defaults to UTC.
        endDT: str, datetime.datetime, np.datetime64, pd.Timestamp, or None
                 optional, default None
        period: str, None
            Observation record for period until current time. Uses ISO 8601
            period time.
        siteStatus: str, optional, default 'active'
            Site status in string format
        params: 
            Additional parameters passed directly to service.

        Returns
        -------
        pandas.DataFrame :
            DataFrame in semi-WRES compatible format

        Examples
        --------
        >>> from evaluation_tools.nwis_client.iv import IVDataService
        >>> data = IVDataService.get(sites='01646500')
        """
        iv_data_service = cls()
        raw_data = iv_data_service.get_raw(
            sites,
            parameterCd=parameterCd,
            startDT=startDT,
            endDT=endDT,
            period=period,
            siteStatus=siteStatus,
            **params,
        )

        def list_to_df_helper(item: dict):
            values = item.pop("values")
            df = pd.DataFrame(values)

            for column_name, value in item.items():
                df[column_name] = value

            return df

        list_of_frames = map(list_to_df_helper, raw_data)

        # Concatenate list in single pd.DataFrame
        dfs = pd.concat(list_of_frames, ignore_index=True)

        # Convert values to numbers
        dfs.loc[:, "value"] = pd.to_numeric(dfs["value"], downcast="float")

        # Convert all times to UTC
        dfs["value_date"] = pd.to_datetime(
            dfs["dateTime"], utc=True, infer_datetime_format=True
        ).dt.tz_localize(None)

        # # Simplify variable name
        dfs["variable_name"] = dfs["variableName"].apply(
            iv_data_service.simplify_variable_name
        )

        # Sort DataFrame
        dfs = dfs.sort_values(["usgs_site_code", "measurement_unit", 
            "value_date"], ignore_index=True)

        # Fill NaNs
        dfs = dfs.fillna("")

        # Convert categories
        cols = (['variable_name', 'usgs_site_code', 'measurement_unit', 
            'qualifiers', 'series'])
        dfs[cols] = dfs[cols].astype(str)
        dfs[cols] = dfs[cols].astype(dtype='category')

        # Downcast floats
        df_float = dfs.select_dtypes(include=['float'])
        converted_float = df_float.apply(pd.to_numeric, downcast='float')
        dfs[converted_float.columns] = converted_float

        # DataFrame in semi-WRES compatible format
        return dfs[
            [
                "value_date",
                "variable_name",
                "usgs_site_code",
                "measurement_unit",
                "value",
                "qualifiers",
                "series",
            ]
        ]

    @classmethod
    def get_as_json(
        cls,
        sites: Union[str, List[str]],
        parameterCd: str = "00060",
        startDT: Union[
            str, datetime.datetime, np.datetime64, pd.Timestamp, None,
        ] = None,
        endDT: Union[str, datetime.datetime, np.datetime64, pd.Timestamp, None,] = None,
        period: Union[str, None] = None,
        siteStatus: str = "active",
        **params,
    ):
        # TODO Docstring
        iv_data_service = cls()
        return iv_data_service.get_raw(
            sites,
            parameterCd=parameterCd,
            startDT=startDT,
            endDT=endDT,
            period=period,
            siteStatus=siteStatus,
            **params,
        )

    def get_raw(
        self,
        sites: Union[str, List[str]],
        parameterCd: str = "00060",
        startDT: Union[
            str, datetime.datetime, np.datetime64, pd.Timestamp, None,
        ] = None,
        endDT: Union[str, datetime.datetime, np.datetime64, pd.Timestamp, None,] = None,
        period: Union[str, None] = None,
        siteStatus: str = "active",
        max_sites_per_request: int = 100,
        **params,
    ) -> List[requests.Response]:
        """
        Return raw requests data from the NWIS IV Rest API in a list.

        Parameters
        ----------
        sites: str, list, pandas.Series, or numpy.Array, required
            Comma separated list of sites in string format or iterable.
        parameterCd: str, optional, default '00060' (Discharge)
            Comma separated list of parameter codes in string format.
        startDT: str, datetime.datetime, np.datetime64, pd.Timestamp, or None
                 optional, default None
            Observation record start time. If timezone information not provided,
            defaults to UTC.
        endDT: str, datetime.datetime, np.datetime64, pd.Timestamp, or None
                 optional, default None
            Observation record end time. If timezone information not provided,
            defaults to UTC.
        period: str, None
            Observation record for period until current time. Uses ISO 8601
            period time.
        siteStatus: str, optional, default 'active'
            Site status in string format
        max_sites_per_request: int, optional, default 100
            Generally should not be changed. Maximum number of sites in any single
            `requests.get` call. Any number greater will cause sites to be divided
            evenly amongst parallel
            `requests.get` calls.
        kwargs: 
            Additional parameters passed directly to service.

        Returns
        -------
        List[requests.Response] :
            A list of retrieved requests objects

        Examples
        --------
        >>> from evaluation_tools.nwis_client import iv
        >>> service = iv.IVDataService()
        >>> data = service.get(sites='01646500')

        """
        # TODO: Add split multiprocessing logic for splitting by dates
        # Update parameters
        params.update(
            {"parameterCd": parameterCd, "siteStatus": siteStatus, "format": "json"}
        )

        # Handle startDT, endDT, and period optional arguments.
        if startDT is None and endDT is None and period is None:
            # Return IV case
            ...

        # startDT and endDT if included logic
        elif startDT and period is None:

            startDT = self._handle_date(startDT)
            params.update({"startDT": startDT})

            if endDT is not None:
                endDT = self._handle_date(endDT)
                params.update({"endDT": endDT})

        # period logic
        elif period and startDT is None and endDT is None:
            # TODO: Add iso 8601 period handler for splitting purposes
            if self._validate_period_string(period):
                params.update({"period": period})
            else:
                error_message = "`period` is not a valid ISO 8601 period string."
                raise TypeError(error_message)

        elif endDT or (period and endDT) or (period and startDT):
            # Throw
            error_message = (
                "Invalid set of datetime related arguments.\n"
                "Cannot only supply:\n"
                "`period` and `startDT` and `endDT`\n"
                "`period` and `startDT`\n"
                "`period` and `endDT`\n"
                "`endDT`"
            )
            raise KeyError(error_message)

        # Handle startDT and endDT if provided.
        # endDT is ignored if startDT not provided
        if startDT is not None:
            startDT = self._handle_date(startDT)
            params.update({"startDT": startDT})

            if endDT is not None:
                endDT = self._handle_date(endDT)
                params.update({"endDT": endDT})

        # Split sites into list
        if isinstance(sites, str):
            sites = sites.split(",")

        sites = self._unique(sites)

        # Threshold for divide and conquer requests
        if len(sites) > max_sites_per_request:
            response = self._multiprocessing_get(
                sites, params, partition_max_size=max_sites_per_request
            )

        else:
            # requests GET call
            response = self._get_and_handle_response(sites, parameters=params)

        return response

    @staticmethod
    def _handle_response(raw_response: requests.Response) -> List[dict]:
        """ From a raw response, return a list of extracted sites in dictionary form.
        Relevant dictionary keys are:
            "usgs_site_code"
            "variableName"
            "measurement_unit"
            "values"
            "series"

        Parameters
        ----------
        raw_response : requests.Response
            Request GET response

        Returns
        -------
        List[dict]
            A list of handled responses
        """
        # TODO: Speed test using orjson instead of native
        deserialized_response = raw_response.json()

        def extract_metadata(json_time_series):
            return {
                # Add site code
                "usgs_site_code": json_time_series["sourceInfo"]["siteCode"][0][
                    "value"
                ],
                # Add variable name
                "variableName": IVDataService.simplify_variable_name(
                    json_time_series["variable"]["variableName"]
                ),
                # Add units
                "measurement_unit": json_time_series["variable"]["unit"]["unitCode"],
            }

        flattened_data = []

        for response_value_timeSeries in deserialized_response["value"]["timeSeries"]:

            for indicies, site_data in enumerate(response_value_timeSeries["values"]):

                # Create general site metadata dictionary
                site_metadata = extract_metadata(response_value_timeSeries)

                # Add site time series values and its index number
                site_metadata.update({"values": site_data["value"], "series": indicies})
                flattened_data.append(site_metadata)

        return flattened_data

    @staticmethod
    def _unique(collection: Union[List, Tuple]) -> List:
        """ Return a list of the unique items from a collection.

        Parameters
        ----------
        collection : Union[List, Tuple]

        Returns
        -------
        List
            Returns a list of unique members of collection in their input type
        """
        if isinstance(collection, set):
            return list(collection)

        return list(set(collection))

    @staticmethod
    def simplify_variable_name(variable_name: str, split_delimiter: str = ",") -> str:
        """
        Split an input string by a delimiter and return only the first split
        result lowered.

        Parameters
        ----------
        variable_name : str
            String to simplify.

        split_delimiter : str
            Delimiter used to split data.

        Returns
        -------
        variable_name : str
            Simplified variable name

        """
        return variable_name.lower().split(split_delimiter)[0]

    def _get_and_handle_response(
        self, sites: List[str], parameters: dict = None,
    ):
        """ Thin procedural convenience wrapper for getting a GET response and
        handling a single response. Useful for multiprocessing.
        """
        if parameters is None:
            parameters = {}

        parameters["sites"] = sites
        # requests GET call
        response = self._restclient.get(
            parameters=parameters, headers=self._headers, parameter_delimeter=",",
        )
        return self._handle_response(response)

    def _multiprocessing_get(
        self, sites: List[str], parameters: dict, partition_max_size: int = 100
    ):
        chunked_sites = self._chunk_request(
            sites, partition_max_size=partition_max_size
        )

        # Retrieve responses
        with Pool(processes=self._procs) as pool:
            responses = pool.map(
                partial(self._get_and_handle_response, parameters=parameters),
                chunked_sites,
            )
            flattend_responses = [item for sublist in responses for item in sublist]

        return flattend_responses

    def _chunk_request(
        self, sites: Union[str, List[str]], partition_max_size: int = 100
    ) -> Tuple[List[str]]:

        # Split sites into list
        if isinstance(sites, str):
            # Split into list, notice variable name reuse
            sites = sites.split(",")

        elif not isinstance(sites, Iterable):
            error_message = (
                f"Did not receive a string or iterable collection.\nSites: {sites}"
            )
            raise TypeError(error_message)

        n_groups = (len(sites) // partition_max_size) + 1
        return np.array_split(sites, n_groups)

    def _handle_date(
        self, date: Union[str, datetime.datetime, np.datetime64, pd.Timestamp],
    ) -> str:
        """ Wrangle dates from a wide range of formats into a standard strftime
        string representation.

        Parameters
        ----------
        date : Union[str, datetime.datetime, np.datetime64, pd.Timestamp]
            Single date

        Returns
        -------
        str
            strftime string
        """

        if isinstance(date, str):
            try:
                handled_date = datetime.datetime.fromisoformat(date)
            except ValueError:
                # Parsing datetime aware strings is deprecated. However zulu strings
                # don't play well with fromisoformat.
                # See https://github.com/numpy/numpy/issues/6453
                handled_date = pd.to_datetime(date)
            except Exception as e:
                error_message = (
                    "Error parsing datetime string. Verify the format is "
                    f"consistent with iso8601.\n {e}"
                )
                BaseException(error_message)

        handled_date = pd.to_datetime(date)

        if handled_date.tz is None:
            # Default timezone to UTC if not provided
            handled_date = handled_date.tz_localize("UTC")

        return handled_date.strftime(self.datetime_format)

    def _validate_period_string(self, period: str) -> bool:
        """Validate if a string adheres to the duration format introduced in
        in ISO 8601.

        Parameters
        ----------
        period : str
            ISO 8601 period string.

        Returns
        -------
        bool
            True if validates against regex
        """
        pattern = "^(-?)P(?=\d|T\d)(?:(\d+)Y)?(?:(\d+)M)?(?:(\d+)([DW]))?(?:T(?:(\d+)H)?(?:(\d+)M)?(?:(\d+(?:\.\d+)?)S)?)?$"
        return True if re.fullmatch(pattern, period) else False

    @property
    def procs(self) -> int:
        """ Number of max multiprocessing processes """
        return self._procs

    @property
    def base_url(self) -> str:
        """ API Baseurl """
        return self._base_url

    @property
    def headers(self) -> dict:
        """ HTTP GET Headers """
        return self._headers

    @property
    def datetime_format(self) -> str:
        """ API's expected datetime format """
        return self._datetime_format

# Testing the T-Route FastAPI Extension

The following folder contains data files that are to be used to test the T-Route FastAPI code within src/app

To use these files, follow the steps below:

1. Copy the `test_compose.yaml` file in the base project dir (/t-route)
2. Run `docker compose -f test_compose.yaml up`
3. visit `localhost:8000/docs` in your browser
4. Enter the following parameters into the `/api/v1/flow_routing/v4` endpoint
- lid=CAGM7
- feature_id=2930769
- hy_id=1074884
- initial_start=0
- start_time=2024-08-24T00:00:00
- num_forecast_days=5
5. Click execute
6. A Status 200 code means the run ran, and test/api/data/troute_output will be populated in the `{lid}/` folder

## Docs:
See `doc/api_docs.md` for specific docs related to the T-Route API

## Running the Lower Colorado Example using the API
#### Written 9/6/24

To run the Lower Colorado example using the T-Route FastAPI instance, you can use the following parameters with the `/api/v1/flow_routing/v4` endpoint
- lid=LowerColorado
- feature_id=2420801
- hy_id=2420801
- initial_start=0
- start_time=2023-04-01T00:00:00
- num_forecast_days=2

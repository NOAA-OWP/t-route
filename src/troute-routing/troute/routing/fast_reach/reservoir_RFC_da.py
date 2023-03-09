


def _validate_RFC_data(lake_number, time_series, synthetic, time_series_file):
    use_RFC = True

    if all(synthetic)==1:
        use_RFC = False
        print(f"WARNING: RFC Forecast Time Series discharges for reservoir {lake_number} \
              using time series file {time_series_file} are all synthetic. \
              This reservoir will use level pool calculations instead.")
    elif any(time_series)<0:
        use_RFC = False
        print(f"WARNING: RFC Forecast Time Series discharges for reservoir {lake_number} \
              using time series file {time_series_file} contains missing or negative values. \
              This reservoir will use level pool calculations instead.")
    elif any(time_series)>=90000:
        use_RFC = False
        print(f"WARNING: RFC Forecast Time Series discharges for reservoir {lake_number} \
              using time series file {time_series_file} contain one or more values greater than or \
              equal to 90,000 Cubic Meters per Second (twice the Mississippi River historical peak flow). \
              This reservoir will use level pool calculations instead.")
    
    return use_RFC

def reservoir_RFC_da(lake_number, time_series, synthetic, time_series_file, routing_period, current_time,
                     rfc_forecast_persist_seconds, reservoir_type, ):
    use_RFC = _validate_RFC_data(lake_number, time_series, synthetic, time_series_file)

    if use_RFC and (routing_period<=3600) and current_time<=rfc_forecast_persist_seconds:
        if 1==1: #TODO: check if current model time is past res update time (is this every hour?)
            # If reservoir_type is 4 for CONUS RFC reservoirs
            if reservoir_type==4:
                # Set outflow to corresponding discharge from array
                outflow = time_series #TODO: figure out how to index the time series

            # Else reservoir_type 5 for for Alaska RFC glacier outflows
            else:
                # Set outflow to sum inflow and corresponding discharge from array
                outflow = inflow + time_series #TODO: figure out how to index the time series
            
            



        


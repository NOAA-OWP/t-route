


def _validate_RFC_data(lake_number, time_series, synthetic,):
    use_RFC = True

    if all(synthetic)==1:
        use_RFC = False
        print(f"WARNING: RFC Forecast Time Series discharges for reservoir {lake_number} \
              are all synthetic. \
              This reservoir will use level pool calculations instead.")
    elif any(time_series)<0:
        use_RFC = False
        print(f"WARNING: RFC Forecast Time Series discharges for reservoir {lake_number} \
              contains missing or negative values. \
              This reservoir will use level pool calculations instead.")
    elif any(time_series)>=90000:
        use_RFC = False
        print(f"WARNING: RFC Forecast Time Series discharges for reservoir {lake_number} \
              contain one or more values greater than or equal to 90,000 Cubic Meters per \
              Second (twice the Mississippi River historical peak flow). \
              This reservoir will use level pool calculations instead.")
    
    return use_RFC

def reservoir_RFC_da(lake_number, time_series, timeseries_idx, synthetic, routing_period, current_time,
                     update_time, DA_time_step, rfc_forecast_persist_seconds, reservoir_type, inflow, 
                     water_elevation, levelpool_outflow, levelpool_water_elevation, lake_area, 
                     max_water_elevation):
    use_RFC = _validate_RFC_data(lake_number, time_series, synthetic)

    if use_RFC and (routing_period<=3600) and current_time<=rfc_forecast_persist_seconds:
        if current_time >= update_time:
            # Advance update_time to the next timestep and time_series_idx to next index
            update_time += DA_time_step
            timeseries_idx += 1

        # If reservoir_type is 4 for CONUS RFC reservoirs
        if reservoir_type==4:
            # Set outflow to corresponding discharge from array
            outflow = time_series[timeseries_idx]

        # Else reservoir_type 5 for for Alaska RFC glacier outflows
        else:
            # Set outflow to sum inflow and corresponding discharge from array
            outflow = inflow + time_series[timeseries_idx]
        
        # Update water elevation
        new_water_elevation = water_elevation + ((inflow - outflow)/lake_area) * routing_period

        # Ensure that the water elevation is within the minimum and maximum elevation
        if new_water_elevation < 0.0:
            new_water_elevation = 0.0

        elif new_water_elevation > max_water_elevation:
            new_water_elevation = max_water_elevation
        
        # Set dynamic_reservoir_type to RFC Forecasts Type
        dynamic_reservoir_type = reservoir_type

        # Set the assimilated_value to corresponding discharge from array
        assimilated_value = time_series[timeseries_idx]

        # Check for outflows less than 0 and cycle backwards in the array until a
        # non-negative value is found. If all previous values are negative, then
        # use level pool outflow.
        if outflow < 0:
            missing_outflow_index = timeseries_idx

            while outflow < 0 and missing_outflow_index > 1:
                missing_outflow_index = missing_outflow_index - 1
                outflow = time_series[missing_outflow_index]
            
            if outflow < 0:
                # If reservoir_type is 4 for CONUS RFC reservoirs
                if reservoir_type == 4:
                    outflow = levelpool_outflow

                # Else reservoir_type 5 for for Alaska RFC glacier outflows
                else:
                    outflow = inflow

                # Update water elevation to levelpool water elevation
                new_water_elevation = levelpool_water_elevation

                # Set dynamic_reservoir_type to levelpool type
                dynamic_reservoir_type = 1

                # Set the assimilated_value to sentinel, -9999.0
                assimilated_value = -9999.0

                # Set the assimilated_source_file to empty string
                assimilated_source_file = ""

    else:
        # If reservoir_type is 4 for CONUS RFC reservoirs
        if reservoir_type == 4:
            outflow = levelpool_outflow

        # Else reservoir_type 5 for for Alaska RFC glacier outflows
        else:
            outflow = inflow

        # Update water elevation to levelpool water elevation
        new_water_elevation = levelpool_water_elevation

        # Set dynamic_reservoir_type to levelpool type
        dynamic_reservoir_type = 1

        # Set the assimilated_value to sentinel, -9999.0
        assimilated_value = -9999.0

        # Set the assimilated_source_file to empty string
        assimilated_source_file = ""
    
    return outflow, new_water_elevation, update_time, timeseries_idx, dynamic_reservoir_type, assimilated_value, assimilated_source_file
            



        


#ifndef RFC_STRUCTS_H
#define RFC_STRUCTS_H
/*
    C Structures
*/
#include "../../reach_structs.h"
typedef struct {
  int lake_number;
  float dam_length, area, max_depth;
  float orifice_area, orifice_coefficient, orifice_elevation;
  float weir_coefficient, weir_elevation, weir_length;
  float initial_fractional_depth, water_elevation;

  int reservoir_type;
  char reservoir_parameter_file;
  char start_date;
  char time_series_path;
  int forecast_lookback_hours;

  //Handle to operate hybrid res
  void* handle;
} _MC_Hybrid;

#endif //RFC_STRUCTS_H
#include <stdlib.h>
#include "../../reach_structs.h"
#include <stdio.h>
/* Level Pool Reservoir Interface */
extern void* get_lp_handle();

extern void init_lp(void* handle, float *water_elevation, float *lake_area, float *weir_elevation,
                    float *weir_coefficient, float *weir_length, float *dam_length, float *orifice_elevation,
                    float *orifice_coefficient, float *orifice_area, float *max_depth, int *lake_number);

extern void run_lp(void* handle, float *inflow, float *lateral_inflow,
                    float *water_elevation, float *outflow, float *routing_period);

extern void free_lp(void* handle);

init_levelpool_reach(_Reach* reach, int lake_number,
                          float dam_length, float area, float max_depth,
                          float orifice_area, float orifice_coefficient, float orifice_elevation,
                          float weir_coefficient, float weir_elevation, float weir_length,
                          float initial_fractional_depth, float water_elevation
)
{
  if( reach != NULL )
  {
    reach->reach.lp.lake_number = lake_number;
    reach->reach.lp.dam_length = dam_length;
    reach->reach.lp.area = area;
    reach->reach.lp.max_depth = max_depth;
    reach->reach.lp.orifice_area = orifice_area;
    reach->reach.lp.orifice_coefficient = orifice_coefficient;
    reach->reach.lp.orifice_elevation = orifice_elevation;
    reach->reach.lp.weir_coefficient = weir_coefficient;
    reach->reach.lp.weir_elevation = weir_elevation;
    reach->reach.lp.weir_length = weir_length;
    reach->reach.lp.initial_fractional_depth = initial_fractional_depth;

    if(water_elevation < 0){
      //Equation below is used in wrf-hydro
      printf("WARNING: LEVELPOOL USING COLDSTART WATER ELEVATION\n");
      fflush(stdout);
      reach->reach.lp.water_elevation = orifice_elevation + ((max_depth - orifice_elevation) * initial_fractional_depth);
    }
    else{
      reach->reach.lp.water_elevation = water_elevation;
    }

    reach->reach.lp.handle = get_lp_handle();
    init_lp(reach->reach.lp.handle, &reach->reach.lp.water_elevation, &reach->reach.lp.area,
                 &reach->reach.lp.weir_elevation, &reach->reach.lp.weir_coefficient, &reach->reach.lp.weir_length,
                 &reach->reach.lp.dam_length, &reach->reach.lp.orifice_elevation, &reach->reach.lp.orifice_coefficient,
                 &reach->reach.lp.orifice_area, &reach->reach.lp.max_depth, &reach->reach.lp.lake_number);
  }
}

void free_levelpool_reach(_Reach* reach)
{
  free_lp(reach->reach.lp.handle);
}

void route(_Reach* reach, float inflow, float lateral_inflow, float routing_period,
           float* outflow, float* water_elevation)
{
  run_lp(reach->reach.lp.handle, &inflow, &lateral_inflow, &reach->reach.lp.water_elevation, outflow, &routing_period);
  *water_elevation = reach->reach.lp.water_elevation;
}

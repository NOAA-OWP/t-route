#include <stdio.h>

extern void* get_lp_handle();

extern void init_lp(void* handle, float *water_elevation, float *lake_area, float *weir_elevation, float *weir_coefficient, float *weir_length, float *dam_length, float *orifice_elevation, float *orifice_coefficient, float *orifice_area, float *max_depth, int *lake_number);

extern void run_lp(void* handle, float *previous_timestep_inflow, float *inflow, float *lateral_inflow, float *water_elevation, float *outflow, float *routing_period, int *dynamic_reservoir_type, float *assimilated_value, char *assimilated_source_file);

extern void free_lp(void* handle);

int main(int argc, char** argv)
{
    float water_elevation = 0.0;
    float lake_area = 1.509490013122558594e+01;
    float weir_elevation = 9.626000022888183238e+00;
    float weir_coefficient = 0.4;
    float weir_length = 1.000000000000000000e+01;
    float dam_length = 10.0;
    float orifice_elevation = 7.733333269755045869e+00;
    float orifice_coefficient = 1.000000000000000056e-01;
    float orifice_area = 1.0;
    float max_depth = 9.960000038146972656e+00;
    int lake_number = 16944276;
    void* test_p = get_lp_handle();

    init_lp(test_p, &water_elevation, &lake_area, &weir_elevation, &weir_coefficient, &weir_length, &dam_length, &orifice_elevation, &orifice_coefficient, &orifice_area, &max_depth, &lake_number);

    float previous_timestep_inflow = 0.0;
    float inflow = 0.0;
    float lateral_inflow = 0.0;
    float outflow = 0.0;
    float routing_period = 300.0;
    int dynamic_reservoir_type = 1;
    float assimilated_value = 0.0;
    char assimilated_source_file[256];
    water_elevation = 9.73733330;
    
    run_lp(test_p, &previous_timestep_inflow, &inflow, &lateral_inflow, &water_elevation, &outflow, &routing_period, &dynamic_reservoir_type, &assimilated_value, assimilated_source_file);
    
    printf ("Outflow: %f\n", outflow);
    printf("Complete \n");
}


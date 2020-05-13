"""NHD Network traversal
A demonstration version of this code is stored in this Colaboratory notebook:
    https://colab.research.google.com/drive/1ocgg1JiOGBUl3jfSUPCEVnW5WNaqLKCD


## command line input order [debuglevel = 0,verbose = true,showtiming = true, supernetworks = supernetwork_options]
"""
import nhd_network_utilities as nnu
import recursive_print
import time
import os
import sys

try:
    arg1 = int(sys.argv[1])
    arg2 = sys.argv[2]
    arg3 = sys.argv[3]
    arg4 = sys.argv[4]
except:
    arg1 = 0
    arg2 = True
    arg3 = True
    arg4 = 'Brazos_LowerColorado_ge5'
    supernetwork_options = {
            'Pocono_TEST1'
            ,'Pocono_TEST2'
            , 'LowerColorado_Conchos_FULL_RES'
            , 'Brazos_LowerColorado_ge5'
            , 'Brazos_LowerColorado_FULL_RES'
            , 'Brazos_LowerColorado_Named_Streams'
            , 'CONUS_ge5'
            , 'Mainstems_CONUS'
            , 'CONUS_Named_Streams'
            , 'CONUS_FULL_RES_v20'
        }
    printout1 = 'command line arguements available: debuglevel = {-3,-2,-1, 0}, verbose = {True/False}, showtiming = {True/False}, supernetwork = ' + str(supernetwork_options)
    
def main():
    # find the path of the test scripts, several levels above the script path
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    test_folder = os.path.join(root, r'test')
    
    supernetworks = {}
    supernetworks.update({arg4:{}})
    # supernetworks.update({'Pocono_TEST1':{}})
    # supernetworks.update({'LowerColorado_Conchos_FULL_RES':{}}) 
    # supernetworks.update({'Brazos_LowerColorado_ge5':{}}) ##NHD Subset (Brazos/Lower Colorado)"""
    # supernetworks.update({'Brazos_LowerColorado_FULL_RES':{}}) 
    # supernetworks.update({'Brazos_LowerColorado_Named_Streams':{}}) 
    # supernetworks.update({'CONUS_ge5':{}}) ##NHD CONUS order 5 and greater"""
    # supernetworks.update({'Mainstems_CONUS':{}})
    # supernetworks.update({'CONUS_Named_Streams':{}})
    # supernetworks.update({'CONUS_FULL_RES_v20':{}}) # = False

    debuglevel = arg1
    verbose = arg2
    showtiming = arg3
    supernetwork = arg4

    for supernetwork in supernetworks:
        supernetworks[supernetwork] = nnu.set_supernetwork_data(
          supernetwork = supernetwork
            , geo_input_folder = os.path.join(test_folder, r'input', r'geo')
            , debuglevel = debuglevel
            , verbose = verbose
        )
        if debuglevel <= -1: 
            if verbose: print(f'\n\nSupernetwork:')
            print(f'{supernetwork}')
        if debuglevel <= -2: 
            if verbose: print(r'All items in the above supernetwork:')
            for k,v in supernetworks[supernetwork].items():
                 print(f"{{'{k}' : {v}}}")
        if showtiming: start_time = time.time()

        network_out_values = \
          nnu.get_nhd_connections(
            supernetworks[supernetwork]
            , debuglevel = debuglevel
            , verbose = verbose
        )

        recursive_print.print_basic_network_info(
          connections = network_out_values[0]
            , headwater_keys = network_out_values[3]
            , junction_keys = network_out_values[7]
            , terminal_keys = network_out_values[4]
            , terminal_code = supernetworks[supernetwork]['terminal_code']
            , verbose = verbose
        )

        if debuglevel <= -3: 
        # THE RECURSIVE PRINT IS NOT A GOOD IDEA WITH LARGE NETWORKS!!!
        # The `Pocono_TEST1` supernetwork is a good test case to run with 
        # the debuglevel set at -3. 
            recursive_print.print_connections(
                        headwater_keys = network_out_values[3]
                        , down_connections = network_out_values[0]
                        , up_connections = network_out_values[0]
                        , terminal_code = supernetworks[supernetwork]['terminal_code']
                        , terminal_keys = network_out_values[4]
                        , terminal_ref_keys = network_out_values[5]
                        , debuglevel = debuglevel
                        )
        if showtiming: print(f"Supernetwork `{supernetwork}` read and traversed\n... in %s seconds.\n\n" % (time.time() - start_time))
        print(printout1)
if __name__ == '__main__':
    main()

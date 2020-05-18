"""NHD Network traversal
A demonstration version of this code is stored in this Colaboratory notebook:
    https://colab.research.google.com/drive/1ocgg1JiOGBUl3jfSUPCEVnW5WNaqLKCD

"""
import nhd_network_utilities as nnu
import recursive_print
import time
import os
import argparse

def _handle_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--debuglevel',
                        help='Set the debuglevel',
                        dest='debuglevel',
                        default=0)
    parser.add_argument('-v, --verbose',
                        help='Verbose output (leave blank for quiet output)',
                        dest='verbose',
                        action='store_true')
    parser.add_argument('-t, --showtiming',
                        help='Set the showtiming - leave blank for False',
                        dest='showtiming',
                        action='store_true')
    parser.add_argument('-n, --supernetwork',
                        help='List of supernetworks (Pocono_TEST1, Pocono_TEST2, LowerColorado_Conchos_FULL_RES, Brazos_LowerColorado_ge5, Brazos_LowerColorado_FULL_RES, Brazos_LowerColorado_Named_Streams, CONUS_ge5, Mainstems_CONUS, CONUS_Named_Streams, CONUS_FULL_RES_v20', 
                        # TODO: accept multiple or a Path (argparse Action perhaps)
                        # action='append',
                        # nargs=1,
                        dest='supernetwork',
                        default='Pocono_TEST1')

    return parser.parse_args()


def main():

    args = _handle_args()

    # find the path of the test scripts, several levels above the script path
    root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    test_folder = os.path.join(root, r'test')
    
    debuglevel = -1 * int(args.debuglevel)
    verbose = args.verbose
    showtiming = args.showtiming
    supernetworks = {str(args.supernetwork):{}}

    # supernetworks = {}
    # supernetworks.update({'Pocono_TEST1':{}})
    # supernetworks.update({'LowerColorado_Conchos_FULL_RES':{}}) 
    # supernetworks.update({'Brazos_LowerColorado_ge5':{}}) ##NHD Subset (Brazos/Lower Colorado)"""
    # supernetworks.update({'Brazos_LowerColorado_FULL_RES':{}}) 
    # supernetworks.update({'Brazos_LowerColorado_Named_Streams':{}}) 
    # supernetworks.update({'CONUS_ge5':{}}) ##NHD CONUS order 5 and greater"""
    # supernetworks.update({'Mainstems_CONUS':{}})
    # supernetworks.update({'CONUS_Named_Streams':{}})
    # supernetworks.update({'CONUS_FULL_RES_v20':{}}) # = False

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
        
if __name__ == '__main__':
    main()

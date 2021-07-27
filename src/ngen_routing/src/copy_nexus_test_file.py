import pandas as pd
import json
import troute.nhd_network as nhd_network
import troute.nhd_io as nhd_io
from shutil import copyfile


#ngen_network_df = nhd_io.read_geopandas( "........../ngen_pybind/crosswalk_mapping/flowpath_data.ge


#mapping_dict = json.read("../../../test/input/next_gen/nexus_id_to_downstream_comid_mapping_subset1.json")


with open("../../../test/input/next_gen/nexus_id_to_downstream_comid_mapping_subset1.json") as f:

    mapping_dict = json.load(f)


print (mapping_dict)


#nex_set = 
nex_list = []

for key in mapping_dict:

   print (key)

   if key not in nex_list:
       nex_list.append(key)
   else:
       print ("above key already in nex_list")



#nex_set = set(nex_list)

#print ("nex_set")
#print (nex_set)

copy_nex_file = "/apd_common/anthro/david_code/ngen/nex-11111_copy.csv"

path = "/apd_common/anthro/david_code/ngen/"

print ("each nex in list ----")
for nex in nex_list:
    print (nex)

    new_file_name = path + "nex-" + str(nex) + "_output.csv"

    print (new_file_name)

    copyfile(copy_nex_file, new_file_name)




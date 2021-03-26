import numpy as np
import pathlib
import sys
import yaml
import os

ENV_IS_CL = False
if ENV_IS_CL:
    root = pathlib.Path("/", "content", "t-route")
elif not ENV_IS_CL:
    root = pathlib.Path("../..").resolve()
    sys.path.append("fast_reach")

## network and reach utilities
import diffusive

# open input data from yaml
with open(os.path.join(root,'src/python_routing_v02_diffusive/diff_inputs.yml')) as file:
    data_list = yaml.load(file, Loader=yaml.Loader)
    
print(data_list["frnw_g"])
    
out_q, out_elv = diffusive.compute_diffusive(data_list["dtini_g"],
                                data_list["t0_g"],
                                data_list["tfin_g"],
                                data_list["saveinterval_g"],
                                data_list["saveinterval_ev_g"],
                                data_list["dt_ql_g"],
                                data_list["dt_ub_g"],
                                data_list["dt_db_g"],
                                data_list["nts_ql_g"],
                                data_list["nts_ub_g"],
                                data_list["nts_db_g"],
                                data_list["mxncomp_g"],
                                data_list["nrch_g"],
                                np.asfortranarray(data_list["z_ar_g"]),
                                np.asfortranarray(data_list["bo_ar_g"]),
                                np.asfortranarray(data_list["traps_ar_g"]),
                                np.asfortranarray(data_list["tw_ar_g"]),
                                np.asfortranarray(data_list["twcc_ar_g"]),
                                np.asfortranarray(data_list["mann_ar_g"]),
                                np.asfortranarray(data_list["manncc_ar_g"]),
                                np.asfortranarray(data_list["so_ar_g"]),
                                np.asfortranarray(data_list["dx_ar_g"]),
                                data_list["nhincr_m_g"],
                                data_list["nhincr_f_g"],
                                np.asfortranarray(data_list["ufhlt_m_g"]),
                                np.asfortranarray(data_list["ufqlt_m_g"]),
                                np.asfortranarray(data_list["ufhlt_f_g"]),
                                np.asfortranarray(data_list["ufqlt_f_g"]),
                                data_list["frnw_col"],
                                np.asfortranarray(data_list["frnw_g"].astype('double')),
                                np.asfortranarray(data_list["qlat_g"]),
                                np.asfortranarray(data_list["ubcd_g"]),
                                np.asfortranarray(data_list["dbcd_g"]),
                                data_list["cfl_g"],
                                data_list["theta_g"],
                                data_list["tzeq_flag_g"],
                                data_list["y_opt_g"],
                                data_list["so_llm_g"],
                                data_list["ntss_ev_g"],
                               )

print(out_elv)

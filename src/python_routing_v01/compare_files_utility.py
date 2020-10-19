import pathlib
import os
import sys 
import argparse

def pr_check(root_folder_PR2,x):
    root = pathlib.Path("../../../../").resolve()
    sys.setrecursionlimit(4000)
    sys.path.append(os.path.join(root, root_folder_PR2, r"t-route", r"src", r"python_framework_v01"))
    sys.path.append(os.path.join(root, root_folder_PR2, r"t-route", r"src", r"python_framework_v02"))
    sys.path.append(os.path.join(root, root_folder_PR2, r"t-route", r"src", r"python_routing_v01"))
    fortran_routing_dir = os.path.join(
        root, root_folder_PR2, r"t-route", r"src", r"fortran_routing", r"mc_pylink_v00", r"MC_singleSeg_singleTS"
    )
    fortran_reservoir_dir = os.path.join(
        root, root_folder_PR2, r"t-route", r"src", r"fortran_routing", r"mc_pylink_v00", r"Reservoir_singleTS"
    )
    routing_directory = os.path.join(root, root_folder_PR2, r"t-route", r"src", r"python_routing_v01")
    routing_directory_v02 = os.path.join(root, root_folder_PR2, r"t-route", r"src", r"python_routing_v02")
    sys.path.append(fortran_routing_dir)

    os.chdir(routing_directory)
    
    debuglevel = x
    COMPILE = True
    if COMPILE:
        try:
            import subprocess
            from subprocess import PIPE, run

            python_compile_call = []
            python_compile_call.append(r"python3")
            python_compile_call.append(r"compute_nhd_routing_SingleSeg.py")
            python_compile_call.append(r"--verbose")
            python_compile_call.append(r'--debuglevel '+ x)
        except Exception as e:
            # print(e)
            if debuglevel <= -4:
                traceback.print_exc()

    result_sseg2 = subprocess.run(python_compile_call, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    # print(result_sseg2.returncode, result_sseg2.stdout, result_sseg2.stderr)
    result_sseg2 = str(subprocess.run(python_compile_call, stdout=PIPE, stderr=PIPE, universal_newlines=True))
    os.chdir(routing_directory_v02)

    debuglevel = x
    COMPILE = True
    if COMPILE:
        try:
            import subprocess
            from subprocess import PIPE, run

            python_compile_call = []
            python_compile_call.append(r"python3")
            python_compile_call.append(r"compute_nhd_routing_SingleSeg_v02.py")
            python_compile_call.append(r"--verbose")
            python_compile_call.append(r'--debuglevel '+ -x)
        except Exception as e:
            # print(e)
            if debuglevel <= -4:
                traceback.print_exc()


    result_sseg_v02_2 = subprocess.run(python_compile_call, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    # print(result_sseg_v02_2.returncode, result_sseg_v02_2.stdout, result_sseg_v02_2.stderr)
    result_sseg_v02_2 = str(subprocess.run(python_compile_call, stdout=PIPE, stderr=PIPE, universal_newlines=True))
    # print(root)
    # print(sys.path)
    return result_sseg2, result_sseg_v02_2
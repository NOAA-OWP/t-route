import pathlib
import os
import sys 
import argparse
import compare_files_utility as cfu

def _handle_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "-m",
        "--root_folder_PR1",
        help="Choose the folder of the top most root location of your first PR request comparison.",
        dest="root_folder_PR1",
        default="github",
    )
    parser.add_argument(
        "-n",
        "--root_folder_PR2",
        help="Choose the folder of the top most root location of your first PR request comparison.",
        dest="root_folder_PR2",
        default="tempforpr",
    )
    return parser.parse_args()


args = _handle_args()
root_folder_PR2 = args.root_folder_PR2

root = pathlib.Path("../../../../").resolve()
args = _handle_args()
root_folder_PR1 = args.root_folder_PR1
sys.setrecursionlimit(4000)
sys.path.append(os.path.join(root, root_folder_PR1, r"t-route", r"src", r"python_framework_v01"))
sys.path.append(os.path.join(root, root_folder_PR1, r"t-route", r"src", r"python_framework_v02"))
sys.path.append(os.path.join(root, root_folder_PR1, r"t-route", r"src", r"python_routing_v01"))
fortran_routing_dir = os.path.join(
    root, root_folder_PR1, r"t-route", r"src", r"fortran_routing", r"mc_pylink_v00", r"MC_singleSeg_singleTS"
)
fortran_reservoir_dir = os.path.join(
    root, root_folder_PR1, r"t-route", r"src", r"fortran_routing", r"mc_pylink_v00", r"Reservoir_singleTS"
)
sys.path.append(fortran_routing_dir)
routing_directory = os.path.join(root, root_folder_PR2, r"t-route", r"src", r"python_routing_v01")
routing_directory_v02 = os.path.join(root, root_folder_PR2, r"t-route", r"src", r"python_routing_v02")

for x in range(0,4):
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


    result_sseg1 = subprocess.run(python_compile_call, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    # print(result_sseg1.returncode, result_sseg1.stdout, result_sseg1.stderr)
    result_sseg1 = str(subprocess.run(python_compile_call, stdout=PIPE, stderr=PIPE, universal_newlines=True))
    # print(root)
    # print(sys.path)
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


    result_sseg_v02_1 = subprocess.run(python_compile_call, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    # print(result_sseg_v02_1.returncode, result_sseg_v02_1.stdout, result_sseg_v02_1.stderr)
    result_sseg_v02_1 = str(subprocess.run(python_compile_call, stdout=PIPE, stderr=PIPE, universal_newlines=True))
    # print(sys.path)
    result_sseg2, result_sseg_v02_2 = cfu.pr_check(root_folder_PR2,x)

    if result_sseg1 == result_sseg2:
        print("PASS - results of compute_nhd_routing_SingleSeg.py match : debuglevel " + str(x) )
    else:
        print("FAIL - results of compute_nhd_routing_SingleSeg.py are not equivalent : debuglevel " + str(x) )

    if result_sseg_v02_1 == result_sseg_v02_2:
        print("PASS - results of compute_nhd_routing_SingleSeg_v02.py match: debuglevel " + str(x) )
    else:
        print("FAIL - results of compute_nhd_routing_SingleSeg_v02.py are not equivalent : debuglevel " + str(x) )



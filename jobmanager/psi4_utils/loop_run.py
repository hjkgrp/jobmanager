import os
import json
from jobmanager.psi4_utils.run import run_b3lyp, run_general

psi4_config = {'bashrc': '/home/crduan/.bashrc',
               'conda_env': '/home/crduan/miniconda/envs/mols_py36'}
with open("psi4_config.json", "r") as f:
    psi4_config.update(json.load(f))
success_count = 0
print("===b3lyp===")
if not os.path.isdir("b3lyp"):
    success = run_b3lyp(psi4_config, rundir="b3lyp")
    print("success: ", success)
    if success:
        success_count += 1
else:
    print("folder exists.")
    files = os.listdir("b3lyp")
    resubed = False
    functional = "b3lyp"
    if not os.path.isfile(functional + "/output.dat"):
        resubed = True
    else:
        with open(functional + "/output.dat", "r") as fo:
            txt = "".join(fo.readlines())
        if "==> Iterations <==" not in txt:
            resubed = True
    if resubed:
        print("previous errored out. resubmitting...")
        success = run_b3lyp(psi4_config, rundir="b3lyp")
        print("success: ", success)
        if success:
            success_count += 1
    else:
        with open(functional + "/output.dat", "r") as fo:
            txt = "".join(fo.readlines())
        if 'PsiException: Could not converge SCF iterations' not in txt and os.path.isfile(functional + "/wfn.180.npy"):
            print("success: ", True)
            success_count += 1
for ii, functional in enumerate(psi4_config["functional"]):
    print("===%d: %s===" % (ii, functional))
    if not os.path.isdir(functional.replace("(", "l-").replace(")", "-r")):
        os.makedirs(functional.replace("(", "l-").replace(")", "-r"))
        success = run_general(psi4_config, functional)
        print("success: ", success)
        if success:
            success_count += 1
    else:
        print("folder exists.")
        files = os.listdir(functional.replace("(", "l-").replace(")", "-r"))
        resubed = False
        if not os.path.isfile(functional.replace("(", "l-").replace(")", "-r") + "/output.dat"):
            resubed = True
        else:
            with open(functional.replace("(", "l-").replace(")", "-r") + "/output.dat", "r") as fo:
                txt = "".join(fo.readlines())
            if "==> Iterations <==" not in txt or (not (("@DF-UKS iter" in txt) or ("@DF-RKS iter" in txt) or ("@DF-UHF iter" in txt) or ("@DF-RHF iter" in txt))):
                resubed = True
        if resubed:
            print("previous errored out. resubmitting...")
            success = run_general(psi4_config, functional)
            print("success: ", success)
            if success:
                success_count += 1
        else:
            with open(functional.replace("(", "l-").replace(")", "-r") + "/output.dat", "r") as fo:
                txt = "".join(fo.readlines())
            if 'PsiException: Could not converge SCF iterations' not in txt and os.path.isfile(functional + "/wfn.180.npy"):
                print("success: ", True)
                success_count += 1
print("total successful jobs : %d/ %d." %
      (success_count, len(psi4_config["functional"]) + 1))

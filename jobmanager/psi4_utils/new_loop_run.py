#Runs the Psi4 workflow. From a directory with the Psi4 config and the TC molden,
#runs Psi4 single points with the specified functionals for each subfolder.
#Calculations are run from the subfolders, specifying geometries to run on,
#and are stored in subsubdirectories corresponding to the functional.

import os
import json
from jobmanager.psi4_utils.psi4_utils import Psi4Utils

with open("psi4_config.json", "r") as f:
    psi4_config = json.load(f)
success_count = 0
psi4_utils = Psi4Utils(psi4_config)

#Run the initial B3LYP calculation from the provided molden
print("===b3lyp===")
if not os.path.isdir("b3lyp"):
    success = psi4_utils.run_b3lyp(rundir="b3lyp")
    print("success: ", success)
    if success:
        success_count += 1
else:
    #If B3LYP calculation already attempted, check if converged
    print("folder exists.")
    files = os.listdir("b3lyp")
    resubed = False
    functional = "b3lyp"
    #Need to resubmit if output not present or if iterations not reached
    if not os.path.isfile(functional + "/output.dat"):
        resubed = True
    else:
        with open(functional + "/output.dat", "r") as fo:
            txt = "".join(fo.readlines())
        if "==> Iterations <==" not in txt:
            resubed = True
    if resubed:
        #Resubmit the B3LYP calculation if it failed
        print("previous errored out. resubmitting...")
        success = psi4_utils.run_b3lyp(rundir="b3lyp")
        print("success: ", success)
        if success:
            success_count += 1
    else:
        #If checks above pass, ensure no SCF error in B3LYP calculation and that .wfn written
        with open(functional + "/output.dat", "r") as fo:
            txt = "".join(fo.readlines())
        if 'PsiException: Could not converge SCF iterations' not in txt and os.path.isfile(functional + "/wfn.180.npy"):
            print("success: ", True)
            success_count += 1

#for all other functionals, run single points from the B3LYP reference
for ii, functional in enumerate(psi4_config["functional"]):
    print("===%d: %s===" % (ii, functional))
    #If the functional does not already have a folder, attempt the calculation
    if not os.path.isdir(functional.replace("(", "l-").replace(")", "-r")):
        os.makedirs(functional.replace("(", "l-").replace(")", "-r"))
        success = psi4_utils.run_general(functional)
        print("success: ", success)
        if success:
            success_count += 1
    #If the calculation already attempted, check for convergence
    else:
        print("folder exists.")
        files = os.listdir(functional.replace("(", "l-").replace(")", "-r"))
        resubed = False
        #Need resubmission if no output or iterations not in the text
        if not os.path.isfile(functional.replace("(", "l-").replace(")", "-r") + "/output.dat"):
            resubed = True
        else:
            with open(functional.replace("(", "l-").replace(")", "-r") + "/output.dat", "r") as fo:
                txt = "".join(fo.readlines())
            if "==> Iterations <==" not in txt or (not (("@DF-UKS iter" in txt) or ("@DF-RKS iter" in txt) or ("@DF-UHF iter" in txt) or ("@DF-RHF iter" in txt))):
                resubed = True
        if resubed:
            #Resubmit job
            print("previous errored out. resubmitting...")
            success = psi4_utils.run_general(functional)
            print("success: ", success)
            if success:
                success_count += 1
        else:
            #Check that no SCF exception and that a .wfn file written
            with open(functional.replace("(", "l-").replace(")", "-r") + "/output.dat", "r") as fo:
                txt = "".join(fo.readlines())
            #wfn file will not be written since run_general has default return_wfn=False
            if 'PsiException: Could not converge SCF iterations' not in txt: # and os.path.isfile(functional + "/wfn.180.npy"):
                print("success: ", True)
                success_count += 1
print("total successful jobs : %d/ %d." %
      (success_count, len(psi4_config["functional"]) + 1))

#Attempts to converge a calculation by converging it at different
#HFX levels. Specifically, increases HFX to 20% to correspond to B3LYP
#and then steps it down to 0 to match non-hybrid functionals.

import os
import json
import shutil
from jobmanager.psi4_utils.psi4_utils import Psi4Utils

with open("psi4_config.json", "r") as f:
    psi4_config = json.load(f)
#List of HFX amounts to use
alphalist = [20, 15, 10, 5, 2]
rescued = 0
success_count = 0
psi4_utils = Psi4Utils(psi4_config)

for ii, functional in enumerate(psi4_config["functional"]):
    failed = False
    print("===%d: %s===" % (ii, functional))
    with open(functional + "/output.dat", "r") as fo:
        txt = "".join(fo.readlines())
    #Needs rescuing if SCF iterations did not converge or empty output
    #Simpler causes of errors (no output, no iterations) already caught in loop_run
    if 'PsiException: Could not converge SCF iterations' in txt or txt == "":
        success = False
        for jj, alpha in enumerate(alphalist):
            #Try starting the calculation the corresponding alpha
            print(alpha)
            if jj == 0:
                #For 20% HFX, use default b3lyp as initial guess
                wfn = "b3lyp/wfn.180.npy"
            else:
                #For non-20%, use the functional converged at the previous HFX step
                wfn = "%s/wfn.180.npy" % (functional + "-%d" %
                                             alphalist[jj-1])
            #If the calculation at this HFX has not been tried before and the previous calculation did not fail
            if not os.path.isdir(functional + "-%d" % alpha) and not failed:
                os.makedirs(functional + "-%d" % alpha)
                success = psi4_utils.run_general_hfx(functional, hfx=alpha, wfn=wfn)
                print("success: ", success)
                if success:
                    success_count += 1
                else:
                    failed = True
            else:
                #If the calculation (at jj HFX) has been tried before (or a previous calculation has failed)
                #Note: if a previous calculation fails, all subsequent calculations will also fail since wfn will not be found
                print("attempted rescue: ", functional)
                resubed = False
                #Need resubmission if no output or if iterations not reached in output
                if not os.path.isfile(functional + "-%d" % alpha + "/output.dat"):
                    resubed = True
                else:
                    with open(functional + "-%d" % alpha + "/output.dat", "r") as fo:
                        txt = "".join(fo.readlines())
                    if "==> Iterations <==" not in txt:
                        resubed = True
                #Resubmit calculation using the (jj-1)th HFX as a reference
                if resubed and os.path.isfile(wfn):
                    print("previously errored out. resubmitting...")
                    success = psi4_utils.run_general_hfx(functional, hfx=alpha, wfn=wfn)
                    print("success: ", success)
                    if success:
                        success_count += 1
                    else:
                        failed = True
                #If the error was not due to SCF error and no wfn file written
                elif ('PsiException: Could not converge SCF iterations') not in txt and (not os.path.isfile(functional + "-%d" % alpha + "/wfn.180.npy")):
                    #Try running the calculation again, move the old output to -timeout so it can be checked later if desired
                    _functional = functional + "-%d" % alpha
                    if not os.path.isfile(_functional + "/output-timeout.dat"):
                        shutil.copy(_functional + "/output.dat",
                                    _functional + "/output-timeout.dat")
                        print("Time out, direct resub...")
                        success = psi4_utils.run_general_hfx(functional, hfx=alpha, wfn=wfn)
                        print("success: ", success)
                        if success:
                            success_count += 1
                        else:
                            failed = True
                    else:
                        #Only try resubmission due to SCF errors once, if already tried, give up
                        failed = True
                        print("Already submit once for timeout.")
                else:
                    print("give up resubmission.")
        if not failed:
            rescued += 1
            print("rescued: ", functional)
    else:
        rescued += 1
        print("%s has already succeeded!" % functional)
print("all finished!")

import os
import json
import shutil
import numpy as np
from jobmanager.psi4_utils.run import run_general_hfx

psi4_config = {'bashrc': '/home/crduan/.bashrc',
               'conda_env': '/home/crduan/miniconda/envs/mols_py36'}
with open("psi4_config.json", "r") as f:
    psi4_config.update(json.load(f))
alphalist = [20, 15, 10, 5, 2]
rescued = 0
success_count = 0
basedir = os.getcwd()
for ii, functional in enumerate(psi4_config["functional"]):
    failed = False
    print("===%d: %s===" % (ii, functional))
    with open(functional + "/output.dat", "r") as fo:
        txt = "".join(fo.readlines())
    if 'PsiException: Could not converge SCF iterations' in txt or txt == "":
        success = False
        e0 = np.nan
        for jj, alpha in enumerate(alphalist):
            print(alpha)
            if jj == 0:
                wfn = "../b3lyp/wfn.180.npy"
                wfn_b = "b3lyp/wfn.180.npy"
            else:
                wfn = "../%s/wfn.180.npy" % (functional + "-%d" %
                                             alphalist[jj-1])
                wfn_b = "%s/wfn.180.npy" % (functional + "-%d" %
                                            alphalist[jj-1])
            if not os.path.isdir(functional + "-%d" % alpha) and not failed:
                os.makedirs(functional + "-%d" % alpha)
                shutil.copyfile(psi4_config['xyzfile'], functional + "-%d" %
                                alpha + '/' + psi4_config['xyzfile'])
                success = run_general_hfx(psi4_config, functional, hfx=alpha, wfn=wfn)
                print("success: ", success)
                if success:
                    success_count += 1
                else:
                    failed = True
            else:
                print("attempted rescue: ", functional)
                resubed = False
                if not os.path.isfile(functional + "-%d" % alpha + "/output.dat"):
                    resubed = True
                else:
                    with open(functional + "-%d" % alpha + "/output.dat", "r") as fo:
                        txt = "".join(fo.readlines())
                    if "==> Iterations <==" not in txt:
                        resubed = True
                if resubed and os.path.isfile(wfn_b):
                    print("previously errored out. resubmitting...")
                    success = run_general_hfx(psi4_config, functional, hfx=alpha, wfn=wfn)
                    print("success: ", success)
                    if success:
                        success_count += 1
                    else:
                        failed = True
                elif ('PsiException: Could not converge SCF iterations') not in txt and (not os.path.isfile(functional + "-%d" % alpha + "/wfn.180.npy")):
                    _functional = functional + "-%d" % alpha
                    if not os.path.isfile(_functional + "/output-timeout.dat"):
                        shutil.copy(_functional + "/output.dat",
                                    _functional + "/output-timeout.dat")
                        print("Time out, direct resub...")
                        success = run_general_hfx(psi4_config, functional, hfx=alpha, wfn=wfn)
                        print("success: ", success)
                        if success:
                            success_count += 1
                        else:
                            failed = True
                    else:
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

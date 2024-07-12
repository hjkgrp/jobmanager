#From a directory, will converge the B3LYP wavefunctions of all derivative jobs.
#The first derivative job is converged from the TC molden, while the rest are 
#converged from the .wfn of the first converged job.
#Should be launched from the parent directory, i.e., the folder containing
#the subfolders containing the trigger for derivative jobs.

import os
import json
from jobmanager.psi4_utils.derivative_utils import DerivativeUtils

#from jobmanager.psi4_utils.run import run_b3lyp, run_general
#from jobmanager.psi4_utils.derivative import derivative_tree, get_wfn_path
#from jobmanager.psi4_utils.stable_run import run_with_check


success_count = 0
with open("psi4_config.json", "r") as f:
    psi4_config = json.load(f)
derivative_utils = DerivativeUtils()

jobs = derivative_utils.derivative_tree(path="./", trigger=psi4_config["trigger"])
print("jobs:", jobs)
# ---run first job from scratch---
success_count = derivative_utils.run_with_check(job=jobs[0], psi4_config=psi4_config,
                               success_count=success_count, run_func='run_b3lyp', error_scf=True)
# ---run other jobs using the wfn from the previous step---
for ii, job in enumerate(jobs[1:]):
    psi4_config["wfnfile"] = derivative_utils.get_wfn_path(jobs, ii+1)
    success_count = derivative_utils.run_with_check(job=job, psi4_config=psi4_config,
                                   success_count=success_count, run_func='run_general',
                                   error_scf=True)
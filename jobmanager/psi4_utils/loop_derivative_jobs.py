import os
import json
from jobmanager.psi4_utils.run import run_b3lyp, run_general
from jobmanager.psi4_utils.derivative import derivative_tree, get_wfn_path
from jobmanager.psi4_utils.stable_run import run_with_check


basedir = os.getcwd()
success_count = 0
psi4_config = {'bashrc': '/home/crduan/.bashrc',
               'conda_env': '/home/crduan/miniconda/envs/mols_py36'}
with open("psi4_config.json", "r") as f:
    psi4_config.update(json.load(f))
jobs = derivative_tree(path="./", trigger=psi4_config["trigger"])
print("jobs:", jobs)
# ---run first job from scratch---
success_count = run_with_check(job=jobs[0], basedir=basedir, psi4_config=psi4_config,
                               success_count=success_count, run_func=run_b3lyp, error_scf=True)
# ---run other jobs using the wfn from the previous step---
for ii, job in enumerate(jobs[1:]):
    psi4_config["wfnfile"] = get_wfn_path(jobs, ii+1)
    success_count = run_with_check(job=job, basedir=basedir, psi4_config=psi4_config,
                                   success_count=success_count, run_func=run_general,
                                   error_scf=True)

import os
import types


def run_with_check(job: str, basedir: str, psi4_config: dict,
                   success_count: int, run_func: types.FunctionType,
                   error_scf: bool = True):
    print("====running on====: ", job)
    os.chdir(job)
    success = False
    if not os.path.isdir("b3lyp"):
        success = run_func(psi4_config, return_wfn=True)
        print("success: ", success)
        if success:
            success_count += 1
    else:
        print("folder exists.")
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
            success = run_func(psi4_config, return_wfn=True)
            print("success: ", success)
            if success:
                success_count += 1
        else:
            with open(functional + "/output.dat", "r") as fo:
                txt = "".join(fo.readlines())
            if 'PsiException: Could not converge SCF iterations' not in txt and os.path.isfile(functional + "/wfn.180.npy"):
                print("success: ", True)
                success = True
                success_count += 1
    os.chdir(basedir)
    if not success and error_scf:
        raise ValueError(
            "Failed on the job: %s. Other derivative jobs won't run." % job)
    return success_count

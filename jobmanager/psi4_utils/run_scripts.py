import os
import json
import shutil
from jobmanager.psi4_utils.psi4_utils import Psi4Utils
from jobmanager.psi4_utils.derivative_utils import DerivativeUtils

class RunScripts:
    """
    Contains scripts that, when launched from a subdirectory (the folder corresponding to a structure),
    will run the Psi4 workflow for that structure.
    """

    def __init__(self) -> None:
        pass

    def loop_run(self, rundir='$SGE_O_WORKDIR'):
        """
        Runs the Psi4 workflow. From a directory with the Psi4 config and the TC molden,
        runs Psi4 single points with the specified functionals for each subfolder.
        Calculations are run from the subfolders, specifying geometries to run on,
        and are stored in subsubdirectories corresponding to the functional.
        """

        #To deal with bash variable paths
        if rundir[0] == '$':
            rundir = os.environ[rundir[1:]]

        with open(rundir + "/../psi4_config.json", "r") as f:
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

    def loop_rescue(self, rundir='$SGE_O_WORKDIR'):
        """
        Attempts to converge a calculation by converging it at different
        HFX levels. Specifically, increases HFX to 20% to correspond to B3LYP
        and then steps it down to 0 to match non-hybrid functionals.
        """

        #To deal with bash variable paths
        if rundir[0] == '$':
            rundir = os.environ[rundir[1:]]

        with open(rundir + "/../psi4_config.json", "r") as f:
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
        
    def loop_derivative_jobs(self, rundir='$SGE_O_WORKDIR'):
        """
        From a directory, will converge the B3LYP wavefunctions of all derivative jobs.
        The first derivative job is converged from the TC molden, while the rest are 
        converged from the .wfn of the first converged job.
        Should be launched from the parent directory, i.e., the folder containing
        the subfolders containing the trigger for derivative jobs.
        """

        #To deal with bash variable paths
        if rundir[0] == '$':
            rundir = os.environ[rundir[1:]]

        basedir = os.getcwd()
        success_count = 0
        with open(rundir + "/psi4_config.json", "r") as f:
            psi4_config = json.load(f)
        derivative_utils = DerivativeUtils()

        jobs = derivative_utils.derivative_tree(path="./", trigger=psi4_config["trigger"])
        print("jobs:", jobs)
        # ---run first job from scratch---
        os.chdir(jobs[0])
        success_count = derivative_utils.run_with_check(job=jobs[0], psi4_config=psi4_config,
                                    success_count=success_count, run_func='run_b3lyp', error_scf=True)
        os.chdir(basedir)
        # ---run other jobs using the wfn from the previous step---
        for ii, job in enumerate(jobs[1:]):
            psi4_config["wfnfile"] = derivative_utils.get_wfn_path(jobs, ii+1)
            os.chdir(job)
            success_count = derivative_utils.run_with_check(job=job, psi4_config=psi4_config,
                                        success_count=success_count, run_func='run_general',
                                        error_scf=True)
            os.chdir(basedir)
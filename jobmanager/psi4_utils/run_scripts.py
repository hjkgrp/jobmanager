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

        Parameters:
            rundir: str
                Path where the calculation is run from.
        """

        #To deal with bash variable paths
        if rundir[0] == '$':
            rundir = os.environ[rundir[1:]]

        with open(rundir + "/../psi4_config.json", "r") as f:
            psi4_config = json.load(f)
        success_count = 0
        psi4_utils = Psi4Utils(psi4_config)

        #Run the initial calculation from the provided molden
        success = psi4_utils.run_with_check('run_initial')
        success_count += success

        #for all other functionals, run single points from the initial calculation
        for ii, functional in enumerate(psi4_config["functional"]):
            success = psi4_utils.run_with_check('run_general', functional=functional)
            success_count += success

        print("total successful jobs : %d/ %d." %
            (success_count, len(psi4_config["functional"]) + 1))

    def loop_rescue(self, rundir='$SGE_O_WORKDIR', alphalist=[20, 15, 10, 5, 2]):
        """
        Attempts to rescue a calculation by converging a series of calculations
        using the same functional with varying HFX percentages. The first
        calculation will be functional with alphalist[0] HFX percent,
        converged from the wave function specified in wfnpath. The next
        calculation will be functional with alphalist[1] HFX percent,
        converged from the result of the previous calculation, etc.
        Should be run after a run_with_check on the same functional has
        already been run. That method will catch errors like no output/iterations,
        while this method addresses SCF nonconvergence or other errors.
        By default, increases HFX to 20% to correspond to B3LYP
        and then steps it down to 0 to match non-hybrid functionals.

        Default alphalist may perform poorly when B3LYP is not the base functional.

        Parameters:
            rundir: str
                Path where the calculation is run from.
            alphalist : list of ints
                List of the HFX percentages to use in the rescue scheme.
                The calculation will be converged using the first value with the base functional result,
                then the result of that will be used to converge the second value, etc.
        """

        #To deal with bash variable paths
        if rundir[0] == '$':
            rundir = os.environ[rundir[1:]]

        with open(rundir + "/../psi4_config.json", "r") as f:
            psi4_config = json.load(f)
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
                    print(f"HFX: {alpha}")
                    if jj == 0:
                        #For first calculation, use base functional as initial guess
                        wfn = psi4_config['wfnfile']
                    else:
                        #For subsequent calculations, use the functional converged at the previous HFX step
                        wfn = "%s/wfn.180.npy" % (functional + "-%d" %
                                                    alphalist[jj-1])
                    #If the prior calculation succeeded (or on first calculation), continue the scheme
                    if success or jj==0:
                        success = psi4_utils.rescue_with_check(functional, wfn, alpha)
                        success_count += success
                    #Rescue scheme fails if any calculation in the series fails
                    else:
                        failed = True
                        print(f'Rescue scheme has failed for functional {functional}, on HFX {alphalist[jj-1]}.')
                        break
                if not failed:
                    rescued += 1
                    print("rescued: ", functional)
            else:
                rescued += 1
                print("%s has already succeeded!" % functional)
        print("all finished!")

    def loop_derivative_jobs(self, rundir='$SGE_O_WORKDIR'):
        """
        From a directory, will converge the base functional wavefunctions of all derivative jobs.
        The first derivative job is converged from the TC molden, while the rest are
        converged from the .wfn of the first converged job.
        Should be launched from the parent directory, i.e., the folder containing
        the subfolders containing the trigger for derivative jobs.

        Parameters:
            rundir: str
                Path where the calculation is run from.
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
                                    success_count=success_count, run_func='run_initial', error_scf=True)
        os.chdir(basedir)
        # ---run other jobs using the wfn from the previous step---
        for ii, job in enumerate(jobs[1:]):
            psi4_config["wfnfile"] = derivative_utils.get_wfn_path(psi4_config, jobs, ii+1)
            os.chdir(job)
            success_count = derivative_utils.run_with_check(job=job, psi4_config=psi4_config,
                                        success_count=success_count, run_func='run_general',
                                        error_scf=True)
            os.chdir(basedir)


    def loop_hfx_jobs(self, rundir='$SGE_O_WORKDIR'):
        """
        For each specified functional, calculates the value at each of the HFX percentages specified.
        Starts each calculation from the converged wavefunction of the previous calculation.

        Runs multiple sweeps: first, a 20% HFX result is done from the base functional wfn,
        and then the program steps down and up from that calculation to the specified values.
        Then, attempts several schemes to address points that did not converge:
        (1) Start from the 0% result and step up, which allows different initial guesses for HFX<20 and retries calculations above.
        (2) Start from the 100% result and step down, which allows different initial guesses for HFX>100 and retries calculations below.

        Parameters:
            rundir: str
                Directory where the calculations are run from.
        """
        #To deal with bash variable paths
        if rundir[0] == '$':
            rundir = os.environ[rundir[1:]]

        with open(rundir + "/../psi4_config.json", "r") as f:
            psi4_config = json.load(f)
        orig_wfn_dir = psi4_config['wfnfile']
        success_count = 0
        psi4_utils = Psi4Utils(psi4_config)
        #get the base functional (called original functional here since the unmodified HFX functionals called base below)
        if "base_functional" not in psi4_config:
            #defaults to B3LYP
            orig_func = 'b3lyp'
        else:
            #remove parentheses from functional names
            functional = psi4_config["base_functional"]
            orig_func = functional.replace("(", "l-").replace(")", "-r")

        #Run the initial calculation from the provided molden
        success = psi4_utils.run_with_check('run_initial')
        success_count += success

        #Determine the amount of HFX used in the original functional
        if orig_func == 'b3lyp':
            orig_hfx = 20
        elif orig_func == 'pbe0':
            orig_hfx = 25
        else:
            raise ValueError("Base functional not yet supported for the HFX workflow.")

        #Get and sort the HFX levels desired in the calculation
        hfx_amounts = sorted(psi4_config['hfx_levels'])
        below = [alpha for alpha in hfx_amounts if alpha < orig_hfx][::-1] #reversed since you want to step down from 20
        above = [alpha for alpha in hfx_amounts if alpha > orig_hfx]

        #for all other functionals (i.e., not the original functional)
        for ii, base_functional in enumerate(psi4_config["functional"]):
            #For each functional, want to start from the original wfn
            #Note: after updating psi4_config, have to reinitialize psi4_utils to get the right parameters
            psi4_config["wfnfile"] = orig_wfn_dir
            psi4_utils = Psi4Utils(psi4_config)
            base_wfn = '' #to store the path of the result with HFX matching the original functional

            print(f'Pass 1: Starting from {orig_hfx} and stepping up/down:')
            #converge result with HFX matching the original functional
            functional = base_functional + '_hfx_' + str(orig_hfx)
            success = psi4_utils.run_with_check('run_general', functional, return_wfn=True, verbose=True, retry_scf=True)
            if success:
                #If success, want the next calculation to be run from the wfn of this calculation
                #Otherwise, run it from the last converged calculation
                base_wfn = functional.replace("(", "l-").replace(")", "-r") + '/wfn.180.npy'
                psi4_config["wfnfile"] = base_wfn
                psi4_utils = Psi4Utils(psi4_config)
                print('Wfn updated!')
            #Converge calculations below the original HFX
            for jj, alpha in enumerate(below):
                functional = base_functional + '_hfx_' + str(alpha)
                success = psi4_utils.run_with_check('run_general', functional, return_wfn=True, verbose=True, retry_scf=True)
                if success:
                    psi4_config["wfnfile"] = functional.replace("(", "l-").replace(")", "-r") + '/wfn.180.npy'
                    psi4_utils = Psi4Utils(psi4_config)
                    print('Wfn updated!')
            #Get orig_hfx wfn for other half of checks (if orig_hfx did not converge, use original functional result)
            psi4_config["wfnfile"] = base_wfn if base_wfn != '' else orig_wfn_dir
            psi4_utils = Psi4Utils(psi4_config)
            #Converge calculations above orig_hfx
            for jj, alpha in enumerate(above):
                functional = base_functional + '_hfx_' + str(alpha)
                success = psi4_utils.run_with_check('run_general', functional, return_wfn=True, verbose=True, retry_scf=True)
                if success:
                    psi4_config["wfnfile"] = functional.replace("(", "l-").replace(")", "-r") + '/wfn.180.npy'
                    psi4_utils = Psi4Utils(psi4_config)
                    print('Wfn updated!')

            print(f'Pass 2: Starting from {hfx_amounts[0]} and stepping up:')
            #Start from the original functional wavefunction if not converged already
            psi4_config["wfnfile"] = orig_wfn_dir
            psi4_utils = Psi4Utils(psi4_config)
            for jj, alpha in enumerate(hfx_amounts):
                functional = base_functional + '_hfx_' + str(alpha)
                success = psi4_utils.run_with_check('run_general', functional, return_wfn=True, verbose=True, retry_scf=True)
                if success:
                    psi4_config["wfnfile"] = functional.replace("(", "l-").replace(")", "-r") + '/wfn.180.npy'
                    psi4_utils = Psi4Utils(psi4_config)
                    print('Wfn updated!')

            print(f'Pass 3: Starting from {hfx_amounts[-1]} and stepping down:')
            #Start from the original functional wavefunction if not converged already
            psi4_config["wfnfile"] = orig_wfn_dir
            psi4_utils = Psi4Utils(psi4_config)
            for jj, alpha in enumerate(hfx_amounts[::-1]):
                functional = base_functional + '_hfx_' + str(alpha)
                success = psi4_utils.run_with_check('run_general', functional, return_wfn=True, verbose=True, retry_scf=True)
                success_count += success #only include on final pass to get accurate count
                if success:
                    psi4_config["wfnfile"] = functional.replace("(", "l-").replace(")", "-r") + '/wfn.180.npy'
                    psi4_utils = Psi4Utils(psi4_config)
                    print('Wfn updated!')

        print("total successful jobs : %d/ %d." %
            (success_count, len(psi4_config["functional"])*len(psi4_config['hfx_levels']) + 1))




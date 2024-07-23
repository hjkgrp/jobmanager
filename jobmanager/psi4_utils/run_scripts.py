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

        #Run the initial B3LYP calculation from the provided molden
        success = psi4_utils.run_with_check('run_b3lyp')
        success_count += success

        #for all other functionals, run single points from the B3LYP reference
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

        Parameters:
            rundir: str
                Path where the calculation is run from.
            alphalist : list of ints
                List of the HFX percentages to use in the rescue scheme.
                The calculation will be converged using the first value with the B3LYP result,
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
                        #For first calculation, use b3lyp as initial guess
                        wfn = "b3lyp/wfn.180.npy"
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
        From a directory, will converge the B3LYP wavefunctions of all derivative jobs.
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


    def loop_hfx_jobs(self, rundir='$SGE_O_WORKDIR'):
        """
        For each specified functional, calculates the value at each of the HFX percentages specified.
        Starts each calculation from the converged wavefunction of the previous calculation.

        Parameters:
            rundir: str
                Directory where the calculations are run from.
        """
        #To deal with bash variable paths
        if rundir[0] == '$':
            rundir = os.environ[rundir[1:]]

        with open(rundir + "/../psi4_config.json", "r") as f:
            psi4_config = json.load(f)
        b3lyp_wfn_dir = psi4_config['wfnfile']
        success_count = 0
        psi4_utils = Psi4Utils(psi4_config)
        
        #Run the initial B3LYP calculation from the provided molden
        success = psi4_utils.run_with_check('run_b3lyp')
        success_count += success

        #Get and sort the HFX levels desired in the calculation
        #Unsorted to allow different orders if desired
        hfx_amounts = psi4_config['hfx_levels']

        #for all other functionals
        for ii, base_functional in enumerate(psi4_config["functional"]):
            #For each functional, want to start from the B3LYP wfn
            #Note: after updating psi4_config, have to reinitialize
            #psi4_utils to get the right parameters
            psi4_config["wfnfile"] = b3lyp_wfn_dir
            psi4_utils = Psi4Utils(psi4_config)
            for jj, alpha in enumerate(hfx_amounts):
                #functional you want to run is base functional plus the HFX level
                functional = base_functional + '_hfx_' + str(alpha)
                success = psi4_utils.run_with_check('run_general', functional)
                success_count += success
                if success:
                    #If success, want the next calculation to be run from the wfn of this calculation
                    #Otherwise, run it from the last converged calculation
                    psi4_config["wfnfile"] = functional.replace("(", "l-").replace(")", "-r") + '/wfn.180.npy'
                    psi4_utils = Psi4Utils(psi4_config)
                    print('Wfn updated!')

        print("total successful jobs : %d/ %d." %
            (success_count, len(psi4_config["functional"])*len(psi4_config['hfx_levels']) + 1))

        


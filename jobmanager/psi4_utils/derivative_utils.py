import os
import operator
from typing import List
from jobmanager.psi4_utils.psi4_utils import Psi4Utils

class DerivativeUtils():
    """
    Contains functions useful for running derivative jobs.
    """

    def __init__(self):
        pass

    def get_subfolders(self, path: str = "./") -> list:
        """
        Gets all subfolders in the specified path.

        Parameters:
            path: str
                Path to get subfolders of.

        Outputs:
            folders: list
                List of subfolders in path.
        """
        files = os.listdir(path)
        folders = []
        for file in files:
            if os.path.isdir(path + "/" + file):
                folders.append(file)
        return folders

    def sanity_check(self, folders: list, trigger: str = "_derivNo_") -> dict:
        """
        Returns a dictionary mapping folders of derivative jobs to their job number.
        Raises errors if the trigger for derivative jobs does not appear or if there
        are multiple base jobs.

        Parameters:
            folders: list
                List of folder names.
            trigger: str
                Substring to be used to denote derivative jobs.
        Outputs:
            jobs: dict
                Maps folder names to which job number they are.
        """
        trigger_appear = False
        basename: List[str] = list()
        folders_ignored = list()
        jobs = dict()
        for folder in folders:
            if trigger in folder:
                trigger_appear = True
                base = folder.split(trigger)[0]
                # print(base, folder)
                if not len(basename):
                    first_base = base
                    basename.append(base)
                if base not in basename:
                    basename.append(base)
                if base == first_base:
                    jobs[folder] = int(folder.split(trigger)[-1])
            else:
                folders_ignored.append(folder)
        if not trigger_appear:
            raise KeyError("<trigger>: %s not appeared in any folders." % trigger)
        if len(basename) > 1:
            raise ValueError("multiple base jobs occur: %s" % basename)
        if len(folders_ignored):
            print("Warning: The following folders are ignored: %s" % str(folders_ignored))
        return jobs

    def derivative_tree(self, path: str = "./", trigger: str = "_derivNo_") -> list:
        '''
        Make a sequence structure for derivative jobs in a path (specified as path).
        Derivative jobs are denoted by containing the substring contained in trigger.

        Parameters:
            folders: list
                List of folder names.
            trigger: str
                Substring to be used to denote derivative jobs.
        Outputs:
            list
                List of folder names in order of when they should be calculated.
        '''
        folders = self.get_subfolders(path)
        jobs = self.sanity_check(folders, trigger=trigger)
        #sort by values, which are the derivative numbers
        jobs = dict(sorted(jobs.items(), key=operator.itemgetter(1)))
        return list(jobs.keys())

    def get_wfn_path(self, psi4_config, jobs, ii):
        """
        Gets the path of the .wfn file from the ii-th item in jobs.

        Assumes calculation run from the parent directory.

        Parameters:
            psi4_config:
                Loaded JSON file giving the settings of the calculation.
            jobs:
                list of folder names that mark jobs.
            ii: int
                index of job one wants to get the wfn from.
        Outputs:
            str
                Path of the wfn file in the specified job.
        """
        assert ii > 0
        #get the base functional
        if "base_functional" not in psi4_config:
            #defaults to B3LYP
            base_func = 'b3lyp'
        else:
            #remove parentheses from functional names
            functional = psi4_config["base_functional"]
            base_func = functional.replace("(", "l-").replace(")", "-r")

        return './' + jobs[ii - 1] + "/" + base_func + "/wfn.180.npy"

    def run_with_check(self, job: str, psi4_config: dict,
                       run_func: str, success_count: int,
                       error_scf: bool = True):
        """
        Checks if the initial calculation is converged for the specified job.
        If not converged, will resubmit the calculation.

        Parameters:
            job
                Gives the name of the folder where the job is being run.
            psi4_config
                Loaded JSON file giving the settings of the calculation.
            run_func
                Either run_initial or run_general, depending on what should be run
            success_count
                Current number of successes.
            error_scf
                Whether one wants to continue doing derivative jobs after a SCF error occurs.
        """
        print("====running for====: ", job)
        success = False
        psi4_utils = Psi4Utils(psi4_config)

        #get the base functional
        if "base_functional" not in psi4_config:
            #defaults to B3LYP
            base_func = 'b3lyp'
        else:
            #remove parentheses from functional names
            functional = psi4_config["base_functional"]
            base_func = functional.replace("(", "l-").replace(")", "-r")

        if run_func == 'run_initial':
            run_function = psi4_utils.run_initial # type: ignore
        elif run_func == 'run_general':
            run_function = psi4_utils.run_general # type: ignore

        #If the initial folder does not exist, run calculation
        if not os.path.isdir(base_func):
            success = run_function(rundir='./', return_wfn=True)
            print("success: ", success)
            if success:
                success_count += 1
        else:
            #Check if initial calculation converged
            print("folder exists.")
            resubed = False
            #need to resubmit if no output or if no iterations in output
            functional = base_func
            if not os.path.isfile(functional + "/output.dat"):
                resubed = True
            else:
                with open('./' + functional + "/output.dat", "r") as fo:
                    txt = "".join(fo.readlines())
                if "==> Iterations <==" not in txt:
                    resubed = True
            #Resubmit initial calculation
            if resubed:
                print("previous errored out. resubmitting...")
                success = run_function(rundir='./', return_wfn=True)
                print("success: ", success)
                if success:
                    success_count += 1
            else:
                #If the output file exists and iterations were reached, check for error messages
                with open(functional + "/output.dat", "r") as fo:
                    txt = "".join(fo.readlines())
                if 'PsiException: Could not converge SCF iterations' not in txt and os.path.isfile(functional + "/wfn.180.npy"):
                    #Success if no SCF error and a wfn file was written
                    print("success: ", True)
                    success = True
                    success_count += 1

        if not success and error_scf:
            #If the initial calculation fails, raise an error since all subsequent calculations will not work
            raise ValueError(
                "Failed on the job: %s. Other derivative jobs won't run." % job)
        return success_count

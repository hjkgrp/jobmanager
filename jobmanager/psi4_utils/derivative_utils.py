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

    def get_wfn_path(self, jobs, ii):
        """
        Gets the path of the .wfn file from the ii-th item in jobs.

        Assumes calculation run from the parent directory.

        Parameters:
            jobs:
                list of folder names that mark jobs.
            ii: int
                index of job one wants to get the wfn from.
        Outputs:
            str
                Path of the wfn file in the specified job.
        """
        assert ii > 0
        return './' + jobs[ii - 1] + "/b3lyp/wfn.180.npy"

    def run_with_check(self, job: str, psi4_config: dict,
                       run_func: str, success_count: int,
                       error_scf: bool = True):
        """
        Checks if a B3LYP calculation is converged for the specified job.
        If not converged, will resubmit the calculation.

        Parameters:
            job
                Gives the name of the folder where the job is being run.
            psi4_config
                Loaded JSON file giving the settings of the calculation.
            run_func
                Either run_b3lyp or run_general, depending on what should be run
            success_count
                Current number of successes.
            error_scf
                Whether one wants to continue doing derivative jobs after a SCF error occurs.
        """
        print("====running for====: ", job)
        success = False
        psi4_utils = Psi4Utils(psi4_config)

        if run_func == 'run_b3lyp':
            run_function = psi4_utils.run_b3lyp # type: ignore
        elif run_func == 'run_general':
            run_function = psi4_utils.run_general # type: ignore

        #If the B3LYP folder does not exist, run from TC molden
        if not os.path.isdir("b3lyp"):
            success = run_function(rundir='./', return_wfn=True)
            print("success: ", success)
            if success:
                success_count += 1
        else:
            #Check if B3LYP calculation converged
            print("folder exists.")
            resubed = False
            #need to resubmit if no output or if no iterations in output
            functional = "b3lyp"
            if not os.path.isfile(functional + "/output.dat"):
                resubed = True
            else:
                with open('./' + functional + "/output.dat", "r") as fo:
                    txt = "".join(fo.readlines())
                if "==> Iterations <==" not in txt:
                    resubed = True
            #Resubmit B3LYP calculation
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
            #If the B3LYP calculation fails, raise an error since all subsequent calculations will not work
            raise ValueError(
                "Failed on the job: %s. Other derivative jobs won't run." % job)
        return success_count

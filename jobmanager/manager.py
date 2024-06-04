import os
import subprocess
from collections import defaultdict

def find_calcs(dirpath, extension='.xyz', original=None):
        """Find calculations that need to be run by extension.
        Based on os.walk https://github.com/python/cpython/blob/a372a7d/Lib/os.py#L344

        Parameters
        ----------
            dirpath : str
                The name of the parent/top directory.
            topdown : bool
                Whether or not to search via topdown or bottom up order.

        Returns
        -------
            gen : generator
                Generator which returns the next calculation to be done.
        """

        walk_dirs = []
        calc_paths = []
        for entry in os.scandir(dirpath):
            # do not check in hidden directories or scratch directories
            if entry.is_dir() and not entry.name.startswith('.') and not entry.name.startswith('scr'):
                walk_dirs.append(entry.path)
            elif entry.name.endswith(extension):
                calc_paths.append(entry.path)
        for new_path in walk_dirs:
            if original:
                yield from find_calcs(new_path, original=original)
            else:
                yield from find_calcs(new_path, original=dirpath)
        for path in calc_paths:
            if original:
                yield path[len(original)+1:]
            else:
                yield path[len(dirpath)+1:]


class Manager:
    """
    Carry out the steps between generating inputs & getting results
        - Manage Job instances and interfaces with queue system
        - Recovery/resubmission attempts
    """
    def __init__(self,directory):
        self.directory = directory

    def find_jobs_to_submit(self):
        """ Find jobs to submit
        Returns
        -------
            jobscripts : dict
            { 'path/to/jobscript': { 'input_file': 'path/to/input', 'xyz_file': 'path/to/xyz'} }

        """
        jobscripts = defaultdict(dict)
        directory = self.directory
        for jobscript in find_calcs(directory, extension='_jobscript'):
            in_working_dir = os.getcwd()
            os.chdir(directory)
            jobscript_abs = os.path.abspath(jobscript)
            inp_file = jobscript.rsplit('_jobscript')[0]+'.in'
            inp_file_abs = os.path.abspath(inp_file)
            coords = inp_file.rsplit('.')[0]+'.xyz'
            coords_abs = os.path.abspath(coords)
            if os.path.exists(inp_file_abs):
                jobscripts[jobscript_abs]["input_file"] = inp_file_abs
            if os.path.exists(coords_abs):
                jobscripts[jobscript_abs]["xyz_file"] = coords_abs
        os.chdir(in_working_dir)
        return jobscripts

    def check_finished(self):
        """
        Check if jobs are finished, terminated, running, or didn't run at all

        Returns
        -------
            no_run : list
                list of jobs that didn't run
            terminated : list
                list of jobs that terminated / failed
            finished : list
                list of jobs that finished
            running : list
                list of jobs that are running at the moment
        """
        no_run = []
        terminated = []
        finished = []
        running = []
        jobscripts = self.find_jobs_to_submit()
        for job in jobscripts:
            outfile = job.rsplit('_jobscript')[0]+'.out'
            if os.path.exists(outfile):
                with open(outfile, 'r', errors='ignore') as f:
                    outfile_lines = f.readlines()
                if 'Job finished' in outfile_lines[-1]:
                    finished.append(job)
                elif 'Job terminated' in outfile_lines[-1]:
                    terminated.append(job)
                else:
                    running.append(job)
            else:
                no_run.append(job)
        return no_run, terminated, finished, running

    def submit(self):
        """
        Submits the jobs that have both input and xyz files
        """
        jobscript_dict = self.find_jobs_to_submit()
        for jobscript in jobscript_dict.keys():
            if "input_file" in jobscript_dict[jobscript] and "xyz_file" in jobscript_dict[jobscript]:
                print('- - - - Submitting the job - - - -')
                subprocess.call(f'qsub {jobscript}', shell=True)

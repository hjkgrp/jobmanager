
import os
class Job:
    """
    Organize calculation(s) to run (independent of queue system)
    - Run history
    - Dependent calculation structure
    - Batch jobs
    Add properties to Job class that correspond to jobscript properties
    (e.g., queue system (SGE or SLURM), max run time, memory,
    number of cores to run job with, etc)
    """
    def __init__(self):
        self.queue_system = None
        self.runtime = None #h_rss
        self.gpu = None # smp
        self.maxtime = None # h_rss
        self.mem_per_cpu = None # mem-per-cpu

        self.cpu = None # cpus-per-task
        self.nodes = None # nodes
        self.ntasks = None # ntasks

    def read_jobscript(self, jobscript):
        """
        reads jobscript file and sets data as a property of class

            Parameters
            ----------
                jobscript : str
                    Absolute path to the jobscript.

        """
        #read the lines
        with open(jobscript, 'r', errors='ignore') as f:
            jobsc_lines = f.readlines()

        ref_sge = ['h_rt','h_rss']
        req_sge_dict = {}
        ref_slurm = ['nodes', 'ntasks', 'cpus-per-task', 'mem-per-cpu', 'time']
        req_slurm_dict = {}

        if '#$' in ''.join(jobsc_lines): #check if SGE:
            self.queue_system = 'SGE'
            for line in jobsc_lines:
                line = line.strip()
                if not line.startswith('#'): continue
                line = line.split(' ')
                if line[-1].split('=')[0] in ref_sge:
                    key,val = line[-1].split('=')
                    req_sge_dict[key] = val
                elif 'smp' in line:
                    self.gpu = line[-1]
            self.runtime = req_sge_dict['h_rt']
            self.maxmemory = req_sge_dict['h_rss']

        elif '#!' in ''.join(jobsc_lines): #check if SLURM
            self.queue_system = 'SLURM'
            for line in jobsc_lines:
                line = line.strip()
                if not line.startswith('#SBATCH') or not '--' in line: continue
                line = line.split(' ')
                if len(line) != 2: continue
                key,val = line[1][2:].split('=')
                if key in ref_slurm:
                    req_slurm_dict[key] = val

            self.nodes = req_slurm_dict['nodes']
            self.ntasks = req_slurm_dict['ntasks']
            self.cpu = req_slurm_dict['cpus-per-task']
            self.mem_per_cpu = req_slurm_dict['mem-per-cpu']
            self.runtime = req_slurm_dict['time']

        else:
            print('- - - given queue system is not supported - - -')

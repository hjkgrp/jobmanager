import os
from importlib_resources import files as resource_files
import subprocess

class RunUtils():
    """
    Contains functions useful for running Psi4 calculations.
    """

    def __init__(self):
        pass

    def ensure_dir(self, dirpath):
        """
        Checks if a directory exists. If not, makes that directory.

        Parameters:
            dirpath: str
                Path of the directory that you want to ensure.
        """
        if not os.path.isdir(dirpath):
            os.makedirs(dirpath)

    def check_success(self, path='./', output_name='output.dat'):
        """
        Checks if a Psi4 calculation is completed.
        Assumes it is being called in the same directory as the Psi4 output file,
        and that the output file is named output.dat.

        Parameters:
            path: str
                Path where the calculation result is stored.
            output_name: str
                Name of the output file.

        Returns:
            success: bool
                Whether or not the calculation succeeded.
        """
        success = False
        #only look at the text from the latest computation
        #Psi4 will call tstart() at the start of each computation
        txt = ""
        with open(path + '/' + output_name, "r") as fo:
            for idx, line in enumerate(reversed(fo.readlines())):
                txt += line
                if 'tstart()' in line:
                    break
        #Psi4 will write Computation Completed once the calculation finishes.
        if 'Computation Completed' in txt:
            success = True
        return success

    def write_jobscript(self, psi4_config):
        """
        From a psi4_config JSON file, writes a jobscript for the appropriate cluster.
        Memory should be specified in MB.

        Parameters:
            psi4_config: str
                Path to the psi4_config.json file containing job details.
        """
        if "cluster" not in psi4_config or psi4_config["cluster"] == "gibraltar":
            mem = int(psi4_config['memory'].split(" ")[0])/1000
            with open("./jobscript.sh", "w") as fo:
                fo.write("#$ -S /bin/bash\n")
                fo.write("#$ -N psi4_dft\n")
                fo.write("#$ -R y\n")
                fo.write("#$ -cwd\n")
                fo.write("#$ -l h_rt=240:00:00\n")
                fo.write("#$ -l h_rss=%dG\n" % (mem))
                fo.write("#$ -q cpus\n")
                fo.write("#$ -l cpus=1\n")
                fo.write("#$ -pe smp %d\n" % psi4_config['num_threads'])
                fo.write("# -fin *\n")

                fo.write(f"source {psi4_config['bashrc']}\n")
                fo.write(f"conda activate {psi4_config['conda_env']}\n")
                fo.write("export PSI_SCRATCH='./'\n")
                fo.write("echo 'psi4 scr: ' $PSI_SCRATCH\n")

                if "trigger" in psi4_config:
                    fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_derivative_jobs()'  > $SGE_O_WORKDIR/deriv_nohup1.out 2> $SGE_O_WORKDIR/deriv_nohup1.err\n")
                elif "hfx_levels" in psi4_config:
                    fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_hfx_jobs()'  > $SGE_O_WORKDIR/nohup1.out 2> $SGE_O_WORKDIR/nohup1.err\n")
                else:
                    fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_run()'  > $SGE_O_WORKDIR/nohup1.out 2> $SGE_O_WORKDIR/nohup1.err\n")
                    fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_run()'  > $SGE_O_WORKDIR/nohup2.out 2> $SGE_O_WORKDIR/nohup1.err\n")
                    fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_run()'  > $SGE_O_WORKDIR/nohup3.out 2> $SGE_O_WORKDIR/nohup1.err\n")
                    if "hfx_rescue" in psi4_config and psi4_config["hfx_rescue"]:
                        fo.write("echo rescuing...\n")
                        fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_rescue()' > $SGE_O_WORKDIR/rescue_nohup1.out\n")
                        fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_rescue()' > $SGE_O_WORKDIR/rescue_nohup2.out\n")
                        fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_rescue()' > $SGE_O_WORKDIR/rescue_nohup3.out\n")
                fo.write("cp -rf * $SGE_O_WORKDIR\n")
                fo.write("sleep 30\n")
        elif psi4_config["cluster"] == "supercloud":
            mem = int(psi4_config['memory'].split(" ")[0])/1000
            with open("./jobscript.sh", "w") as fo:
                fo.write("#!/bin/bash\n")
                fo.write("#SBATCH --job-name=psi4_multiDFA\n")
                fo.write("#SBATCH --nodes=1\n")
                fo.write("#SBATCH --time=96:00:00\n")
                fo.write("#SBATCH --ntasks-per-node=%d\n" % psi4_config['num_threads'])
                if "queue" in psi4_config and psi4_config["queue"] == "normal":
                    fo.write("#SBATCH --partition=normal\n")
                fo.write("#SBATCH --mem=%dG\n\n" % mem)

                fo.write("source /etc/profile\n")
                fo.write("source ~/.profile\n")
                fo.write("source ~/.bashrc\n")
                fo.write(f"conda activate {psi4_config['conda_env']}\n")
                fo.write("export PSI_SCRATCH='./'\n\n")
                fo.write("subdir=$PWD\n")
                fo.write("echo subdir: $subdir\n")
                fo.write("echo tmpdir: $TMPDIR\n")
                fo.write("cp -rf * $TMPDIR\n")
                fo.write("cd $TMPDIR\n\n")

                if "trigger" in psi4_config:
                    fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_derivative_jobs()'  > $SGE_O_WORKDIR/deriv_nohup1.out 2> $SGE_O_WORKDIR/deriv_nohup1.err\n")
                    fo.write("rm */*/psi.* */*/dfh.* */*-*/*.npy */b3lyp/*.molden */b3lyp/*1step*\n")
                else:
                    fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_run()'  > $SGE_O_WORKDIR/nohup1.out 2> $SGE_O_WORKDIR/nohup1.err\n")
                    fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_run()'  > $SGE_O_WORKDIR/nohup2.out 2> $SGE_O_WORKDIR/nohup1.err\n")
                    fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_run()'  > $SGE_O_WORKDIR/nohup3.out 2> $SGE_O_WORKDIR/nohup1.err\n")
                    if "hfx_rescue" in psi4_config and psi4_config["hfx_rescue"]:
                        fo.write("echo rescuing...\n")
                        fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_rescue()' > $SGE_O_WORKDIR/rescue_nohup1.out\n")
                        fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_rescue()' > $SGE_O_WORKDIR/rescue_nohup2.out\n")
                        fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_rescue()' > $SGE_O_WORKDIR/rescue_nohup3.out\n")
                    fo.write("rm */psi.* */dfh.* *-*/*.npy b3lyp/*.molden b3lyp/*1step*\n")
                fo.write("cp -rf * $subdir\n")
        elif psi4_config["cluster"] == "expanse":
            mem = int(psi4_config['memory'].split(" ")[0])/psi4_config['num_threads']/1000
            with open("./jobscript.sh", "w") as fo:
                fo.write("#!/bin/sh\n")
                fo.write("#SBATCH -A mit136\n")
                fo.write("#SBATCH --job-name=psi4_multiDFA\n")
                fo.write("#SBATCH --partition=shared\n")
                fo.write("#SBATCH -t 48:00:00\n")
                fo.write("#SBATCH --nodes=1\n")
                fo.write("#SBATCH --ntasks-per-node=16\n")
                fo.write("#SBATCH --error=job.%J.err\n")
                fo.write("#SBATCH --output=job.%J.out\n")
                fo.write("#SBATCH --export=ALL\n")
                fo.write("#SBATCH --mem=64G\n")

                fo.write("source /home/crduan/.bashrc\n")
                fo.write("conda activate mols_psi4\n")
                fo.write("export PSI_SCRATCH='./'\n")
                fo.write("echo 'psi4 scr: ' $PSI_SCRATCH\n")

                if "trigger" in psi4_config:
                    fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_derivative_jobs()'  > $SGE_O_WORKDIR/deriv_nohup1.out 2> $SGE_O_WORKDIR/deriv_nohup1.err\n")
                else:
                    fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_run()'  > $SGE_O_WORKDIR/nohup1.out 2> $SGE_O_WORKDIR/nohup1.err\n")
                    fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_run()'  > $SGE_O_WORKDIR/nohup2.out 2> $SGE_O_WORKDIR/nohup1.err\n")
                    fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_run()'  > $SGE_O_WORKDIR/nohup3.out 2> $SGE_O_WORKDIR/nohup1.err\n")
                    if "hfx_rescue" in psi4_config and psi4_config["hfx_rescue"]:
                        fo.write("echo rescuing...\n")
                        fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_rescue()' > $SGE_O_WORKDIR/rescue_nohup1.out\n")
                        fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_rescue()' > $SGE_O_WORKDIR/rescue_nohup2.out\n")
                        fo.write("python -c 'from jobmanager.psi4_utils.run_scripts import RunScripts; RunScripts().loop_rescue()' > $SGE_O_WORKDIR/rescue_nohup3.out\n")
        elif psi4_config["cluster"] == "mustang":
            #TODO: Update
            with open("./jobscript.sh", "w") as fo:
                fo.write("#!/bin/bash\n")
                fo.write("#PBS -N psi4_multiDFA\n")
                fo.write("#PBS -A ONRDC42143511\n")
                fo.write("#PBS -o psi4.out\n")
                fo.write("#PBS -e psi4.err\n")
                fo.write("#PBS -l select=1:ncpus=48:mpiprocs=48\n")
                fo.write("#PBS -l walltime=6:00:00\n")
                fo.write("#PBS -q standard\n")
                fo.write("#PBS -j oe\n")
                fo.write("#PBE -V\n")

                fo.write("source $HOME/.personal.bashrc\n")
                fo.write("export PSI_SCRATCH='./'\n")
                fo.write("rundir=$PBS_O_WORKDIR\n")
                fo.write("cd $rundir\n")
                fo.write("homekey='home'\n")
                fo.write("homedir=${rundir/work1/$homekey}\n")
                fo.write("echo homedir: $homedir\n")
                fo.write("echo rundir: $rundir\n")
                fo.write("mkdir -p $homedir\n")
                fo.write("python -u loop_run.py  > nohup1.out\n")
                fo.write("zip -r outdat.zip */*.dat > zip.out\n")
                fo.write("cp outdat.zip $homedir\n")
                fo.write("python -u loop_run.py  > nohup2.out\n")
                fo.write("zip -r outdat.zip */*.dat > zip.out\n")
                fo.write("cp outdat.zip $homedir\n")
                fo.write("python -u loop_run.py  > nohup3.out\n")
                fo.write("zip -r outdat.zip */*.dat > zip.out\n")
                fo.write("cp outdat.zip $homedir\n")
                if "hfx_rescue" in psi4_config and psi4_config["hfx_rescue"]:
                    fo.write("echo rescuing...\n")
                    fo.write("python -u loop_rescue.py > rescue_nohup1.out\n")
                    fo.write("zip -r outdat.zip */*.dat > zip.out\n")
                    fo.write("cp outdat.zip $homedir\n")
                    fo.write("python -u loop_rescue.py > rescue_nohup2.out\n")
                    fo.write("zip -r outdat.zip */*.dat > zip.out\n")
                    fo.write("cp outdat.zip $homedir\n")
                    fo.write("python -u loop_rescue.py > rescue_nohup3.out\n")
                    fo.write("zip -r outdat.zip */*.dat > zip.out\n")
                    fo.write("cp outdat.zip $homedir\n")
                fo.write("echo all done.\n")

    def run_bash(self, cmd):
        """
        Runs the specified bash command from the current directory.

        Parameters:
            cmd: str
                Bash command you want to run.
        """
        subprocess.call(cmd, shell=True)
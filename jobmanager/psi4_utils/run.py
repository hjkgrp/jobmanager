#Functionality being moved to run_utils.py, psi4_utils.py

import psi4
import os
import numpy as np
import shutil
import subprocess
from jobmanager.io.molden import load_molden
import json
from jobmanager.psi4_utils.molden2psi4wfn import tcmolden2psi4wfn_ao_mapping
from jobmanager.psi4_utils.molden2psi4wfn_spherical import tcmolden2psi4wfn_ao_mapping_spherical
from importlib_resources import files as resource_files

#moved to psi4_utils
def lacvps(mol, role):
    '''
    Define the LACVP* basis set.
    lanl2dz for metals and 6-31g* for others.
    '''
    basstrings = {}
    mol.set_basis_all_atoms("6-31g*", role=role)
    mol.set_basis_by_symbol("fe", "lanl2dz", role=role)
    mol.set_basis_by_symbol("co", "lanl2dz", role=role)
    mol.set_basis_by_symbol("cr", "lanl2dz", role=role)
    mol.set_basis_by_symbol("mn", "lanl2dz", role=role)
    mol.set_basis_by_symbol("mo", "lanl2dz", role=role)
    mol.set_basis_by_symbol("tc", "lanl2dz", role=role)
    mol.set_basis_by_symbol("ru", "lanl2dz", role=role)
    mol.set_basis_by_symbol("rh", "lanl2dz", role=role)
    mol.set_basis_by_symbol("I", "lanl2dz", role=role)
    mol.set_basis_by_symbol("Br", "lanl2dz", role=role)
    mol.set_basis_by_symbol("hf", "lanl2dz", role=role)
    return basstrings

#moved to psi4_utils
def get_molecule(xyzfile, charge, spin, sym='c1'):
    '''
    Assemble a molecule object from xyzfile, charge and spin.
    '''
    wholetext = "%s %s\n" % (charge, spin)
    if os.path.isfile(xyzfile):
        with open(xyzfile, "r") as fo:
            natoms = int(fo.readline().split()[0])
            fo.readline()
            for ii in range(natoms):
                wholetext += fo.readline()
    wholetext += "\nsymmetry %s\nnoreorient\nnocom\n" % sym
    mol = psi4.geometry("""%s""" % wholetext)
    return mol

#moved to psi4_utils
def setup_dft_parameters(psi4_config):
    psi4.set_memory(psi4_config["memory"])
    psi4.set_num_threads(psi4_config["num_threads"])
    if psi4_config["basis"] == "lacvps":
        psi4.qcdb.libmintsbasisset.basishorde['LACVPS'] = lacvps
        psi4.set_options({"puream": False})
    elif psi4_config["basis"] == "def2-tzvp":
        psi4.set_options({"puream": True})
    else:
        psi4.set_options({"puream": True})
        # pass
        # raise ValueError("Only lacvps is supported!")
    psi4.set_options({
        'reference': psi4_config["ref"],
        "DF_SCF_GUESS": False,
        "scf_type": "df",
        "dft_pruning_scheme": "robust",
        "basis": psi4_config["basis"],
        "DFT_BASIS_TOLERANCE": 1e-10,
        "INTS_TOLERANCE": 1e-10,
        "PRINT_MOS": False,
        "dft_spherical_points": 590,
        "dft_radial_points": 99,
        # "SOSCF": True,
        # "DAMPING_PERCENTAGE": 20,
        # "BASIS_GUESS": True,
        "guess": "read", })

#moved to run_utils
def ensure_dir(dirpath):
    if not os.path.isdir(dirpath):
        os.makedirs(dirpath)

#moved to run_utils
def check_sucess():
    success = False
    with open("output.dat", "r") as fo:
        txt = "".join(fo.readlines())
    if 'Computation Completed' in txt:
        success = True
    return success

#redundant with get_hfx_functional, removed in psi4_utils
def b3lyp_hfx():
    b3lyp_d = {}
    for hfx in [0, 5, 10, 15, 20, 25, 30, 35, 40]:
        hfx_func = {
            "name": "hfx_func",
            "x_functionals": {
                "GGA_X_B88": {"alpha": 0.9*(1-hfx*0.01)},
                "LDA_X": {"alpha": 0.1*(1-hfx*0.01)}
                    },
            "x_hf": {"alpha": hfx*0.01},
            "c_functionals": {
                "GGA_C_LYP": {"alpha": 0.81},
                "LDA_C_VWN_RPA": {"alpha": 0.19}
            }
        }
        b3lyp_d["b3lyp_" + str(hfx)] = hfx_func
    return b3lyp_d

#moved to psi4_utils, edited there so that not running in the b3lyp folder, but rather in the parent folder
def run_b3lyp(psi4_config, rundir="./b3lyp", return_wfn=True):
    b3lyp_d = b3lyp_hfx()
    psi4_scr = './'
    filename = "output"
    basedir = os.getcwd()
    with open(psi4_config["charge-spin-info"], "r") as f:
        d = json.load(f)
    psi4_config.update(d)
    ensure_dir(rundir)
    shutil.copyfile(psi4_config["xyzfile"], os.path.join(rundir, psi4_config["xyzfile"]))
    sym = 'c1' if 'sym' not in psi4_config else psi4_config['sym']
    mol = get_molecule(psi4_config["xyzfile"], psi4_config["charge"], psi4_config["spin"], sym)
    setup_dft_parameters(psi4_config)
    pid = str(os.getpid())
    if os.path.isfile(psi4_config["moldenfile"]):
        shutil.copyfile(psi4_config["moldenfile"], rundir + '/' + psi4_config["moldenfile"])
        os.chdir(rundir)
        psi4.core.set_output_file(filename + '.dat', False)
        # 1-step SCF
        psi4.set_options({
            "maxiter": 5,
            "D_CONVERGENCE": 1e5,
            "E_CONVERGENCE": 1e5,
            "fail_on_maxiter": False})
        if psi4_config["basis"] == "def2-tzvp":
            psi4.set_options({"basis": "def2-sv(p)"})
        if "b3lyp_hfx" in psi4_config and psi4_config["b3lyp_hfx"] != 20:
            print("customized b3lyp with different HFX: ", psi4_config["b3lyp_hfx"])
            functional = "b3lyp_" + str(psi4_config["b3lyp_hfx"])
            e, wfn = psi4.energy("scf", dft_functional=b3lyp_d[functional],  molecule=mol, return_wfn=True)
        else:
            e, wfn = psi4.energy('b3lyp', molecule=mol, return_wfn=True)
        wfn.to_file("wfn-1step.180")
        # Get converged WFN
        d_molden = load_molden(psi4_config["moldenfile"])
        restricted = True if any(x in psi4_config["ref"] for x in ["r", "R"]) else False
        if not psi4_config["basis"] == "def2-tzvp":
            Ca, Cb, mapping = tcmolden2psi4wfn_ao_mapping(d_molden, restricted=restricted)
        else:
            Ca, Cb, mapping = tcmolden2psi4wfn_ao_mapping_spherical(d_molden, restricted=restricted)
        wfn_minimal_np = np.load("wfn-1step.180.npy", allow_pickle=True)
        wfn_minimal_np[()]['matrix']["Ca"] = Ca
        if not restricted:
            wfn_minimal_np[()]['matrix']["Cb"] = Cb
        else:
            wfn_minimal_np[()]['matrix']["Cb"] = Ca
        np.save("wfn-1step-tc.180.npy", wfn_minimal_np)
        # Copy wfn file to the right place with a right name
        pid = str(os.getpid())
        targetfile = psi4_scr + filename + '.default.' + pid + '.180.npy'
        shutil.copyfile("wfn-1step-tc.180.npy", targetfile)
        # Final scf---
        psi4.set_options({
            "maxiter": 50,
            "D_CONVERGENCE": 3e-5,
            "E_CONVERGENCE": 3e-5,
            "fail_on_maxiter": True})
    else:
        os.chdir(rundir)
        psi4.core.set_output_file(filename + '.dat', False)
        print("Warning: no Molden file is used to initialize this calculation!")
        psi4.set_options({
            "maxiter": 250 if "maxiter" not in psi4_config else psi4_config["maxiter"],
            # "guess": "GWH",
            "D_CONVERGENCE": 3e-5,
            "E_CONVERGENCE": 3e-5,
            "fail_on_maxiter": True})
        if psi4_config["basis"] == "def2-tzvp":
            psi4.set_options({"basis": "def2-sv(p)"})
    sucess = False
    try:
        if "b3lyp_hfx" in psi4_config and psi4_config["b3lyp_hfx"] != 20:
            print("customized b3lyp with different HFX: ", psi4_config["b3lyp_hfx"])
            functional = "b3lyp_" + str(psi4_config["b3lyp_hfx"])
            e, wfn = psi4.energy("scf", dft_functional=b3lyp_d[functional],  molecule=mol, return_wfn=True)
        else:
            e, wfn = psi4.energy('b3lyp', molecule=mol, return_wfn=True)
        wfn.to_file("wfn.180")
        sucess = True
    except:
        print("This calculation does not converge.")
    if psi4_config["basis"] == "def2-tzvp" and sucess:
        psi4.set_options({"basis": "def2-tzvp", "maxiter": 200 if "maxiter" not in psi4_config else psi4_config["maxiter"]})
        try:
            if "b3lyp_hfx" in psi4_config and psi4_config["b3lyp_hfx"] != 20:
                print("customized b3lyp with different HFX: ", psi4_config["b3lyp_hfx"])
                functional = "b3lyp_" + str(psi4_config["b3lyp_hfx"])
                e, wfn = psi4.energy("scf", dft_functional=b3lyp_d[functional],  molecule=mol, return_wfn=True)
            else:
                e, wfn = psi4.energy('b3lyp', molecule=mol, return_wfn=True)
            wfn.to_file("wfn.180")
        except:
            print("This calculation does not converge.")
    success = check_sucess()
    for filename in os.listdir("./"):
        if ("psi." in filename) or ("default" in filename):
            print("removing: :", filename)
            os.remove(filename)
    os.chdir(basedir)
    return success

#moved to psi4_utils, edited there to run in the parent folder instead of in the functional subfolder
def run_general(psi4_config, functional="b3lyp", return_wfn=False):
    b3lyp_d = b3lyp_hfx()
    psi4_scr = './'
    filename = "output"
    basedir = os.getcwd()
    rundir = "./" + functional.replace("(", "l-").replace(")", "-r")
    with open(psi4_config["charge-spin-info"], "r") as f:
        d = json.load(f)
    psi4_config.update(d)
    ensure_dir(rundir)
    shutil.copyfile(psi4_config["xyzfile"], functional.replace("(", "l-").replace(")", "-r") + '/' + psi4_config["xyzfile"])
    os.chdir(rundir)
    psi4.core.set_output_file(filename + '.dat', False)
    sym = 'c1' if 'sym' not in psi4_config else psi4_config['sym']
    mol = get_molecule(psi4_config["xyzfile"], psi4_config["charge"], psi4_config["spin"], sym)
    setup_dft_parameters(psi4_config)
    # Copy wfn file to the right place with a right name---
    pid = str(os.getpid())
    targetfile = psi4_scr + filename + '.default.' + pid + '.180.npy'
    if os.path.isfile(psi4_config["wfnfile"]):
        shutil.copyfile(psi4_config["wfnfile"], targetfile)
    # Final scf---
    psi4.set_options({
        "maxiter": 50 if "maxiter" not in psi4_config else psi4_config["maxiter"],
        "D_CONVERGENCE": 3e-5,
        "E_CONVERGENCE": 3e-5,
        "fail_on_maxiter": True})
    if not (("ccsd" in functional) or ("mp2" in functional) or ("scf" in functional)):
        try:
            if (functional not in b3lyp_d) and ("hfx_" not in functional) and ("ccsd" not in functional):
                e, wfn = psi4.energy(functional, molecule=mol, return_wfn=True)
            elif "hfx_" in functional:
                basefunc, hfx = functional.split("_")[0], int(functional.split("_")[-1])
                print("HFX sampling: ", basefunc, hfx)
                e, wfn = psi4.energy("scf", dft_functional=get_hfx_functional(basefunc, hfx),  molecule=mol, return_wfn=True)
            else:
                print("customized b3lyp with different HFX: ", functional)
                e, wfn = psi4.energy("scf", dft_functional=b3lyp_d[functional],  molecule=mol, return_wfn=True)
            if return_wfn:
                wfn.to_file("wfn.180")
        except:
            print("This calculation does not converge.")
    else:
        print("running CC: ", functional)
        psi4.set_options({
            'reference': d['ref'].replace("ks", "hf"),
            'R_CONVERGENCE': 1e-5,
            'E_CONVERGENCE': 5e-5,
            'D_CONVERGENCE': 5e-5,
            "mp2_type": "df",
            "cc_type": "conv",
            "scf_type": "df",
            'nat_orbs': True,
            'FREEZE_CORE': True,
            "GUESS": "SAD",
            })
        e, wfn = psi4.energy(functional, molecule=mol, return_wfn=True)
        wfn.to_file("wfn.180")
    success = check_sucess()
    for filename in os.listdir("./"):
        if ("psi." in filename) or ("default" in filename):
            print("removing: :", filename)
            os.remove(filename)
    os.chdir(basedir)
    return success

#moved to psi4_utils, edited there to run in the parent directory.
def run_general_hfx(psi4_config, functional, hfx, wfn):
    psi4_scr = './'
    filename = "output"
    basedir = os.getcwd()
    rundir = "./" + functional + "-%d" % hfx
    with open(psi4_config["charge-spin-info"], "r") as f:
        d = json.load(f)
    psi4_config.update(d)
    shutil.copyfile(psi4_config["xyzfile"], os.path.join(rundir, psi4_config["xyzfile"]))
    ensure_dir(rundir)
    os.chdir(rundir)
    psi4.core.set_output_file(filename + '.dat', False)
    sym = 'c1' if 'sym' not in psi4_config else psi4_config['sym']
    mol = get_molecule(psi4_config["xyzfile"], psi4_config["charge"], psi4_config["spin"], sym)
    setup_dft_parameters(psi4_config)
    # Copy wfn file to the right place with a right name---
    pid = str(os.getpid())
    targetfile = psi4_scr + filename + '.default.' + pid + '.180.npy'
    if not os.path.isfile(wfn):
        print("Previous calculation failed... This one is skipped.")
        return False
    shutil.copyfile(wfn, targetfile)
    # Final scf---
    psi4.set_options({
        "maxiter": 50 if "maxiter" not in psi4_config else psi4_config["maxiter"],
        "D_CONVERGENCE": 3e-5,
        "E_CONVERGENCE": 3e-5,
        "fail_on_maxiter": True})
    try:
        e, wfn_o = psi4.energy("scf", molecule=mol, return_wfn=True, dft_functional=get_hfx_functional(functional, hfx))
        wfn_o.to_file("wfn.180")
        # os.remove(wfn)
    except:
        print("This calculation does not converge.")
    success = check_sucess()
    for filename in os.listdir("./"):
        if ("psi." in filename) or ("default" in filename):
            print("removing: :", filename)
            os.remove(filename)
    os.chdir(basedir)
    return success

#moved to psi4_utils
def get_hfx_functional(functional, hfx):
    fmap = {"tpss": "TPSS", "scan": "SCAN", "m06-l": "M06_L", "mn15-l": "MN15_L"}
    if functional == "bp86":
        hfx_func = {
            "name": "hfx_func",
            "x_functionals": {"GGA_X_B88": {"alpha": 1-hfx*0.01}},
            "x_hf": {"alpha": hfx*0.01},
            "c_functionals": {"GGA_C_P86": {}}
        }
    elif functional == "blyp":
        hfx_func = {
            "name": "hfx_func",
            "x_functionals": {"GGA_X_B88": {"alpha": 1-hfx*0.01}},
            "x_hf": {"alpha": hfx*0.01},
            "c_functionals": {"GGA_C_LYP": {}}
        }
    elif functional == "b3lyp":
        hfx_func = {
            "name": "hfx_func",
            "x_functionals": {
                "GGA_X_B88": {"alpha": 0.9*(1-hfx*0.01)},
                "LDA_X": {"alpha": 0.1*(1-hfx*0.01)}
                    },
            "x_hf": {"alpha": hfx*0.01},
            "c_functionals": {
                "GGA_C_LYP": {"alpha": 0.81},
                "LDA_C_VWN_RPA": {"alpha": 0.19}
            }
        }
    elif functional == "pbe":
        hfx_func = {
            "name": "hfx_func",
            "x_functionals": {"GGA_X_PBE": {"alpha": 1-hfx*0.01}},
            "x_hf": {"alpha": hfx*0.01},
            "c_functionals": {"GGA_C_PBE": {}}
        }
    elif functional in ["m06-l", "mn15-l", "scan", "tpss"]:
        mega = "" if "PBE" in functional else "M"
        hfx_func = {
            "name": "hfx_func",
            "x_functionals": {"%sGGA_X_%s" % (mega, fmap[functional]): {"alpha": 1-hfx*0.01}},
            "x_hf": {"alpha": hfx*0.01},
            "c_functionals": {"%sGGA_C_%s" % (mega, fmap[functional]): {}}
        }
    else:
        raise ValueError("This functional has not been implemented with HFX resampling yet: ", functional)
    return hfx_func

#moved to run_utils
def write_jobscript(psi4_config):
    if "cluster" not in psi4_config:
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
                fo.write("python -u loop_derivative_jobs.py  > $SGE_O_WORKDIR/deriv_nohup1.out\n")
            else:
                fo.write("python -u loop_run.py  > $SGE_O_WORKDIR/nohup1.out 2> $SGE_O_WORKDIR/nohup1.err\n")
                fo.write("python -u loop_run.py  > $SGE_O_WORKDIR/nohup2.out 2> $SGE_O_WORKDIR/nohup1.err\n")
                fo.write("python -u loop_run.py  > $SGE_O_WORKDIR/nohup3.out 2> $SGE_O_WORKDIR/nohup1.err\n")
                if "hfx_rescue" in psi4_config and psi4_config["hfx_rescue"]:
                    fo.write("echo rescuing...\n")
                    fo.write("python -u loop_rescue.py > $SGE_O_WORKDIR/rescue_nohup1.out\n")
                    fo.write("python -u loop_rescue.py > $SGE_O_WORKDIR/rescue_nohup2.out\n")
                    fo.write("python -u loop_rescue.py > $SGE_O_WORKDIR/rescue_nohup3.out\n")
            fo.write("cp -rf * $SGE_O_WORKDIR\n")
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
                fo.write("python -u loop_derivative_jobs.py  > $subdir/deriv_nohup1.out\n")
                fo.write("rm */*/psi.* */*/dfh.* */*-*/*.npy */b3lyp/*.molden */b3lyp/*1step*\n")
            else:
                fo.write("python -u loop_run.py  > $subdir/nohup1.out\n")
                fo.write("python -u loop_run.py  > $subdir/nohup2.out\n")
                fo.write("python -u loop_run.py  > $subdir/nohup3.out\n")
                if "hfx_rescue" in psi4_config and psi4_config["hfx_rescue"]:
                    fo.write("echo rescuing...\n")
                    fo.write("python -u loop_rescue.py > $subdir/rescue_nohup1.out\n")
                    fo.write("python -u loop_rescue.py > $subdir/rescue_nohup2.out\n")
                    fo.write("python -u loop_rescue.py > $subdir/rescue_nohup3.out\n")
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
                fo.write("python -u loop_derivative_jobs.py  > deriv_nohup1.out\n")
            else:
                fo.write("python -u loop_run.py  > nohup1.out\n")
                fo.write("python -u loop_run.py  > nohup2.out\n")
                fo.write("python -u loop_run.py  > nohup3.out\n")
                if "hfx_rescue" in psi4_config and psi4_config["hfx_rescue"]:
                    fo.write("echo rescuing...\n")
                    fo.write("python -u loop_rescue.py > rescue_nohup1.out\n")
                    fo.write("python -u loop_rescue.py > rescue_nohup2.out\n")
                    fo.write("python -u loop_rescue.py > rescue_nohup3.out\n")
    elif psi4_config["cluster"] == "mustang":
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

#moved to run_utils.py. Edited so that the files are copied to the parent directory, so the files/commands may need to be modified to match.
def run_bash(cmd, basedir, rundir):
    os.chdir(rundir)
    infile = resource_files("jobmanager").joinpath("psi4_utils/loop_run.py")
    shutil.copy(infile, "./")
    infile_rescue = resource_files("jobmanager").joinpath("psi4_utils/loop_rescue.py")
    shutil.copy(infile_rescue, "./")
    infile_deriv = resource_files("jobmanager").joinpath("psi4_utils/loop_derivative_jobs.py")
    shutil.copy(infile_deriv, "./")
    print("Executing: ", cmd, "at: ", rundir)
    subprocess.call(cmd, shell=True)
    os.chdir(basedir)

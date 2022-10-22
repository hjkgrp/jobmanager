import glob
import numpy as np
import os
import shutil
import subprocess

from molSimplifyAD.utils.pymongo_tools import connect2db, query_lowestE_converged

from .tools import manager_io


def isCSD(job):
    try:
        iscsd = True
        for ii in range(6):
            if not job[ii].isupper():
                iscsd = False
                break
        return iscsd
    except:
        return False


def call_molsimplify(geodir, job, jobname):
    liglist = ",".join(job["ligstr"].split("_"))  # can be a single SMILES string, or list of database ligands (e.g. water_water_water_water_water_water)
    tmp_name = str(np.random.randint(10 ** 18))  # assign a temporary name so that the results are findable
    temp_rundir = os.path.join(os.path.expanduser('~'), 'Runs')
    bash_command = " ".join(["molsimplify ", '-core ' + job["metal"],
                             '-lig ' + str(liglist),
                             '-ligloc ' + 'yes', '-calccharge yes',
                             '-spin ' + str(job["spin"]), '-oxstate ' + str(job["ox"]),
                             "-ffoption " + "b", ' -ff UFF',
                             "-name", tmp_name, "-rundir", temp_rundir])
    if "geometry" in job:
        bash_command = " ".join([bash_command, "-geometry", job["geometry"]])
    if "coord" in job:
        bash_command = " ".join([bash_command, "-coord", str(job["coord"])])
    if "ligocc" in job:  # must be a string, e.g. "6" or "1,1,1,1,1,1"
        bash_command = " ".join([bash_command, "-ligocc", job["ligocc"]])
    else:
        bash_command = " ".join([bash_command, "-ligocc", "1,1,1,1,1,1"])
    if "keepHs" in job:
        bash_command = " ".join([bash_command, "-keepHs", job["keepHs"]])
    else:
        bash_command = " ".join([bash_command, "-keepHs", 'yes,yes,yes,yes,yes,yes'])
    if "smicat" in job:  # must be a string, e.g. "1"
        bash_command = " ".join([bash_command, "-smicat", job["smicat"]])
    if "skipANN" in job:
        bash_command = " ".join([bash_command, "-skipANN", job["skipANN"]])
    print(("call: ", bash_command))
    bash_command = bash_command.split()
    subprocess.call(bash_command)

    file_name = os.path.join(temp_rundir, tmp_name)
    print(("file_name: ", file_name))
    print(("geodir: ", geodir))
    inner_folder_path = glob.glob(os.path.join(file_name, '*'))[0]
    xyz_path = glob.glob(os.path.join(inner_folder_path, '*.xyz'))[0]
    charge = False
    with open(inner_folder_path + '/terachem_input', "r") as fo:
        for line in fo:
            if "charge" == line[:len("charge")]:
                charge = int(line.split()[-1])
    if not charge:
        raise ValueError("No charge is extracted from terachem input: ", inner_folder_path + '/terachem_input')
    shutil.copyfile(xyz_path, geodir + '/' + jobname + '.xyz')
    print((xyz_path, geodir + '/' + jobname + '.xyz'))
    return charge


def write_xyz_from_db(geodir, jobname, optgeo):
    natoms = len(optgeo.split('\n')) - 1
    with open(geodir + '/' + jobname + ".xyz", "w") as fo:
        fo.write("%d\n" % natoms)
        fo.write("====Geometry adopted from the database====\n")
        fo.write(optgeo)


def generate_fake_results_from_db(rundir, jobname, tmcdoc):
    scrdir = rundir + '/scr'
    if not os.path.isdir(scrdir):
        os.makedirs(scrdir)
    _ = write_xyz_from_db(scrdir, 'optim', tmcdoc["opt_geo"])
    if tmcdoc['wavefunction']:
        if int(tmcdoc['spin']) == 1:
            try:
                shutil.copy(tmcdoc['wavefunction']['c0'], scrdir + '/c0')
            except FileNotFoundError:
                pass
        else:
            try:
                shutil.copy(tmcdoc['wavefunction']['ca0'], scrdir + '/ca0')
                shutil.copy(tmcdoc['wavefunction']['cb0'], scrdir + '/cb0')
            except FileNotFoundError:
                pass
    inpath = rundir + '/' + jobname + '.in'
    with open(inpath, "w") as fo:
        fo.write("=======This outfile is FAKE and generated artificialy======\n")
        fo.write("charge %d\n" % int(tmcdoc["charge"]))
        fo.write("spinmult %d\n" % int(tmcdoc["spin"]))
        fo.write("method %s\n" % str(tmcdoc["functional"]))
        fo.write("basis %s\n" % str(tmcdoc["basis"]))
        fo.write("coordinates %s.xyz\n" % jobname)
        fo.write("run minimize\n")
    outpath = rundir + '/' + jobname + '.out'
    with open(outpath, "w") as fo:
        fo.write("=======This outfile is FAKE and generated artificialy======\n")
        fo.write('Startfile from command line: %s.in\n' % jobname)
        fo.write('Hartree-Fock exact exchange:          %.2f\n' % (float(tmcdoc['alpha']) * 0.01))
        fo.write('DFT Functional requested: %s\n' % str(tmcdoc['functional']))
        fo.write('*                    TeraChem %s            *\n' % str(tmcdoc['terachem_version']))
        fo.write('Alpha level shift: 0.25\n')
        fo.write('Beta level shift: 0.25\n')
        fo.write('Using basis set: %s\n' % str(tmcdoc["basis"]))
        fo.write("Total charge:    %d\n" % int(tmcdoc["charge"]))
        fo.write("Spin multiplicity: %d\n" % int(tmcdoc["spin"]))
        fo.write("SPIN S-SQUARED: %f (exact: %f)\n" % (float(tmcdoc["ss_act"]), float(tmcdoc["ss_target"])))
        fo.write("-=#=-      Optimization Cycle     1   -=#=-\n")
        fo.write("FINAL ENERGY: %.10f a.u.\n" % float(tmcdoc["energy"]))
        fo.write("-=#=-     Optimization Converged.     -=#=-\n")
        fo.write("Total processing time: 0.00 sec\n")
        fo.write("Job finished: A time of mystery\n")
    return outpath


def populate_single_job(basedir, job, db, safe_filenames=True):
    geodir = basedir + "/initial_geometry/"
    if not os.path.isdir(geodir):
        os.makedirs(geodir)
    iscsd = isCSD(job['ligstr'])
    query_constraints = {"metal": job['metal'], "spin": job["spin"], "ligstr": job["ligstr"], "alpha": 20,
                         "wavefunction": {"$exists": True}}
    if not iscsd:
        query_constraints.update({"ox": job["ox"]})
        jobname = "_".join([job['metal'], str(job['ox']), str(job['spin']), job['ligstr']])
    else:
        jobname = "_".join([job['ligstr'], job['metal'], str(job['spin'])])
    tmcdoc, recover = None, True
    if db is not None:
        tmcdoc = query_lowestE_converged(db, collection='oct', constraints=query_constraints)
        if tmcdoc is not None:
            print(("Bingo! Optimized geometry found in db: ", query_constraints))
            try:
                charge = int(tmcdoc["charge"])
                ss_act, ss_target = float(tmcdoc["ss_act"]), float(tmcdoc["ss_target"])
                if abs(ss_act - ss_target) > 1 and tmcdoc["ss_flag"] == 1:
                    recover = False
                    print("Spin contamination for singlets! (used ub3lyp)")
                write_xyz_from_db(geodir, jobname, tmcdoc["opt_geo"])
                wfnpath = "/home/data/wfn/" + str(tmcdoc['unique_name'])
                wfnfiles = os.listdir(wfnpath)
                if not any(x in ["c0", "ca0", "cb0"] and os.stat(wfnpath+'/%s' % x).st_size > 1000 for x in wfnfiles):
                    recover = False
                    print("No WFN file found at /home/data/wfn/.")
            except:
                recover = False
        else:
            print("generate initial geometry from molsimplify...")
            charge = call_molsimplify(geodir, job, jobname)
    else:
        if not iscsd:
            print("NO db connection! Generate initial geometry from molsimplify...")
            jobname_safe = jobname.replace("#", "3").replace("(", "[").replace(")", "]")
            charge = call_molsimplify(geodir, job, jobname_safe)  # this used to say rundir for some reason...
        else:
            raise ValueError("Cannot generate initial geometry for CSD complex...")

    rundir = basedir + '/' + jobname
    try:
        rundir_p3 = basedir + '/' + jobname_safe
    except NameError:  # jobname_safe not defined
        rundir_p3 = basedir + '/' + jobname.replace('#', '3')

    # p3 option
    if safe_filenames:
        rundir = rundir_p3
        try:
            jobname = jobname_safe
        except NameError:  # jobname_safe not defined
            jobname = jobname.replace('#', '3')

    populated = True
    if not os.path.isdir(rundir) and not os.path.isdir(rundir_p3) and recover:
        os.makedirs(rundir)
        shutil.copyfile(geodir + '/' + jobname + '.xyz', rundir + "/" + jobname + ".xyz")
        os.chdir(rundir)
        # Add fake files etc for a smooth carry-on in job manager for further dependent jobs.
        if tmcdoc is not None:
            generate_fake_results_from_db(rundir, jobname, tmcdoc)
        else:
            manager_io.write_input(name=jobname, coordinates=jobname + '.xyz', charge=charge, spinmult=int(job["spin"]), run_type='minimize', solvent=False)
            manager_io.write_jobscript(jobname)
    elif os.path.isdir(rundir) or os.path.isdir(rundir_p3):
        print("folder exist.")
        populated = False
    else:
        print(('WARNING: cannot recover %s' % jobname))
        populated = False
    os.chdir(basedir)
    return jobname, populated


def populate_list_of_jobs(basedir, jobs, db_communicate=True):
    if db_communicate:
        db = connect2db(user="readonly_user", pwd="readonly", host="localhost",
                        port=27017, database="tmc", auth=True)
    else:
        db = None
    for job in jobs:
        populate_single_job(basedir, job, db)

#!/usr/bin/env python
import os
import numpy as np
import time
import sys
import json
import jobmanager.tools as tools
import jobmanager.moltools as moltools
import jobmanager.recovery as recovery
from jobmanager.io import io
from jobmanager.psi4_utils.run import write_jobscript, run_bash
import jobmanager.classes as classes

def kill_jobs(kill_names, message1='Killing job: ', message2=' early'):
    """This function takes a list of job names and kills the jobs associated with them, if the jobs are active

        Parameters
        ----------
            kill_names : list
                List of jobs to kill.
            message1 : str
                Message prefix to report to stdout.
            message2 : str
                Message suffix to report to stdout.

    """
    # This function takes a list of job names and kills the jobs associated with them, if the jobs are active
    if type(kill_names) != list:
        kill_names = [kill_names]
    machine = tools.get_machine()

    active_jobs, active_ids = tools.list_active_jobs(ids=True)
    active_jobs = list(zip(active_jobs, active_ids))

    jobs_to_kill = [[name, id_] for name, id_ in active_jobs if name in kill_names]

    for name, id_ in jobs_to_kill:
        print(message1 + name + message2)
        if machine in ['gibraltar']:
            tools.call_bash('qdel ' + str(id_))
        elif machine in ['comet', 'bridges','expanse']:
            tools.call_bash('scancel '+str(id_))
        else:
            raise ValueError('Sardines.')


def prep_derivative_jobs(directory, list_of_outfiles):
    """This function takes a directory and output files and spawns derivative jobs.

        Parameters
        ----------
            directory : str
                Directory of interest to analyze.
            list_of_outfiles : list
                List of output files that are read to spawn derivative jobs.

    """
    for job in list_of_outfiles:
        configure_dict = io.read_configure(directory, job)
        if configure_dict['solvent']:
            tools.prep_solvent_sp(job, configure_dict['solvent'])
        if configure_dict['functionalsSP']:
            tools.prep_functionals_sp(job, configure_dict['functionalsSP'])
        if configure_dict['vertEA']:
            tools.prep_vertical_ea(job)
        if configure_dict['vertIP']:
            tools.prep_vertical_ip(job)
        if configure_dict['thermo']:
            tools.prep_thermo(job)
        if configure_dict['hfx_resample']:
            tools.prep_hfx_resample(job)
        if configure_dict['dissociation']:
            moltools.prep_ligand_breakdown(job, dissociated_ligand_charges=configure_dict['dissociated_ligand_charges'],
                                           dissociated_ligand_spinmults=configure_dict['dissociated_ligand_spinmults'])
        if configure_dict['mbe']:
            moltools.prep_mbe_calc(job)  # needs to be generalized, not just for Fe
            # moltools.prep_mbe_calc(job, metal_charge = configure_dict['metal_charge'])
        if configure_dict['spinSplitting']:
            tools.prep_ad_spin(job)
        if bool(configure_dict['general_sp']):
            tools.prep_general_sp(job, general_config=configure_dict['general_sp'])


def resub(directory='in place'):
    """This function takes a directory and submits calculations.

        Parameters
        ----------
            directory : str
                Directory of interest to analyze.

    """
    # Takes a directory, resubmits errors, scf failures, and spin contaminated cases
    configure_dict = io.read_configure(directory, None)
    max_resub = configure_dict['max_resub']
    max_jobs = configure_dict['max_jobs']
    hard_job_limit = configure_dict['hard_job_limit']
    hit_queue_limit = False  # Describes if this run has limitted the number of jobs submitted to work well with the queue
    # Get the state of all jobs being managed by this instance of the job manager
    completeness = moltools.check_completeness(directory, max_resub, configure_dict=configure_dict)
    # print("completeness: ", completeness)
    errors = completeness['Error']  # These are calculations which failed to complete
    scf_errors = completeness['SCF_Error']  # These are calculations which failed to complete, appear to have an scf error, and hit wall time
    oscillating_scf_errors = completeness['oscillating_scf_errors']  # These are calculations which failed to complete, appear to have an oscillaing scf error,
    need_resub = completeness['Needs_resub']  # These are calculations with level shifts changed or hfx exchange changed
    spin_contaminated = completeness['Spin_contaminated']  # These are finished jobs with spin contaminated solutions
    thermo_grad_error = completeness['Thermo_grad_error']  # These are thermo jobs encountering the thermo grad error
    waiting = completeness['Waiting']  # These are jobs which are or were waiting for another job to finish before continuing.
    bad_geos = completeness['Bad_geos']  # These are jobs which finished, but converged to a bad geometry.
    finished = completeness['Finished']
    molscontrol_kills = completeness['molscontrol_kills']
    nactive = tools.get_number_active()  # number of active jobs, counting bundled jobs as a single job
    # Kill SCF errors in progress, which are wasting computational resources
    all_scf_errors = completeness['SCF_Errors_Including_Active']  # These are all jobs which appear to have scf error, including active ones
    scf_errors_to_kill = [scf_err for scf_err in all_scf_errors if scf_err not in scf_errors]
    names_to_kill = [os.path.split(scf_err)[-1].rsplit('.', 1)[0] for scf_err in scf_errors_to_kill]
    names_to_kill.extend([name + '_jobscript' for name in names_to_kill])
    kill_jobs(names_to_kill, message1='Job: ', message2=' appears to have an scf error. Killing this job early')
    # Prep derivative jobs such as thermo single points, vertical IP, and ligand dissociation energies
    needs_derivative_jobs = list(filter(tools.check_original, finished))
    # print("needs_derivative_jobs: ", needs_derivative_jobs)
    prep_derivative_jobs(directory, needs_derivative_jobs)
    resubmitted = []  # Resubmitted list gets True if the job is submitted or False if not. Contains booleans, not job identifiers.
    counter = 0
    for job in molscontrol_kills:
        counter += 1
        print("killed by molscontrol: ", job, counter)

    # List jobs which might need human intervention
    for error in completeness['Terminated_No_SCF_cycle']:
        print(f"Job terminated without any SCF cycles, check for validity: {error}")

    # Resub unidentified errors
    for error in errors:
        if ((nactive + np.sum(resubmitted)) >= max_jobs) or (
                (tools.get_total_queue_usage() + np.sum(resubmitted)) >= hard_job_limit):
            hit_queue_limit = True
            continue
        resub_tmp = recovery.simple_resub(error)
        if resub_tmp:
            print(('Unidentified error in job: ' + os.path.split(error)[-1] + ' -Resubmitting'))
            print('')
        resubmitted.append(resub_tmp)

    # Resub oscillating_scf convergence errors
    for error in oscillating_scf_errors:
        if ((nactive + np.sum(resubmitted)) >= max_jobs) or (
                (tools.get_total_queue_usage() + np.sum(resubmitted)) >= hard_job_limit):
            hit_queue_limit = True
            continue
        local_configure = io.read_configure(directory, None)
        if 'scf' in local_configure['job_recovery']:
            resub_tmp = recovery.resub_oscillating_scf(error)
            if resub_tmp:
                print(('Oscillating SCF error identified in job: ' + os.path.split(error)[
                    -1] + ' -Resubmitting with adjusted precision and grid.'))
                print('')
            resubmitted.append(resub_tmp)

    # Resub scf convergence errors
    for error in scf_errors:
        if ((nactive + np.sum(resubmitted)) >= max_jobs) or (
                (tools.get_total_queue_usage() + np.sum(resubmitted)) >= hard_job_limit):
            hit_queue_limit = True
            continue
        local_configure = io.read_configure(directory, None)
        if 'scf' in local_configure['job_recovery']:
            resub_tmp = recovery.resub_scf(error)
            if resub_tmp:
                print(('SCF error identified in job: ' + os.path.split(error)[
                    -1] + ' -Resubmitting with adjusted levelshifts'))
                print('')
            resubmitted.append(resub_tmp)

    # Resub jobs which converged to bad geometries with additional constraints
    for error in bad_geos:
        if ((nactive + np.sum(resubmitted)) >= max_jobs) or (
                (tools.get_total_queue_usage() + np.sum(resubmitted)) >= hard_job_limit):
            hit_queue_limit = True
            continue
        local_configure = io.read_configure(directory, None)
        if 'bad_geo' in local_configure['job_recovery']:
            resub_tmp = recovery.resub_bad_geo(error, directory)
            if resub_tmp:
                print(('Bad final geometry in job: ' + os.path.split(error)[
                    -1] + ' -Resubmitting from initial structure with additional constraints'))
                print('')
            resubmitted.append(resub_tmp)

    # Resub spin contaminated cases
    for error in spin_contaminated:
        if ((nactive + np.sum(resubmitted)) >= max_jobs) or (
                (tools.get_total_queue_usage() + np.sum(resubmitted)) >= hard_job_limit):
            hit_queue_limit = True
            continue
        local_configure = io.read_configure(directory, None)
        if 'spin_contaminated' in local_configure['job_recovery']:
            resub_tmp = recovery.resub_spin(error)
            if resub_tmp:
                print(('Spin contamination identified in job: ' + os.path.split(error)[
                    -1] + ' -Resubmitting with adjusted HFX'))
                print('')
            resubmitted.append(resub_tmp)

    # Resub jobs with atypical parameters used to aid convergence
    for error in need_resub:
        if ((nactive + np.sum(resubmitted)) >= max_jobs) or (
                (tools.get_total_queue_usage() + np.sum(resubmitted)) >= hard_job_limit):
            hit_queue_limit = True
            continue
        resub_tmp = recovery.clean_resub(error)
        if resub_tmp:
            print(('Job ' + os.path.split(error)[-1] + ' needs to be rerun with typical paramters. -Resubmitting'))
            print('')
        resubmitted.append(resub_tmp)

    # Create a job with a tighter convergence threshold for failed thermo jobs
    for error in thermo_grad_error:
        if ((nactive + np.sum(resubmitted)) >= max_jobs) or (
                (tools.get_total_queue_usage() + np.sum(resubmitted)) >= hard_job_limit):
            hit_queue_limit = True
            continue
        local_configure = io.read_configure(directory, None)
        if 'thermo_grad_error' in local_configure['job_recovery']:
            resub_tmp = recovery.resub_tighter(error)
            if resub_tmp:
                print(('Job ' + os.path.split(error)[
                    -1] + ' needs a better initial geo. Creating a geometry run with tighter convergence criteria'))
                print('')
            resubmitted.append(resub_tmp)

    # Look at jobs in "waiting," resume them if the job they were waiting for is finished
    # Currently, this should only ever be thermo jobs waiting for an ultratight job
    for waiting_dict in waiting:
        if ((nactive + np.sum(resubmitted)) >= max_jobs) or (
                (tools.get_total_queue_usage() + np.sum(resubmitted)) >= hard_job_limit):
            hit_queue_limit = True
            continue
        if len(list(waiting_dict.keys())) > 1:
            raise Exception('Waiting job list improperly constructed')
        job = list(waiting_dict.keys())[0]
        waiting_for = waiting_dict[job]
        if waiting_for in finished:
            history = recovery.load_history(job)
            history.waiting = None
            history.save()
            results_for_this_job = io.read_outfile(job)
            if results_for_this_job['thermo_grad_error']:
                resubmitted.append(recovery.resub_thermo(job))
            else:
                raise Exception('A method for resuming job: ' + job + ' is not defined')
        else:
            resubmitted.append(False)

    # Submit jobs which haven't yet been submitted
    if (((nactive + np.sum(resubmitted)) < max_jobs) and (
            (tools.get_total_queue_usage() + np.sum(resubmitted)) < hard_job_limit)):
        to_submit = []
        jobscripts = tools.find('*_jobscript')
        active_jobs = tools.list_active_jobs(home_directory=directory, parse_bundles=True)
        for job in jobscripts:
            jobdir, jobname = os.path.split(job.rsplit('_', 1)[0])
            outfile = os.path.join(jobdir, jobname + '.out')
            if not os.path.isfile(outfile) and not any([x.startswith(jobname) for x in active_jobs]):
                to_submit.append(job)

        short_jobs_to_submit = [i for i in to_submit if tools.check_short_single_point(i)]
        long_jobs_to_submit = [i for i in to_submit if i not in short_jobs_to_submit]
        if len(short_jobs_to_submit) > 0:
            bundled_jobscripts = tools.bundle_jobscripts(os.getcwd(), short_jobs_to_submit)
        else:
            bundled_jobscripts = []
        to_submit = bundled_jobscripts + long_jobs_to_submit

        submitted = []
        invalid_jobs = []

        for job in to_submit:
            # update the number of running + submitted jobs (from this instance of jobmanager)
            dynamic_nactive = (len(submitted) + nactive + np.sum(resubmitted))
            # update the number of jobs the user is currently running
            user_nactive = tools.get_total_queue_usage() + len(submitted) + np.sum(resubmitted)
            # make sure we don't exceed set job limits
            if (dynamic_nactive >= max_jobs) or (user_nactive >= hard_job_limit):
                hit_queue_limit = True
                continue

            #get the path to input file
            inp_file_path = job.rsplit('_',1)[0]+'.in'
            inp_file = classes.textfile(inp_file_path)

            #read lines of input file
            inp_lines = inp_file.lines


            #check if Terachem or ORCA
            terachem = True
            for line in inp_lines:
                if not line.strip():
                    continue
                if line.strip()[0] == '!':
                    terachem = False
                    break

            #if terachem, check validity: if valid, submit, else add to invalid jobs list
            if terachem:
                tc_dict = io.read_terachem_input(inp_file)
                if io.spinchargeChecker(tc_dict,inp_file_path):
                    print(('Initial submission for job: ' + os.path.split(job)[-1]))
                    tools.qsub(job)
                    submitted.append(True)
                else:
                    invalid_jobs.append(job)
                    print(f'Invalid job: {os.path.split(job)[-1]}')
            else:  # have not implemented spin/charge checker for ORCA
                print(('Initial submission for job: ' + os.path.split(job)[-1]))
                tools.qsub(job)
                submitted.append(True)

    else:
        print('==== Hit the queue limit for the user, not submitting any more jobs. ====')
        hit_queue_limit = True
        submitted = []

    number_resubmitted = np.sum(np.array(resubmitted + submitted))
    # ~ print str(number_resubmitted)+' Jobs submitted'
    return int(number_resubmitted), int(len(completeness['Active'])), hit_queue_limit


def resub_psi4(psi4_config):
    basedir = os.getcwd()
    if "trigger" in psi4_config:
        write_jobscript(psi4_config)
        if "cluster" not in psi4_config or psi4_config["cluster"] == "mustang":
            cmd = "qsub jobscript.sh"
        else:
            cmd = "sbatch jobscript.sh"
        run_bash(cmd=cmd,
                 basedir=basedir,
                 rundir=basedir)
        time.sleep(3)
    else:
        for path in os.listdir(basedir):
            if os.path.isdir(basedir + "/" + path):
                print("at: ", basedir + "/" + path)
                os.chdir(basedir + "/" + path)
                with open("psi4_config.json", "w") as fo:
                    json.dump(psi4_config, fo)
                write_jobscript(psi4_config)
                os.chdir(basedir)
                if "cluster" not in psi4_config or psi4_config["cluster"] == "mustang":
                    cmd = "qsub jobscript.sh"
                else:
                    cmd = "sbatch jobscript.sh"
                run_bash(cmd=cmd,
                         basedir=basedir,
                         rundir=basedir + "/" + path)
                time.sleep(3)


def main():
    counter = 0
    configure_dict = io.read_configure('in place', None)
    print("configure_dict: ", configure_dict)
    if not configure_dict["run_psi4"]:
        while True:
            print('**********************************')
            print("****** Assessing Job Status ******")
            print('**********************************')
            time1 = time.time()
            with open('complete', 'w') as fil:
                fil.write('Active')

            dir = os.getcwd()
            tools.check_queue(dir)

            number_resubmitted, number_active, hit_queue_limit = resub()

            print('**********************************')
            print(("******** " + str(number_resubmitted) + " Jobs Submitted ********"))
            print('**********************************')

            print(('job cycle took: ' + str(time.time() - time1)))
            print(('sleeping for: ' + str(configure_dict['sleep'])))
            sys.stdout.flush()
            time.sleep(configure_dict[
                           'sleep'])  # sleep for time specified in configure. If not specified, default to 7200 seconds (2 hours)

            # Terminate the script if it is no longer submitting jobs
            if number_resubmitted == 0 and number_active == 0 and not hit_queue_limit:
                counter += 1
            else:
                counter = 0
            if counter >= 3:
                break

        print('**********************************')
        print("****** Normal Terminatation ******")
        print('**********************************')
        with open('complete', 'w') as fil:
            fil.write('True')
    else:
        resub_psi4(configure_dict["psi4_config"])


if __name__ == '__main__':
    main()

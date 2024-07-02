# This module is equivalent to the tools.py module, but includes dependency on molSimplify
# updates 8/2/19

import os
import copy
import numpy as np
import pandas as pd
import jobmanager.tools as tools
from jobmanager.io import io
from molSimplify.Classes.mol3D import mol3D
from molSimplify.Classes.ligand import ligand_breakdown


def read_run(outfile_PATH):
    """Read a DFT output.

    Parameters
    ----------
        outfile_PATH : str
            Path to quantum chemistry output file.

    Returns
    -------
        results : dict
            Analyzed results for job.

    """
    # Evaluates all aspects of a run using the outfile and derivative files
    results = io.read_outfile(outfile_PATH, long_output=True)
    infile_dict = io.read_infile(outfile_PATH)
    results['levela'], results['levelb'] = infile_dict['levelshifta'], infile_dict['levelshiftb']
    results['method'], results['hfx'] = infile_dict['method'], infile_dict['hfx']
    results['constraints'] = infile_dict['constraints']

    mullpop_path = os.path.join(os.path.split(outfile_PATH)[0], 'scr', 'mullpop')
    if os.path.exists(mullpop_path):
        mullpops = io.read_mullpop(outfile_PATH)
        metal_types = ['Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Mo', 'Tc', 'Ru', 'Rh']
        metals = [i for i in mullpops if i.split()[0] in metal_types]
        if len(metals) > 1:
            results['metal_spin'] = np.nan
        elif len(metals) == 0:
            pass
        else:
            results['metal_spin'] = float(metals[0].split()[-1])
    else:
        results['metal_spin'] = np.nan

    optim_path = os.path.join(os.path.split(outfile_PATH)[0], 'scr', 'optim.xyz')
    initial_xyz_path = outfile_PATH.rsplit('.', 1)[0]+'.xyz'
    if not os.path.isfile(initial_xyz_path):
        raise Exception('No initial xyz found at: '+initial_xyz_path)

    check_geo = False
    if os.path.isfile(optim_path):
        with open(optim_path, 'r') as fil:
            lines = fil.readlines()
        if len(lines) > 0:
            check_geo = True  # Only apply geo check if an optimized geometry exists

    if check_geo:
        tools.extract_optimized_geo(optim_path)
        optimized_path = os.path.join(os.path.split(optim_path)[0], 'optimized.xyz')

        mol = mol3D()
        mol.readfromxyz(optimized_path)

        IsOct, flag_list, oct_check = mol.IsOct(dict_check=mol.dict_oct_check_st,
                                                silent=True)

        if IsOct:
            IsOct = True
        else:
            IsOct = False

        results['Is_Oct'] = IsOct
        results['Flag_list'] = flag_list
        results['Oct_check_details'] = oct_check

    else:
        results['Is_Oct'] = None
        results['Flag_list'] = None
        results['Oct_check_details'] = None

    return results


def create_summary(directory=None):
    """Create a summary file.

    Parameters
    ----------
        directory : str, optional
            Directory that contains all of the jobs to be analyzed. Default is in place.

    Returns
    -------
        summary : pd.DataFrame
            Summary of full directory.

    """
    if directory is None:
        directory = os.getcwd()
    # Returns a pandas dataframe which summarizes all outfiles in the directory, defaults to cwd

    outfiles = tools.find_calcs(directory, extension='.out')
    outfiles = list(filter(tools.check_valid_outfile, outfiles))
    results = list(map(read_run, outfiles))
    summary = pd.DataFrame(results)

    return summary


def apply_geo_check(job_outfile_path, geometry):
    """Create a summary file.

    Parameters
    ----------
        job_outfile_path : str
            Path for output file.
        geometry : str
            Type of geometry to do the geometry check.

    Returns
    -------
        geo_flag : bool
            Flag describing geometry. True if good geometry.

    """
    if geometry:  # The geometry variable is set to False if no geo check is requested for this job

        optim_path = os.path.join(os.path.split(job_outfile_path)[0], 'scr', 'optim.xyz')

        if os.path.isfile(optim_path):
            lines = tools.extract_optimized_geo(optim_path)
            optimized_path = os.path.join(os.path.split(optim_path)[0], 'optimized.xyz')

            if len(lines)==0:
                return True
            else:
                mol = mol3D()
                mol.readfromxyz(optimized_path)
        else:
            # If the optim.xyz doesn't exist, assume that it's a single point and should pass geo check
            return True

        if geometry.capitalize() in ['Oct', 'Octahedral']:
            geo_check_dict = mol.dict_oct_check_st
            IsOct, flag_list, oct_check = mol.IsOct(dict_check=geo_check_dict, silent=True)
            if IsOct:
                return True
            else:
                return False

        elif geometry.capitalize() == 'Bidentate_oct':
            # Loosened geo check dict appropriate for bidentates
            geo_check_dict = {'num_coord_metal': 6,
                              'rmsd_max': 3, 'atom_dist_max': 0.45,
                              'oct_angle_devi_max': 15, 'max_del_sig_angle': 30,
                              'dist_del_eq': 0.35, 'dist_del_all': 1,
                              'devi_linear_avrg': 20, 'devi_linear_max': 28}
            outer_dict_flags = list(mol.dict_oct_check_st.keys())
            final_dict = dict()
            for key in outer_dict_flags:
                final_dict[key] = geo_check_dict

            IsOct, flag_list, oct_check = mol.IsOct(dict_check=final_dict, silent=True)
            if IsOct:
                return True
            else:
                return False

        elif geometry.capitalize() == 'Tridentate_oct':
            # Extra loosened geo check dict appropriate for Tridentates
            geo_check_dict = {'num_coord_metal': 6,
                              'rmsd_max': 3, 'atom_dist_max': 0.45,
                              'oct_angle_devi_max': 25, 'max_del_sig_angle': 50,
                              'dist_del_eq': 0.35, 'dist_del_all': 1,
                              'devi_linear_avrg': 20, 'devi_linear_max': 28}
            outer_dict_flags = list(mol.dict_oct_check_st.keys())
            final_dict = dict()
            for key in outer_dict_flags:
                final_dict[key] = geo_check_dict

            IsOct, flag_list, oct_check = mol.IsOct(dict_check=final_dict, silent=True)
            if IsOct:
                return True
            else:
                return False

        elif geometry.capitalize() in ['Tetrahedral', 'Tetra']:
            if len(mol.getBondedAtoms(mol.findMetal()[0])) == 4:
                return True
            else:
                return False

        else:
            raise Exception('A check has not been implemented for geometry: ' + geometry)
    else:
        # print('No geomery check requested for job: '+job_outfile_path)
        # print('Passing job without a geometry check')
        return True


def get_metal_and_bonded_atoms(job_outfile, geometry=None):
    """Get metal and bonded atoms of complex.

    Parameters
    ----------
        job_outfile : str
            Path for output file.
        geometry : str, optional
            Type of geometry to define metal bonding. Default is None.

    Returns
    -------
        geo_flag : bool
            Flag describing geometry. True if good geometry.

    """
    # given the path to the outfile of a job, returns a the metal atom index and a list of indices for the metal bonded atoms
    # indices are zero-indexed...Terachem uses 1 indexed lists

    xyz_path = job_outfile.rsplit('.', 1)[0] + '.xyz'
    mol = mol3D()
    mol.readfromxyz(xyz_path)
    metal_index = mol.findMetal()[0]

    if geometry in ['Oct', 'oct', 'Octahedral', 'octahedral']:
        bonded_atom_indices = mol.getBondedAtomsOct(metal_index)
    else:
        print('Warning, generic getBondedAtoms() used for: ' + job_outfile + '. Check behavior')
        bonded_atom_indices = mol.getBondedAtoms(metal_index)

    return metal_index, bonded_atom_indices


def check_completeness(directory=None, max_resub=5, configure_dict=False, verbose=False, finished_prev=None):
    """Get metal and bonded atoms of complex.

    Parameters
    ----------
        directory : str, optional
            Directory where the jobs are running. Default is in place.
        max_resub : int, optional
            Number of resubmissions allowed. Default is 5.
        configure_dict : dict, optional
            Configure file. Default is False.

    Returns
    -------
        completeness : dict
            Completeness dictionary for a given directory.

    """
    if directory is None:
        directory = os.getcwd()
    completeness = tools.check_completeness(directory, max_resub, configure_dict=configure_dict, verbose=verbose, finished_prev=finished_prev)
    # print("=======")
    # print("completeness: ", completeness)
    # The check_completeness() function in tools doesn't check the geometries (because it's molSimplify dependent)
    # Apply the check here to finished and spin contaminated geometries, then update the completeness dictionary
    finished = completeness['Finished']
    spin_contaminated = completeness['Spin_contaminated']
    needs_resub = completeness['Needs_resub']
    unfinished = completeness['Error']
    if verbose:
        print(f"{len(finished)} finished.", flush=True)
        print(f"{len(needs_resub)} need resubmission.")
        print(f"{len(unfinished)} errors found.")
    bad_geos = []
    new_finished = []
    new_spin_contaminated = []
    new_needs_resub = []
    new_unfinished = []
    new_molscontrol_kills = []
    for job in finished:
        goal_geo = io.read_configure(directory, job)['geo_check']
        if apply_geo_check(job, goal_geo):
            new_finished.append(job)
        else:
            bad_geos.append(job)
    for job in spin_contaminated:
        if not check_molscontrol_log(job):
            goal_geo = io.read_configure(directory, job)['geo_check']
            if apply_geo_check(job, goal_geo):
                new_spin_contaminated.append(job)
            else:
                bad_geos.append(job)
        else:
            new_molscontrol_kills.append(job)
    for job in needs_resub:
        if not check_molscontrol_log(job):
            goal_geo = io.read_configure(directory, job)['geo_check']
            if apply_geo_check(job, goal_geo):
                new_needs_resub.append(job)
            else:
                bad_geos.append(job)
        else:
            new_molscontrol_kills.append(job)
    for job in unfinished:
        if not check_molscontrol_log(job):
            goal_geo = io.read_configure(directory, job)['geo_check']
            if apply_geo_check(job, goal_geo):
                new_unfinished.append(job)
            else:
                bad_geos.append(job)
        else:
            new_molscontrol_kills.append(job)

    completeness['Finished'] = new_finished
    completeness['Spin_contaminated'] = new_spin_contaminated
    completeness['Resub'] = new_needs_resub
    completeness['Error'] = new_unfinished
    completeness['Bad_geos'] = bad_geos
    completeness["molscontrol_kills"] = new_molscontrol_kills
    return completeness


def prep_ligand_breakdown(outfile_path, dissociated_ligand_charges={}, dissociated_ligand_spinmults={}):
    """Prep ligand breakdown.

    Parameters
    ----------
        outfile_path : str
            Path to output file.
        dissociated_ligand_charges : dict, optional
            Charges for dissociated ligands. Default is empty.
        dissociated_ligand_spinmults : dict, optional
            Spin multiplicity for dissociated ligands. Default is empty.

    Returns
    -------
        jobscripts : list
            List of jobscripts for ligand breakdown jobs.

    """
    # Given a path to the outfile of a finished run, this preps the files for rigid ligand dissociation energies of all ligands
    # Returns a list of the PATH(s) to the jobscript(s) to start the rigid ligand calculations

    home = os.getcwd()
    machine = tools.get_machine()
    outfile_path = tools.convert_to_absolute_path(outfile_path)

    results = io.read_outfile(outfile_path)
    if not results['finished']:
        raise Exception('This calculation does not appear to be complete! Aborting...')

    infile_dict = io.read_infile(outfile_path)
    charge = int(infile_dict['charge'])
    spinmult = int(infile_dict['spinmult'])

    base = os.path.split(outfile_path)[0]
    name = os.path.split(outfile_path)[-1][:-4]

    breakdown_folder = os.path.join(base, name + '_dissociation')

    if os.path.isdir(breakdown_folder):
        return ['Ligand dissociation directory already exists']

    optimxyz = os.path.join(base, 'scr', 'optim.xyz')
    tools.extract_optimized_geo(optimxyz)

    mol = mol3D()
    mol.readfromxyz(os.path.join(base, 'scr', 'optimized.xyz'))

    ligand_idxs, _, _ = ligand_breakdown(mol, silent=True)

    ligand_syms = []
    for ii in ligand_idxs:
        ligand_syms.append([mol.getAtom(i).symbol() for i in ii])

    ligand_names = name_ligands(ligand_syms)

    if not os.path.isdir(breakdown_folder):
        os.mkdir(breakdown_folder)
    os.chdir(breakdown_folder)

    jobscripts = []
    for ligand in zip(ligand_names, ligand_idxs):

        # Assign charges to use during the breakdown for special cases specified in the configure file
        # All other ligands are assigned charge 0
        if ligand[0] in list(dissociated_ligand_charges.keys()):
            ligand_charge = dissociated_ligand_charges[ligand[0]]
        else:
            ligand_charge = 0
        metal_charge = charge - ligand_charge

        # Assign spin, which always remains with the metal except when the dissociated ligand is defined to have spin (O2 for example)
        if spinmult == 1:  # If the whole complex is restricted, it's components must be restricted as well
            ligand_spin, metal_spin = 1, 1
        else:
            if ligand[0] in list(dissociated_ligand_spinmults.keys()):
                ligand_spin = dissociated_ligand_spinmults[ligand[0]]
            else:
                ligand_spin = 1

            metal_spin = spinmult - ligand_spin + 1  # Derived from spinmult = (2S+1) where S=1/2 per electron

        # Create the necessary files for the metal complex single point
        local_name = name + '_rm_' + ligand[0]
        if os.path.isdir('rm_' + ligand[0]):
            pass
        else:
            os.mkdir('rm_' + ligand[0])
            os.chdir('rm_' + ligand[0])

            local_mol = mol3D()
            local_mol.copymol3D(mol)

            local_mol.deleteatoms(ligand[1])
            local_mol.writexyz(local_name + '.xyz')

            local_infile_dict = copy.copy(infile_dict)
            local_infile_dict['name'] = local_name
            local_infile_dict['coordinates'] = local_name+'.xyz'
            local_infile_dict['charge'], local_infile_dict['spinmult'] = metal_charge, metal_spin
            local_infile_dict['run_type'] = 'energy'
            local_infile_dict['constraints'], local_infile_dict['convergence_thresholds'] = False, False
            local_infile_dict['machine'] = machine

            io.write_input(local_infile_dict)
            io.write_jobscript(local_name, time_limit='12:00:00', machine=machine)
            jobscripts.append(local_name + '.in')
            os.chdir('..')

        # Create the necessary files for the dissociated ligand single point
        local_name = name + '_kp_' + ligand[0]
        if os.path.isdir('kp_' + ligand[0]):
            pass
        else:
            os.mkdir('kp_' + ligand[0])
            os.chdir('kp_' + ligand[0])

            local_mol = mol3D()
            local_mol.copymol3D(mol)
            deletion_indices = list(set(range(local_mol.natoms)) - set(ligand[1]))
            local_mol.deleteatoms(deletion_indices)
            local_mol.writexyz(local_name + '.xyz')

            local_infile_dict = copy.copy(infile_dict)
            local_infile_dict['name'] = local_name
            local_infile_dict['coordinates'] = local_name+'.xyz'
            local_infile_dict['charge'], local_infile_dict['spinmult'] = ligand_charge, ligand_spin
            local_infile_dict['run_type'] = 'energy'
            local_infile_dict['constraints'], local_infile_dict['convergence_thresholds'] = False, False
            local_infile_dict['machine'] = machine

            io.write_input(local_infile_dict)
            io.write_jobscript(local_name, time_limit='12:00:00', machine=machine)
            jobscripts.append(local_name + '.in')
            os.chdir('..')
    os.chdir(home)

    return jobscripts


def prep_mbe_calc(outfile_path, metal_charge=0):
    """Prep ligand breakdown.

    Parameters
    ----------
        outfile_path : str
            Path to output file.
        metal_charge : int, optional
            Charge for removed metal. Default is 0.

    Returns
    -------
        jobscripts : list
            List of jobscripts for metal binding energy jobs.

    """
    # Given a path to the outfile of a finished run, this preps the files for rigid ligand dissociation energies of all ligands
    # Returns a list of the PATH(s) to the jobscript(s) to start the rigid ligand calculations

    home = os.getcwd()
    machine = tools.get_machine()
    outfile_path = tools.convert_to_absolute_path(outfile_path)

    results = io.read_outfile(outfile_path)
    if not results['finished']:
        raise Exception('This calculation does not appear to be complete! Aborting...')

    infile_dict = io.read_infile(outfile_path)

    # special case hack for Fe test ~Freya
    fe_spin = int(infile_dict['spinmult'])
    if fe_spin == 5:
        fe_charge = 2
    else:
        fe_charge = 3
    charge = int(infile_dict['charge']) - fe_charge
    spinmult = 1  # always have a singlet

    base = os.path.split(outfile_path)[0]
    name = os.path.split(outfile_path)[-1][:-4]

    optimxyz = os.path.join(base, 'scr', 'optim.xyz')
    tools.extract_optimized_geo(optimxyz)

    breakdown_folder = os.path.join(base, name + '_mbe')

    if os.path.isdir(breakdown_folder):
        return ['Metal binding energy directory already exists']

    mol = mol3D()
    mol.readfromxyz(os.path.join(base, 'scr', 'optimized.xyz'))

    if not os.path.isdir(breakdown_folder):
        os.mkdir(breakdown_folder)
    os.chdir(breakdown_folder)

    jobscripts = []
    metal_indices = mol.findMetal()
    # Create the necessary files for the metal-removed complex single point
    mol.deleteatoms(metal_indices)
    mol.writexyz(name + '_no_metal.xyz')

    local_infile_dict = copy.copy(infile_dict)
    local_infile_dict['name'] = name + "_no_metal"
    local_infile_dict['coordinates'] = name + '_no_metal.xyz'
    local_infile_dict['charge'], local_infile_dict['spinmult'] = charge, spinmult
    local_infile_dict['run_type'] = 'energy'
    local_infile_dict['constraints'], local_infile_dict['convergence_thresholds'] = False, False
    local_infile_dict['machine'] = machine

    io.write_input(local_infile_dict)
    io.write_jobscript(name + '_no_metal', time_limit='12:00:00', machine=machine)
    jobscripts.append(name + '_no_metal.in')
    os.chdir('..')
    os.chdir(home)
    return jobscripts


def name_ligands(nested_list):
    """Takes a nested list of atom symbols and converts it to a list of unique chemical names based on the molecular formulas

    Parameters
    ----------
        nested_list : list
            List of lists of atom symbols

    Returns
    -------
        ligand_formulas : list
            List of ligand formulas.

    """
    # takes a nested list of atom symbols and converts it to a list of unique chemical names based on the molecular formulas

    def convert_to_formula(list_of_atom_symbols):
        atom_types = list(set(list_of_atom_symbols))
        atom_types.sort()

        formula = ''
        for element in atom_types:
            counter = 0
            for atom in list_of_atom_symbols:
                if element == atom:
                    counter += 1
            formula += element
            formula += str(counter)
        return formula

    ligand_formulas = list(map(convert_to_formula, nested_list))

    duplicates = []
    for i in ligand_formulas:
        number_of_duplicates = 0
        for ii in ligand_formulas:
            if i == ii:
                number_of_duplicates += 1
        duplicates.append(number_of_duplicates)

    duplication_index = 1
    for counter, i in enumerate(duplicates):
        if i > 1:
            ligand_formulas[counter] = ligand_formulas[counter] + '_' + str(duplication_index)
            duplication_index += 1

    return ligand_formulas


def check_molscontrol_log(job):
    """Check molscontrol log.

    Parameters
    ----------
        job : str
            Path to the output file.

    Returns
    -------
        killed : bool
            True if job killed.

    """
    molscontrol_logfile = "/".join(job.split("/")[:-1]) + '/molscontrol.log'
    if os.path.isfile(molscontrol_logfile):
        with open(molscontrol_logfile, "r") as fo:
            for line in fo:
                if "job killed at step" in line:
                    return True
    return False

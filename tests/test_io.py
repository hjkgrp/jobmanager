from jobmanager.io import io
from jobmanager import classes

def test_read_terachem_infile_constraint(resource_path_root):
    in_file = str(resource_path_root / "inputs" / "infiles" / "tc_constraint.in")
    in_file = classes.textfile(in_file)
    d = io.read_terachem_input(in_file)
    ref_d = {
        "run" : "minimize",
        "gpus" : "1",
        "method" : "ub3lyp",
        "dispersion" : "no",
        "basis" : "lacvps_ecp",
        "dftgrid" : "2",
        "new_minimizer" : "yes",
        "nstep" : "500",
        "scf" : "diis+a",
        "maxit" : "500",
        "watcheindiis" : "yes",
        "start_diis" : "0.001",
        "levelshift" : "yes",
        "levelshiftvala" : "0.25",
        "levelshiftvalb" : "0.25",
        "precision" : "double",
        "scrdir" : "scr",
        "timings" : "yes",
        "coordinates" : "init_geo.xyz",
        "charge" : "2",
        "spinmult" : "2",
        "$constraint_freeze" : ["bond 1_22"]
    }
    assert ref_d == d


def test_spinchargeChecker(resource_path_root):
    in_file_path = str(resource_path_root / 'inputs' / 'infiles' / 'terachem.in')
    in_file = classes.textfile(in_file_path)
    tc_dict = io.read_terachem_input(in_file) #takes file as an input
    spin_charge_bool = io.spinchargeChecker(tc_dict, in_file_path) #takes path as an input
    assert spin_charge_bool

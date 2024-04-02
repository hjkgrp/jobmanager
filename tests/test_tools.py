import jobmanager.tools as tools
from jobmanager.io import io
import os


def test_check_terminated_no_scf(resource_path_root):
    path = os.path.join(resource_path_root,'inputs/failed.out')
    results_tmp = [io.read_outfile(path, short_ouput=True)]
    results_tmp = list(zip([path], results_tmp))
    results_dict = dict()
    for outfile, tmp in results_tmp:
        results_dict[outfile] = tmp
    result = tools.check_terminated_no_scf(path,results_dict)
    reference = {'/home/davutm/jobmanager/tests/testresources/inputs/failed.out':
                 {'name': 'failed', 'charge': None, 'finalenergy': None, 'time': None, 's_squared': None,
                  's_squared_ideal': None, 'finished': False, 'min_energy': None, 'scf_error': False,
                  'thermo_grad_error': False, 'solvation_energy': None, 'optimization_cycles': None,
                  'thermo_vib_energy': None, 'thermo_vib_free_energy': None, 'thermo_suspect': False,
                  'orbital_occupation': None, 'oscillating_scf_error': False,
                  'outfile_path': '/home/davutm/jobmanager/tests/testresources/inputs/failed.out',
                  'terminated': True, 'no_scf': True}}
    assert result == reference[path]['terminated'] and result == reference[path]['no_scf']


# path = os.path.join('/home/davutm/jobmanager/tests','testresources/inputs/failed.out')
# results_tmp = [io.read_outfile(path, short_ouput=True)]
# results_tmp = list(zip([path], results_tmp))
# results_dict = dict()
# for outfile, tmp in results_tmp:
#     results_dict[outfile] = tmp
# result = tools.check_terminated_no_scf(path,results_dict)

# print(f"Assertion result: {result}")
# print(results_dict)

# outfiles = tools.find('*.out', '/home/davutm/jobmanager/tests/testresources/inputs/')
# outfiles = list(filter(tools.check_valid_outfile, outfiles))

# results_tmp = [io.read_outfile(outfile, short_ouput=True) for outfile in outfiles]
# print(results_tmp)
# results_tmp = list(zip(outfiles, results_tmp))
# print(results_tmp)
# results_dict = dict()
# for outfile, tmp in results_tmp:
#     results_dict[outfile] = tmp
# results_tmp = list(zip(path, results_tmp))
# print(results_dict)
# path = os.path.join('/home/davutm/jobmanager/tests','testresources/inputs/failed.out')

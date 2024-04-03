import jobmanager.tools as tools
from jobmanager.io import io
import os


def test_check_terminated_no_scf(resource_path_root):
    path = os.path.join(resource_path_root, "inputs/failed.out")
    results_tmp = [io.read_outfile(path, short_ouput=True)]
    results_tmp = list(zip([path], results_tmp))
    results_dict = dict()
    for outfile, tmp in results_tmp:
        results_dict[outfile] = tmp
    result = tools.check_terminated_no_scf(path, results_dict)
    assert result

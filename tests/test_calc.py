from jobmanager.calc import Calc

def test_read_input(resource_path_root):
    sample = Calc()
    sample.read_input(str(resource_path_root / 'inputs' / 'sample_input.in'))
    assert (sample.run == 'minimize' and sample.spinmult == '5' and
        sample.charge == '3' and sample.method == 'B3LYP' and
        sample.basis == 'lacvps_ecp' and
        sample.opt_dict == {'timings': 'yes', 'maxit': '500', 'scrdir': './scr', 'gpus': '1', 'new_minimizer': 'yes'})

from jobmanager.read_outfile import read_outfile

def test_read_outfile_1(resource_path_root):
    file = str(resource_path_root/'inputs'/'read_outfile'/'fifth.out')
    sample = read_outfile(file)
    dict = sample.read_from(end_line = 690)
    ref_dict = {'linenum': 689,
                'current_scf': 4,
                'currently_running': False,
                'energies': -602.4607426882}

    bool = ((dict['linenum'] == ref_dict['linenum']) and
            (dict['current_scf'] == ref_dict['current_scf']) and
            (dict['currently_running'] == ref_dict['currently_running']) and
            (dict['energies'] == ref_dict['energies']))
    assert bool

def test_read_outfile_2(resource_path_root):
    file = str(resource_path_root/'inputs'/'read_outfile'/'fifth.out')
    sample = read_outfile(file, jsonfile=file.rsplit('.',1)[0]+'.json')
    dict = sample.read_from(end_line = 695)
    ref_dict =  {'linenum': 694,
                 'current_scf': 5,
                 'currently_running': True,
                 'energies': [[]]}
    bool = ((dict['linenum'] == ref_dict['linenum']) and
            (dict['current_scf'] == ref_dict['current_scf']) and
            (dict['currently_running'] == ref_dict['currently_running']) and
            (dict['energies'] == ref_dict['energies']))
    assert bool

def test_read_outfile_3(resource_path_root):
    file = str(resource_path_root/'inputs'/'read_outfile'/'fifth.out')
    sample = read_outfile(file, jsonfile=file.rsplit('.',1)[0]+'.json')
    dict = sample.read_from(end_line = 718)
    ref_dict = {'linenum': 717,
                'current_scf': 5,
                'currently_running': True,
                'energies': [[-602.456973586, -602.2724649804,
                              -602.4452898263, -602.4615364022,
                              -602.4617684258, -602.4618023219,
                              -602.461810432, -602.4618403048,
                              -602.4618436049]]}
    bool = ((dict['linenum'] == ref_dict['linenum']) and
            (dict['current_scf'] == ref_dict['current_scf']) and
            (dict['currently_running'] == ref_dict['currently_running']) and
            (dict['energies'] == ref_dict['energies']))
    assert bool

def test_read_outfile_4(resource_path_root):
    file = str(resource_path_root/'inputs'/'read_outfile'/'fifth.out')
    sample = read_outfile(file, jsonfile=file.rsplit('.',1)[0]+'.json')
    dict = sample.read_from(end_line = 724)
    ref_dict = {'linenum': 723,
                'current_scf': 5,
                'currently_running': True,
                'energies': [[-602.456973586, -602.2724649804,
                              -602.4452898263, -602.4615364022,
                              -602.4617684258, -602.4618023219,
                              -602.461810432, -602.4618403048,
                              -602.4618436049, -602.4618444624,
                              -602.4618456718, -602.4618466261,
                              -602.4618471806, -602.4618473329,
                              -602.4618475092]]}
    bool = ((dict['linenum'] == ref_dict['linenum']) and
            (dict['current_scf'] == ref_dict['current_scf']) and
            (dict['currently_running'] == ref_dict['currently_running']) and
            (dict['energies'] == ref_dict['energies']))
    assert bool

def test_read_outfile_5(resource_path_root):
    file = str(resource_path_root/'inputs'/'read_outfile'/'fifth.out')
    sample = read_outfile(file, jsonfile=file.rsplit('.',1)[0]+'.json')
    dict = sample.read_from(end_line = 731)
    ref_dict = {'linenum': 730,
                'current_scf': 5,
                'currently_running': False,
                'energies': -602.4618478648}
    bool = ((dict['linenum'] == ref_dict['linenum']) and
            (dict['current_scf'] == ref_dict['current_scf']) and
            (dict['currently_running'] == ref_dict['currently_running']) and
            (dict['energies'] == ref_dict['energies']))
    assert bool


def test_read_outfile_6(resource_path_root):
    file = str(resource_path_root/'inputs'/'read_outfile'/'fifth.out')
    sample = read_outfile(file, jsonfile = file.rsplit('.',1)[0]+'.json')
    dict = sample.read_from(end_line = 731)
    ref_dict = {'linenum': 730,
                'current_scf': 5,
                'currently_running': False,
                'energies': -602.4618478648}
    bool = ((dict['linenum'] == ref_dict['linenum']) and
            (dict['current_scf'] == ref_dict['current_scf']) and
            (dict['currently_running'] == ref_dict['currently_running']) and
            (dict['energies'] == ref_dict['energies']))
    assert bool

def test_read_outfile_7(resource_path_root):
    file = str(resource_path_root/'inputs'/'read_outfile'/'fifth.out')
    sample = read_outfile(file, jsonfile = file.rsplit('.',1)[0]+'.json')
    dict = sample.read_from(end_line = None)
    ref_dict = {'linenum': 730,
                'current_scf': 5,
                'currently_running': False,
                'energies': -602.4618478648}
    bool = ((dict['linenum'] == ref_dict['linenum']) and
            (dict['current_scf'] == ref_dict['current_scf']) and
            (dict['currently_running'] == ref_dict['currently_running']) and
            (dict['energies'] == ref_dict['energies']))
    assert bool

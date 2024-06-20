from jobmanager.Classes.read_outfile import read_outfile

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

def test_read_outfile_8(resource_path_root):
    file = str(resource_path_root/'inputs'/'read_outfile'/'co_octahedral_3_water_6_s_5_conf_1_supercloud.out')
    sample = read_outfile(file)
    dict = sample.read_from(end_line = 641)
    ref_dict = {'linenum': 640,
                'current_scf': 4,
                'currently_running': False,
                'energies': -602.4607419824}

    bool = ((dict['linenum'] == ref_dict['linenum']) and
            (dict['current_scf'] == ref_dict['current_scf']) and
            (dict['currently_running'] == ref_dict['currently_running']) and
            (dict['energies'] == ref_dict['energies']))
    assert bool

def test_read_outfile_9(resource_path_root):
    file = str(resource_path_root/'inputs'/'read_outfile'/'co_octahedral_3_water_6_s_5_conf_1_supercloud.out')
    sample = read_outfile(file, jsonfile=file.rsplit('.',1)[0]+'.json')
    dict = sample.read_from(end_line = 646)
    ref_dict =  {'linenum': 645,
                 'current_scf': 5,
                 'currently_running': True,
                 'energies': [[]]}
    bool = ((dict['linenum'] == ref_dict['linenum']) and
            (dict['current_scf'] == ref_dict['current_scf']) and
            (dict['currently_running'] == ref_dict['currently_running']) and
            (dict['energies'] == ref_dict['energies']))
    assert bool

def test_read_outfile_10(resource_path_root):
    file = str(resource_path_root/'inputs'/'read_outfile'/'co_octahedral_3_water_6_s_5_conf_1_supercloud.out')
    sample = read_outfile(file, jsonfile=file.rsplit('.',1)[0]+'.json')
    dict = sample.read_from(end_line = 663)
    ref_dict = {'linenum': 662,
                'current_scf': 5,
                'currently_running': True,
                'energies': [[-602.4454339051, -602.4616824351,
                              -602.4612302697, -602.461789476,
                              -602.4618278078, -602.4618361543,
                              -602.4618390762, -602.4618423109,
                              -602.4618449753]]}
    bool = ((dict['linenum'] == ref_dict['linenum']) and
            (dict['current_scf'] == ref_dict['current_scf']) and
            (dict['currently_running'] == ref_dict['currently_running']) and
            (dict['energies'] == ref_dict['energies']))
    assert bool

def test_read_outfile_11(resource_path_root):
    file = str(resource_path_root/'inputs'/'read_outfile'/'co_octahedral_3_water_6_s_5_conf_1_supercloud.out')
    sample = read_outfile(file, jsonfile=file.rsplit('.',1)[0]+'.json')
    dict = sample.read_from(end_line = 669)
    ref_dict = {'linenum': 668,
                'current_scf': 5,
                'currently_running': True,
                'energies': [[-602.4454339051, -602.4616824351,
                              -602.4612302697, -602.461789476,
                              -602.4618278078, -602.4618361543,
                              -602.4618390762, -602.4618423109,
                              -602.4618449753, -602.4618459587,
                              -602.4618465099, -602.4618466348,
                              -602.4618466814, -602.4618467491,
                              -602.4618468835]]}
    bool = ((dict['linenum'] == ref_dict['linenum']) and
            (dict['current_scf'] == ref_dict['current_scf']) and
            (dict['currently_running'] == ref_dict['currently_running']) and
            (dict['energies'] == ref_dict['energies']))
    assert bool

def test_read_outfile_12(resource_path_root):
    file = str(resource_path_root/'inputs'/'read_outfile'/'co_octahedral_3_water_6_s_5_conf_1_supercloud.out')
    sample = read_outfile(file, jsonfile=file.rsplit('.',1)[0]+'.json')
    dict = sample.read_from(end_line = 674)
    ref_dict = {'linenum': 673,
                'current_scf': 5,
                'currently_running': False,
                'energies': -602.4618471062}
    bool = ((dict['linenum'] == ref_dict['linenum']) and
            (dict['current_scf'] == ref_dict['current_scf']) and
            (dict['currently_running'] == ref_dict['currently_running']) and
            (dict['energies'] == ref_dict['energies']))
    assert bool


def test_read_outfile_13(resource_path_root):
    file = str(resource_path_root/'inputs'/'read_outfile'/'co_octahedral_3_water_6_s_5_conf_1_supercloud.out')
    sample = read_outfile(file, jsonfile = file.rsplit('.',1)[0]+'.json')
    dict = sample.read_from(end_line = 674)
    ref_dict = {'linenum': 673,
                'current_scf': 5,
                'currently_running': False,
                'energies': -602.4618471062}
    bool = ((dict['linenum'] == ref_dict['linenum']) and
            (dict['current_scf'] == ref_dict['current_scf']) and
            (dict['currently_running'] == ref_dict['currently_running']) and
            (dict['energies'] == ref_dict['energies']))
    assert bool

def test_read_outfile_14(resource_path_root):
    file = str(resource_path_root/'inputs'/'read_outfile'/'co_octahedral_3_water_6_s_5_conf_1_supercloud.out')
    sample = read_outfile(file, jsonfile = file.rsplit('.',1)[0]+'.json')
    dict = sample.read_from(end_line = None)
    ref_dict = {'linenum': 3722,
                'current_scf': 33,
                'currently_running': False,
                'energies': -602.4766999337}
    bool = ((dict['linenum'] == ref_dict['linenum']) and
            (dict['current_scf'] == ref_dict['current_scf']) and
            (dict['currently_running'] == ref_dict['currently_running']) and
            (dict['energies'] == ref_dict['energies']))
    assert bool

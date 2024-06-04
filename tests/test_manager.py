from jobmanager.manager import Manager

def test_find_jobs_to_submit(resource_path_root):
    sample = Manager(str(resource_path_root/'inputs'/'jobscripts'))
    jobs_to_submit = sample.find_jobs_to_submit()
    ref_dict = {str(resource_path_root/'inputs/jobscripts/co_octahedral_3_water_6_s_5_conf_1_jobscript'):
                {'input_file': str(resource_path_root/'inputs/jobscripts/co_octahedral_3_water_6_s_5_conf_1.in'),
                 'xyz_file': str(resource_path_root/'inputs/jobscripts/co_octahedral_3_water_6_s_5_conf_1.xyz')},
                 str(resource_path_root/'inputs/jobscripts/fe_oct_2_water_6_s_1_conf_1_sse_s3_jobscript'):
                 {'input_file': str(resource_path_root/'inputs/jobscripts/fe_oct_2_water_6_s_1_conf_1_sse_s3.in')}}
    assert jobs_to_submit == ref_dict

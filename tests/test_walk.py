from jobmanager import tools


def test_find_calcs(resource_path_root):
    in_dir = str(resource_path_root / "example_project")
    xyzs = set([x for x in tools.find_calcs(in_dir)])
    expected = {'co_octahedral_3_water_6_s_5_conf_1/co_octahedral_3_water_6_s_5_conf_1.xyz', 'fe_oct_2_water_6_s_1/fe_oct_2_water_6_s_1_conf_1_sse/3/fe_oct_2_water_6_s_1_conf_1_sse_s3.xyz', 'fe_oct_2_water_6_s_1/fe_oct_2_water_6_s_1_conf_1_sse/5/fe_oct_2_water_6_s_1_conf_1_sse_s5.xyz', 'fe_oct_2_water_6_s_1/fe_oct_2_water_6_s_1_conf_1.xyz'}
    assert xyzs == expected
from jobmanager.Classes.job import Job

def test_read_jobscript(resource_path_root):
    sample = Job()
    sample.read_jobscript(str(resource_path_root/'inputs'/'jobscripts'/'orca_jobscript.sh'))
    bool = ((sample.queue_system == 'SLURM') and
            (sample.runtime=='96:00:00') and
            (sample.mem_per_cpu == '6G') and
            (sample.cpu == '1') and
            (sample.nodes == '1') and
            (sample.ntasks == '16'))


    sample_2 = Job()
    sample_2.read_jobscript(str(resource_path_root/'inputs'/'jobscripts'/'fe_oct_2_water_6_s_1_conf_1_sse_s3_jobscript'))
    bool_2 = ((sample_2.queue_system == 'SGE') and
              (sample_2.gpu == '1') and
              (sample_2.maxmemory == '8G') and
              (sample_2.runtime=='16:00:00'))
    assert bool and bool_2

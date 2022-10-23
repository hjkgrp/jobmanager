# Jobmanager
[![CI](https://github.com/hjkgrp/jobmanager/actions/workflows/CI.yaml/badge.svg)](https://github.com/hjkgrp/jobmanager/actions/workflows/CI.yaml)
[![codecov](https://codecov.io/gh/hjkgrp/jobmanager/branch/main/graph/badge.svg?token=2UCDQM8BJ0)](https://codecov.io/gh/hjkgrp/jobmanager)

## Installation

We currently recommend installation via the [Conda](https://conda.io/docs/) package management system.
1. Prerequisite: have [Anaconda or miniconda](https://www.anaconda.com/distribution/) installed on your system. **For M1 Macs, please use [Miniforge](https://github.com/conda-forge/miniforge) for Mac OSX arm64.** (We do not recommend simultaneously installing Anaconda and Miniforge - only install Miniforge.)

2. Clone source from github.

   ```bash
   git clone https://github.com/hjkgrp/jobmanager.git
   ```

3. Go to the folder root folder for jobmanager, create the conda environment from the yaml file (`tools/conda_envs/jobmanager.yaml`). You can specify a different name for this environment at the first line of the yaml file. If you would like to use the most recent version of molSimplify, delete the line with `molsimplify` in the yaml file and follow the [molSimplify installation instructions](https://github.com/hjkgrp/molSimplify#installation)

   ```bash
   cd jobmanager/tools/conda_envs
   conda env create -f jobmanager.yaml
   ```

## Usage
Job manager reads the configure at the job run time as well as the root directory.

Must be in the directory of the configure, and then call `jobmanager` -- this will follow the configure file and execute on the working tree as described below.

The folder hierarchy is as follows:
```
base_directory/
   |--configure
   |--unique_name1/
   |--unique_name2/
   |--unique_name3/
      |--unique_name3.in
      |--unique_name3.xyz
      |--unique_name3_jobscript
```
*Note*: The output file generated via the jobscript must also follow the naming convention (i.e. must be `unique_name3.out` for the above example).

### Configuration
The configure file (labeled `configure`, NOT `configure.txt`) contains all of the info about what needs to happen to the runs.

The configure file can live within subdirectories, if derivative jobs are only necessary for a subset of jobs.

It can be modified on the fly to add info (if you want to add thermo or solvent, or HFX resample). 

If the geo_check:oct is in the configure file, then jobs that are completed are checked for good geometries before spawning any derivative jobs.

If job_recovery is not desired, that line should not be there. 

### In Python scripts
If loading the module (in a python script), do the following:

```python
import jobmanager.resub as resub
resub.main()
```

### Notes on running
The job manager uses the unique name as a queue identifier. Thus, it will not submit a job with a specific unique name, if that unique name is already in the queue. This will prevent the same job from being resubmitted while it is still running. 

Killing the job manager by using ctrl-c should be ok to stop the manager from running. By default, it cycles with a sleep period of 2 hr.
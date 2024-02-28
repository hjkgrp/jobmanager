import os
from jobmanager.psi4_utils.run import run_bash


def test_run_bash(tmpdir):
    command = "touch test.txt"
    rundir = os.mkdir(tmpdir / "run")
    run_bash(command, tmpdir, rundir)
    # Check that the command "touch test.txt" was executed
    assert os.path.isfile(rundir / "test.txt")
    # Check that all files were copied
    for file in ["loop_run.py", "loop_rescue.py", "loop_derivative_jobs.py"]:
        assert os.path.isfile(rundir / file)

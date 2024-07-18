import os
import psi4
from jobmanager.psi4_utils.run_utils import RunUtils


def test_psi4():
    psi4.set_memory('500 MB')
    h2o = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

    e = psi4.energy('scf/cc-pvdz')

    assert abs(e + 76.0266327351) < 0.00001

def test_run_bash(tmpdir):
    cwd = os.getcwd()
    command = "touch test.txt"
    rundir = tmpdir / "run"
    # Create the run dir
    rundir.mkdir()
    os.chdir(rundir)
    # Actually execute run_bash
    run_utils = RunUtils()
    run_utils.run_bash(command, tmpdir, rundir)
    os.chdir(cwd)
    # Check that the command "touch test.txt" was executed
    assert os.path.isfile(rundir / "test.txt")

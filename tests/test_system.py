from jobmanager.Classes.system import System
from molSimplify.Scripts.rmsd import align_rmsd


def test_build(resource_path_root):
    s1 = System(spinmult=5, charge=2)
    s1.set_xyzpath(resource_path_root / 'fe_2_h2o_6_s_5.xyz')
    s1.build_from_xyz()
    s2 = System(spinmult=5, charge=2)
    s2.build_from_molS("fe",["h2o"]*6,"oct")
    assert align_rmsd(s1.molecule, s2.molecule) < 1e-5

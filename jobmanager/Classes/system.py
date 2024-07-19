from molSimplify.Classes.mol3D import mol3D
from molSimplify.Scripts.generator import startgen_pythonic
import os

class System:
    """
    Interface with molSimplify to generate and check geometries.
    Contains information about a chemical system (e.g., spin, charge, geometry, metal, ligands)
    """
    def __init__(self,spinmult=None, charge=None):
        if spinmult is None:
            self.spinmult = 1
        else:
            self.spinmult = spinmult
        if charge is None:
            self.charge = 0
        else:
            self.charge = charge
        self.xyz = None
        self.molecule = None


    def build_from_molS(self, metal, ligands, geometry, ligand_charges = None):
        if ligand_charges is None:
            metal_ox = self.charge
        else:
            metal_ox = self.charge - sum(ligand_charges)
        input_dict = {
            '-core':metal.lower(),
            '-lig':','.join(ligands),
            '-geometry':geometry,
            '-oxstate':str(metal_ox),
            '-spin':str(self.spinmult)
        }
        files, emsg, rundiag = startgen_pythonic(input_dict)
        self.molecule = rundiag.mol
    
    
    def set_xyzpath(self, path):
        self.xyz = path
    
    def build_from_xyz(self):
        if not (self.xyz and os.path.exists(self.xyz)):
            raise ValueError(f'{self.xyz} is not a valid filename')
        self.molecule = mol3D()
        self.molecule.readfromxyz(self.xyz)
        
import psi4
import os
import json
import shutil
import numpy as np
from jobmanager.psi4_utils.run_utils import RunUtils
from jobmanager.io.molden import load_molden
from jobmanager.psi4_utils.TCtoPsi4 import TCtoPsi4

class Psi4Utils:
    """
    Contains functions useful for interfacing with Psi4.
    """

    def __init__(self, psi4_config):
        self.config = psi4_config

    def lacvps(self, mol, role):
        """
        LACVP* basis set: H-Ar are 6-31G*, rest are LANL2DZ.

        mol is a psi4.core.Molecule object

        Sets the basis set for all atoms to 6-31g* and then
        sets later elements to lanl2dz.

        In the format required for psi4.qcdb.libmintsbasisset.basishorde

        Returns an empty dictionary since setting the basis to LACVP* does
        not introduce any new basis set that a basis string would be needed for
        """
        basstrings = {}
        mol.set_basis_all_atoms("6-31g*", role=role)
        mol.set_basis_by_symbol("fe", "lanl2dz", role=role)
        mol.set_basis_by_symbol("co", "lanl2dz", role=role)
        mol.set_basis_by_symbol("cr", "lanl2dz", role=role)
        mol.set_basis_by_symbol("mn", "lanl2dz", role=role)
        mol.set_basis_by_symbol("mo", "lanl2dz", role=role)
        mol.set_basis_by_symbol("tc", "lanl2dz", role=role)
        mol.set_basis_by_symbol("ru", "lanl2dz", role=role)
        mol.set_basis_by_symbol("rh", "lanl2dz", role=role)
        mol.set_basis_by_symbol("I", "lanl2dz", role=role)
        mol.set_basis_by_symbol("Br", "lanl2dz", role=role)
        mol.set_basis_by_symbol("hf", "lanl2dz", role=role)
        return basstrings
    
    def get_molecule(self, xyzfile, charge, spin, sym='c1'):
        """
        From an xyzfile, charge, and spin, construct a Psi4 molecule object.
        If multiple structures present in the xyz file, will read the first one.
        """
        #Psi4 expects charge, spin as 2 integers on first line
        wholetext = "%s %s\n" % (charge, spin)
        #append the coordinates in the xyz to wholetext
        if os.path.isfile(xyzfile):
            with open(xyzfile, "r") as fo:
                natoms = int(fo.readline().split()[0])
                for ii in range(natoms):
                    wholetext += fo.readline()
        else:
            raise ValueError("Provide a file!")
        #add the symmetry, tell Psi4 to not reordient or recenter
        wholetext += "\nsymmetry %s\nnoreorient\nnocom\n" % sym
        #build the Psi4 geometry
        mol = psi4.geometry("""%s""" % wholetext)
        #store the molecule in the class for reference later
        self.mol = mol
        return mol
    
    def setup_dft_parameters(self):
        """
        Sets up the Psi4 parameters before the DFT calculation is run,
        according to the psi4_config file and some default values.
        """
        psi4_config = self.config
        psi4.set_memory(psi4_config["memory"])
        psi4.set_num_threads(psi4_config["num_threads"])
        #If LACVP*, define the custom mix of basis sets
        #and run in Cartesian coordinates
        if psi4_config["basis"] == "lacvps":
            psi4.qcdb.libmintsbasisset.basishorde['LACVPS'] = self.lacvps
            psi4.set_options({"puream": False})
        #If a def2 basis set, run in spherical coordinates
        elif psi4_config["basis"] == "def2-tzvp" or 'def2' in psi4_config["basis"]:
            psi4.set_options({"puream": True})
        else:
            #psi4.set_options({"puream": True})
            raise ValueError("Only lacvps and def2 basis sets are supported!")
        '''
        The following options do:
        reference: whether to run restricted or unrestricted HF or KS
        DF_SCF_GUESS: does density fitting to converge orbitals before doing exact integrals in SCF
        -turned off to ensure maximal correspondence with the molden initial guess
        scf_type: algorithm to use for SCF computation.
        -recommended algorithm is DF, density-fitted algorithm
        dft_pruning_scheme: controls pruning of quadrature grid
        -recommended setting is robust
        DFT_BASIS_TOLERANCE: determines cutoff radius for each shell of basis functions
        -default is 1e-12, but can be increased to speed up, minimal accuracy penalty
        INTS_TOLERANCE: threshold for screening method, below which two-electron integrals ignored
        -increased to 1e-10 from 1e-12
        PRINT_MOS: whether or not to print molecular orbitals
        dft_spherical_points: number of spherical points for quadrature
        -increased from 302
        sft_radial_points: number of radial points for quadrature
        -increased from 75
        guess: set to read the orbitals from a wfn file
        '''
        psi4.set_options({
            'reference': psi4_config["ref"],
            "DF_SCF_GUESS": False,
            "scf_type": "df",
            "dft_pruning_scheme": "robust",
            "basis": psi4_config["basis"],
            "DFT_BASIS_TOLERANCE": 1e-10,
            "INTS_TOLERANCE": 1e-10,
            "PRINT_MOS": False,
            "dft_spherical_points": 590,
            "dft_radial_points": 99,
            # "SOSCF": True,
            # "DAMPING_PERCENTAGE": 20,
            # "BASIS_GUESS": True,
            "guess": "read", })
        
    def b3lyp_hfx(self, hfx):
        """
        Given a percentage HFX, returns a string that can specify a custom B3LYP with
        that HFX.

        B3LYP originally uses 20% HFX with the 80% non-HF exchange being split between
        0.08 LDA and 0.72 B88. This retains the 90:10 split between B and LDA while allowing
        for scaling of the absolute non-HF exchange fraction to values != 20.
        Retains the correlation mixing, 81% LYP and 19% VWN.
        """
        hfx_func = {
            "name": "hfx_func",
            "x_functionals": {
                "GGA_X_B88": {"alpha": 0.9*(1-hfx*0.01)},
                "LDA_X": {"alpha": 0.1*(1-hfx*0.01)}
                    },
            "x_hf": {"alpha": hfx*0.01},
            "c_functionals": {
                "GGA_C_LYP": {"alpha": 0.81},
                "LDA_C_VWN_RPA": {"alpha": 0.19}
            }
        }
        return hfx_func
    

    def run_b3lyp(self, rundir="./b3lyp", return_wfn=True,
                  psi4_scr = './', filename='output'):
        """
        Runs a Psi4 single point calculation using B3LYP in Psi4, from a TC molden.

        Will write three files:
        -wfn-1step.180.npy: wfn file for a 1-step SCF, used only to get a wfn object.
        -wfn-1step-tc.180.npy: wfn file populated with the coefficients from the molden
        -wfn.180.npy: wfn file that has been converged in Psi4 with the requested parameters,
        starting from the TC initial guess.

        Results will be stored in rundir as a subdirectory of the current directory.
        """
        psi4_config = self.config
        #update the Psi4 config with the charge, spin, and method
        #in the charge_spin_info json file.
        with open(psi4_config["charge-spin-info"], "r") as f:
            d = json.load(f)
        psi4_config.update(d)
        #Ensure that the run directory exists
        run_utils = RunUtils()
        run_utils.ensure_dir(rundir)
        #set up the molecule and calculation parameters
        sym = 'c1' if 'sym' not in psi4_config else psi4_config['sym']
        mol = self.get_molecule(psi4_config["xyzfile"], psi4_config["charge"], psi4_config["spin"], sym)
        self.setup_dft_parameters(psi4_config)
        if os.path.isfile(psi4_config["moldenfile"]):
            #logic to transfer the TC coefficients to Psi4
            psi4.core.set_output_file(rundir + '/' + filename + '.dat', False)
            # 1-step SCF
            #run very crudely, only goal is to get the wfn object, will later populate with the molden results
            psi4.set_options({
                "maxiter": 5,
                "D_CONVERGENCE": 1e5,
                "E_CONVERGENCE": 1e5,
                "fail_on_maxiter": False})
            #If def2-TZVP, converge in def2-SV(P) first, since expecting to project up from def2-SV(P)
            if psi4_config["basis"] == "def2-tzvp":
                psi4.set_options({"basis": "def2-sv(p)"})
            #If the initial calculation was done with B3LYP and a HFX not equal to 20, use functional defined in b3lyp_hfx
            if "b3lyp_hfx" in psi4_config and psi4_config["b3lyp_hfx"] != 20:
                print("customized b3lyp with different HFX: ", psi4_config["b3lyp_hfx"])
                e, wfn = psi4.energy("scf", dft_functional=self.b3lyp_hfx(psi4_config["b3lyp_hfx"]),  molecule=mol, return_wfn=True)
            else:
                e, wfn = psi4.energy('b3lyp', molecule=mol, return_wfn=True)
            wfn.to_file(rundir + "/wfn-1step.180")
            # Get converged WFN
            d_molden = load_molden(psi4_config["moldenfile"])
            #if calculation restricted or not, checking for rks or rhf references
            restricted = True if any(x in psi4_config["ref"] for x in ["r", "R"]) else False
            #initialize TCtoPsi4 class and use to convert molden to wfn
            if psi4_config["basis"] == "lacvps" or '6-31g' in psi4_config['basis']:
                #currently only support Cartesian calculation in Psi4 with 6-31g*
                converter = TCtoPsi4(psi4_config['basis'], 'Cartesian')
                Ca, Cb, mapping = converter.tcmolden2psi4wfn_ao_mapping(d_molden, restricted=restricted)
            elif 'def2' in psi4_config['basis']:
                #currrently only support Spherical calculation in Psi4 with def2 basis sets
                converter = TCtoPsi4(psi4_config['basis'], 'Spherical')
                Ca, Cb, mapping = converter.tcmolden2psi4wfn_ao_mapping(d_molden, restricted=restricted)
            else:
                raise ValueError('Unsuppoorted Basis Set encountered! Please use 6-31G, LACVP*, or def2 basis sets.')
            #populate the 1-step SCF wfn with the coefficients from the molden
            wfn_minimal_np = np.load(rundir + "/wfn-1step.180.npy", allow_pickle=True)
            wfn_minimal_np[()]['matrix']["Ca"] = Ca
            if not restricted:
                wfn_minimal_np[()]['matrix']["Cb"] = Cb
            else:
                wfn_minimal_np[()]['matrix']["Cb"] = Ca
            np.save(rundir + "/wfn-1step-tc.180.npy", wfn_minimal_np)
            # Copy wfn file to the right place with a right name
            #Psi4 requires format /scratch/output.moleculename.pid.180.npy
            pid = str(os.getpid())
            targetfile = psi4_scr + filename + '.default.' + pid + '.180.npy'
            shutil.copyfile(rundir + "/wfn-1step-tc.180.npy", targetfile)
            # Final scf---
            #only allow for 50 iterations since the initial guess should converge quickly
            psi4.set_options({
                "maxiter": 50,
                "D_CONVERGENCE": 3e-5,
                "E_CONVERGENCE": 3e-5,
                "fail_on_maxiter": True})
        else:
            #if no molden provided, allow for more iterations to ensure convergence
            psi4.core.set_output_file(rundir + '/' + filename + '.dat', False)
            print("Warning: no Molden file is used to initialize this calculation!")
            psi4.set_options({
                "maxiter": 250 if "maxiter" not in psi4_config else psi4_config["maxiter"],
                # "guess": "GWH",
                "D_CONVERGENCE": 3e-5,
                "E_CONVERGENCE": 3e-5,
                "fail_on_maxiter": True})
            if psi4_config["basis"] == "def2-tzvp":
                psi4.set_options({"basis": "def2-sv(p)"})
        sucess = False
        try:
            #final SCF calculation
            if "b3lyp_hfx" in psi4_config and psi4_config["b3lyp_hfx"] != 20:
                #if calculation was done with a HFX different from 20, use custom B3LYP
                print("customized b3lyp with different HFX: ", psi4_config["b3lyp_hfx"])
                e, wfn = psi4.energy("scf", self.b3lyp_hfx(psi4_config["b3lyp_hfx"]),  molecule=mol, return_wfn=True)
            else:
                e, wfn = psi4.energy('b3lyp', molecule=mol, return_wfn=True)
            wfn.to_file(rundir + "/wfn.180")
            sucess = True
        except:
            print("This calculation does not converge.")
        #If using def2-TZVP, did the first calculation in def2-SV(P), so do another calculation to get TZVP.
        if psi4_config["basis"] == "def2-tzvp" and sucess:
            #allow for more iterations since doing a projection
            #TODO: seems to not use the converged sv(p) guess, write over the TC guess with the Psi4 result?
            psi4.set_options({"basis": "def2-tzvp", "maxiter": 200 if "maxiter" not in psi4_config else psi4_config["maxiter"]})
            try:
                if "b3lyp_hfx" in psi4_config and psi4_config["b3lyp_hfx"] != 20:
                    print("customized b3lyp with different HFX: ", psi4_config["b3lyp_hfx"])
                    e, wfn = psi4.energy("scf", self.b3lyp_hfx(psi4_config["b3lyp_hfx"]),  molecule=mol, return_wfn=True)
                else:
                    e, wfn = psi4.energy('b3lyp', molecule=mol, return_wfn=True)
                wfn.to_file("wfn.180")
            except:
                print("This calculation does not converge.")
        success = run_utils.check_sucess()
        #remove copied wfn file, psi4 log files
        for filename in os.listdir('./'):
            if ("psi." in filename) or ("default" in filename):
                print("removing: :", filename)
                os.remove(filename)
        return success
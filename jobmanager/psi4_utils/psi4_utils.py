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

        Sets the basis set for all atoms to 6-31g* and then
        sets later elements to lanl2dz.

        Parameters:
            mol: psi4.core.Molecule object
                Molecule you are setting the basis set for
            role: (?)
                Required for psi4.qcdb.libmintsbasisset.basishorde format
        Outputs:
            bassstrings: dict
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

        Parameters:
            xyzfile: str
                Path to the xyz file used to construct the molecule.
            charge, spin: int
                Charge and spin multiplicity of the molecule.
            sym: str
                Symmetry of the molecule.
        Outputs:
            mol: psi4.geometry
                Object that the Psi4 calculations can be performed on.
        """
        #Psi4 expects charge, spin as 2 integers on first line
        wholetext = "%s %s\n" % (charge, spin)
        #append the coordinates in the xyz to wholetext
        if os.path.isfile(xyzfile):
            with open(xyzfile, "r") as fo:
                natoms = int(fo.readline().split()[0])
                #skip the comment line, append all remaining lines
                fo.readline()
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

    def get_hfx_functional(self, functional, hfx):
        """
        Returns a functional with adjusted HFX percentage.

        Parameters:
            functional: str
                Functional one wants to adjust the HFX for.
            hfx: int
                Percentage HFX one wants to use in the functional.
        Outputs:
            hfx_func: dict
                Functional with adjusted HFX in the format desired for Psi4.
        """
        fmap = {"tpss": "TPSS", "scan": "SCAN", "m06-l": "M06_L", "mn15-l": "MN15_L"}
        if functional == "bp86":
            hfx_func = {
                "name": "hfx_func",
                "x_functionals": {"GGA_X_B88": {"alpha": 1-hfx*0.01}},
                "x_hf": {"alpha": hfx*0.01},
                "c_functionals": {"GGA_C_P86": {}}
            }
        elif functional == "blyp":
            hfx_func = {
                "name": "hfx_func",
                "x_functionals": {"GGA_X_B88": {"alpha": 1-hfx*0.01}},
                "x_hf": {"alpha": hfx*0.01},
                "c_functionals": {"GGA_C_LYP": {}}
            }
        elif functional == "b3lyp":
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
        elif functional == "pbe":
            hfx_func = {
                "name": "hfx_func",
                "x_functionals": {"GGA_X_PBE": {"alpha": 1-hfx*0.01}},
                "x_hf": {"alpha": hfx*0.01},
                "c_functionals": {"GGA_C_PBE": {}}
            }
        elif functional in ["m06-l", "mn15-l", "scan", "tpss"]:
            mega = "" if "PBE" in functional else "M"
            hfx_func = {
                "name": "hfx_func",
                "x_functionals": {"%sGGA_X_%s" % (mega, fmap[functional]): {"alpha": 1-hfx*0.01}},
                "x_hf": {"alpha": hfx*0.01},
                "c_functionals": {"%sGGA_C_%s" % (mega, fmap[functional]): {}}
            }
        else:
            raise ValueError("This functional has not been implemented with HFX resampling yet: ", functional)
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

        Parameters:
            rundir: str
                Path where the B3LYP calculation results and wavefunctions will be stored.
            return_wfn: bool
                Whether or not the wfn file should be written after the calculation.
            psi4_scr: str
                Path of the Psi4 scratch directory.
            filename: str
                Name out the output .dat file containing results.
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
        self.setup_dft_parameters()
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
            #If the initial calculation was done with B3LYP and a HFX not equal to 20, use functional defined in get_hfx_functional
            if "b3lyp_hfx" in psi4_config and psi4_config["b3lyp_hfx"] != 20:
                print("customized b3lyp with different HFX: ", psi4_config["b3lyp_hfx"])
                e, wfn = psi4.energy("scf", dft_functional=self.get_hfx_functional('b3lyp', psi4_config["b3lyp_hfx"]),  molecule=mol, return_wfn=True)
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
                e, wfn = psi4.energy("scf", self.get_hfx_functional('b3lyp', psi4_config["b3lyp_hfx"]),  molecule=mol, return_wfn=True)
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
                    e, wfn = psi4.energy("scf", self.get_hfx_functional('b3lyp', psi4_config["b3lyp_hfx"]),  molecule=mol, return_wfn=True)
                else:
                    e, wfn = psi4.energy('b3lyp', molecule=mol, return_wfn=True)
                wfn.to_file(rundir + "/wfn.180")
            except:
                print("This calculation does not converge.")
        #Check if the B3LYP calculation succeeded
        success = run_utils.check_sucess(path=rundir)
        #remove copied wfn file, psi4 log files
        for filename in os.listdir('./'):
            if ("psi." in filename) or ("default" in filename):
                print("removing: :", filename)
                os.remove(filename)
        return success

    def run_general(self, functional="b3lyp", rundir="./", return_wfn=False,
                    psi4_scr='./', filename='output'):
        """
        From a directory, launches calculations with other functionals from the .wfn specified in the functional input.
        Does so in subdirectories with names corresponding to the functional names.

        Parameters:
            functional: str
                Name of the functional one wants to run.
            rundir: str
                Path where all functional calculations are run from.
                Results will be stored in the path rundir+functional.
            return_wfn: bool
                Whether or not the wfn file should be written after the calculation.
            psi4_scr: str
                Path of the Psi4 scratch directory.
            filename: str
                Name out the output .dat file containing results.
        """
        #Make the subdirectory, load relevant information
        psi4_config = self.config
        rundir = rundir + functional.replace("(", "l-").replace(")", "-r")
        with open(psi4_config["charge-spin-info"], "r") as f:
            d = json.load(f)
        psi4_config.update(d)
        run_utils = RunUtils()
        run_utils.ensure_dir(rundir)

        #Set up Psi4 parameters
        psi4.core.set_output_file(rundir + '/' + filename + '.dat', False)
        sym = 'c1' if 'sym' not in psi4_config else psi4_config['sym']
        mol = self.get_molecule(psi4_config["xyzfile"], psi4_config["charge"], psi4_config["spin"], sym)
        self.setup_dft_parameters()

        # Copy wfn file to the right place with a right name---
        pid = str(os.getpid())
        targetfile = psi4_scr + filename + '.default.' + pid + '.180.npy'
        if os.path.isfile(psi4_config["wfnfile"]):
            shutil.copyfile(psi4_config["wfnfile"], targetfile)
        else:
            print("wfn file not found! Check that the proper path is specified.")

        # Final scf---
        psi4.set_options({
            "maxiter": 50 if "maxiter" not in psi4_config else psi4_config["maxiter"],
            "D_CONVERGENCE": 3e-5,
            "E_CONVERGENCE": 3e-5,
            "fail_on_maxiter": True})

        if not (("ccsd" in functional) or ("mp2" in functional) or ("scf" in functional)):
            try:
                if "hfx_" not in functional:
                    #Run Psi4 calculation without HFX adjustment
                    e, wfn = psi4.energy(functional, molecule=mol, return_wfn=True)
                else:
                    #Define custom functional with adjusted HFX
                    basefunc, hfx = functional.split("_")[0], int(functional.split("_")[-1])
                    print("HFX sampling: ", basefunc, hfx)
                    e, wfn = psi4.energy("scf", dft_functional=self.get_hfx_functional(basefunc, hfx),  molecule=mol, return_wfn=True)
                if return_wfn:
                    wfn.to_file(rundir + "/wfn.180")
            except:
                print("This calculation does not converge.")
        else:
            #If a CC, MP2, or HF (SCF) calculation
            #Does not use the molden as the initial guess
            print("running CC: ", functional)
            psi4.set_options({
                'reference': d['ref'].replace("ks", "hf"),
                'R_CONVERGENCE': 1e-5,
                'E_CONVERGENCE': 5e-5,
                'D_CONVERGENCE': 5e-5,
                "mp2_type": "df",
                "cc_type": "conv",
                "scf_type": "df",
                'nat_orbs': True,
                'FREEZE_CORE': True,
                "GUESS": "SAD",
                })
            e, wfn = psi4.energy(functional, molecule=mol, return_wfn=True)
            if return_wfn:
                wfn.to_file(rundir + "/wfn.180")
        #Check success, Remove temporary files
        success = run_utils.check_sucess(path=rundir)
        for filename in os.listdir("./"):
            if ("psi." in filename) or ("default" in filename):
                print("removing: :", filename)
                os.remove(filename)
        return success

    def run_general_hfx(self, functional, hfx, wfn, return_wfn=True,
                        psi4_scr='./', filename='output', rundir='./'):
        """
        From a converged calculation (given in wfn), calculates the energy using
        functional and hfx.

        Parameters:
            functional: str
                Name of the functional one wants to run.
            hfx: int
                Percentage of HFX to use in the calculation.
            wfn: str
                Path to the .wfn file used for reference.
            return_wfn: bool
                Whether or not the wfn file should be written after the calculation.
            psi4_scr: str
                Path of the Psi4 scratch directory.
            filename: str
                Name out the output .dat file containing results.
            rundir: str
                Path where all functional calculations are run from.
                Results will be stored in the path rundir+functional-hfx.
        """
        #Set up folders and configuration
        psi4_config = self.config
        rundir = rundir + functional + "-%d" % hfx
        with open(psi4_config["charge-spin-info"], "r") as f:
            d = json.load(f)
        psi4_config.update(d)
        run_utils = RunUtils()
        run_utils.ensure_dir(rundir)

        #Set up parameters
        psi4.core.set_output_file(rundir + '/' + filename + '.dat', False)
        sym = 'c1' if 'sym' not in psi4_config else psi4_config['sym']
        mol = self.get_molecule(psi4_config["xyzfile"], psi4_config["charge"], psi4_config["spin"], sym)
        self.setup_dft_parameters()

        # Copy wfn file to the right place with a right name---
        pid = str(os.getpid())
        targetfile = psi4_scr + filename + '.default.' + pid + '.180.npy'
        if not os.path.isfile(wfn):
            #If the referenced wfn is not found, skip
            print("Previous calculation failed... This one is skipped.")
            return False
        shutil.copyfile(wfn, targetfile)
        # Final scf---
        psi4.set_options({
            "maxiter": 50 if "maxiter" not in psi4_config else psi4_config["maxiter"],
            "D_CONVERGENCE": 3e-5,
            "E_CONVERGENCE": 3e-5,
            "fail_on_maxiter": True})
        try:
            e, wfn_o = psi4.energy("scf", molecule=mol, return_wfn=True, dft_functional=self.get_hfx_functional(functional, hfx))
            if return_wfn:
                wfn_o.to_file(rundir + "/wfn.180")
            # os.remove(wfn)
        except:
            print("This calculation does not converge.")
        success = run_utils.check_sucess(path=rundir)
        #remove temporary files
        for filename in os.listdir("./"):
            if ("psi." in filename) or ("default" in filename):
                print("removing: :", filename)
                os.remove(filename)
        return success
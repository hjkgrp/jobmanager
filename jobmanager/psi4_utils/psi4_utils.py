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
        """
        Stores the dictionary of options for the Psi4 workflow.
        """

        #define the base functional that other calculations will be started from
        if "base_functional" not in psi4_config:
            #defaults to B3LYP
            psi4_config["base_functional"] = 'b3lyp'
        else:
            #remove parentheses from functional names
            functional = psi4_config["base_functional"]
            psi4_config["base_functional"] = functional.replace("(", "l-").replace(")", "-r")
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

        Supported functionals: bp86, blyp, b3lyp, pbe, m06-l, mn15-l, scan, tpss

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


    def run_initial(self, rundir="./", return_wfn=True,
                  psi4_scr = './', filename='output'):
        """
        Runs a Psi4 single point calculation in Psi4, starting from a TC molden.
        Will use the functional in psi4_config['base_functional'] if provided, otherwise B3LYP.

        Will write three files:
        -wfn-1step.180.npy: wfn file for a 1-step SCF, used only to get a wfn object.
        -wfn-1step-tc.180.npy: wfn file populated with the coefficients from the TC molden.
        -wfn.180.npy: wfn file that has been converged in Psi4 with the requested parameters,
        starting from the TC initial guess.

        Results will be stored in 'rundir/' + psi4_config['base_functional'] .

        Parameters:
            rundir: str
                Path where the calculation will be run from.
                Results and wavefunctions will be stored in a subdirectory with the name in psi4_config['base_functional'].
            return_wfn: bool
                Whether or not the wfn file should be written after the calculation.
            psi4_scr: str
                Path of the Psi4 scratch directory.
            filename: str
                Name out the output .dat file containing results.
        """
        psi4_config = self.config
        base_func = psi4_config['base_functional']
        #update the Psi4 config with the charge, spin, and method
        #in the charge_spin_info json file.
        with open(psi4_config["charge-spin-info"], "r") as f:
            d = json.load(f)
        psi4_config.update(d)
        #Ensure that the run directory exists
        run_utils = RunUtils()
        run_utils.ensure_dir(rundir + base_func)
        #set up the molecule and calculation parameters
        sym = 'c1' if 'sym' not in psi4_config else psi4_config['sym']
        mol = self.get_molecule(psi4_config["xyzfile"], psi4_config["charge"], psi4_config["spin"], sym)
        self.setup_dft_parameters()
        if os.path.isfile(psi4_config["moldenfile"]):
            #logic to transfer the TC coefficients to Psi4
            psi4.core.set_output_file(rundir + base_func + '/' + filename + '.dat', False)
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
            '''
            #Deprecated from when only B3LYP was supported as a base functional, updated to now include other functionals
            #If the initial calculation was done with B3LYP and a HFX not equal to 20, use functional defined in get_hfx_functional
            if "b3lyp_hfx" in psi4_config and psi4_config["b3lyp_hfx"] != 20:
                print("customized b3lyp with different HFX: ", psi4_config["b3lyp_hfx"])
                e, wfn = psi4.energy("scf", dft_functional=self.get_hfx_functional('b3lyp', psi4_config["b3lyp_hfx"]),  molecule=mol, return_wfn=True)
            else:
                e, wfn = psi4.energy('b3lyp', molecule=mol, return_wfn=True)
            wfn.to_file(rundir + "b3lyp/wfn-1step.180")
            '''
            if "_hfx_" in base_func:
                #if hfx in functional, expects in the form functional_hfx_XX, where XX is the percent of HFX included.
                func = base_func.split('_')[0]
                hfx = base_func.split('_')[-1]
                e, wfn = psi4.energy("scf", dft_functional=self.get_hfx_functional(func, hfx),  molecule=mol, return_wfn=True)
            else:
                e, wfn = psi4.energy(base_func, molecule=mol, return_wfn=True)
            wfn.to_file(rundir + base_func + "/wfn-1step.180")
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
            wfn_minimal_np = np.load(rundir + base_func + "/wfn-1step.180.npy", allow_pickle=True)
            wfn_minimal_np[()]['matrix']["Ca"] = Ca
            if not restricted:
                wfn_minimal_np[()]['matrix']["Cb"] = Cb
            else:
                wfn_minimal_np[()]['matrix']["Cb"] = Ca
            np.save(rundir + base_func + "/wfn-1step-tc.180.npy", wfn_minimal_np)
            # Copy wfn file to the right place with a right name
            #Psi4 requires format /scratch/output.moleculename.pid.180.npy
            pid = str(os.getpid())
            targetfile = psi4_scr + filename + '.default.' + pid + '.180.npy'
            shutil.copyfile(rundir + base_func + "/wfn-1step-tc.180.npy", targetfile)
            # Final scf---
            #only allow for 50 iterations since the initial guess should converge quickly
            psi4.set_options({
                "maxiter": 50,
                "D_CONVERGENCE": 3e-5,
                "E_CONVERGENCE": 3e-5,
                "fail_on_maxiter": True})
        else:
            #if no molden provided, allow for more iterations to ensure convergence
            psi4.core.set_output_file(rundir + base_func + '/' + filename + '.dat', False)
            print("Warning: no Molden file is used to initialize this calculation!")
            psi4.set_options({
                "maxiter": 250 if "maxiter" not in psi4_config else psi4_config["maxiter"],
                # "guess": "GWH",
                "D_CONVERGENCE": 3e-5,
                "E_CONVERGENCE": 3e-5,
                "fail_on_maxiter": True})
            if psi4_config["basis"] == "def2-tzvp":
                psi4.set_options({"basis": "def2-sv(p)"})
        success = False
        try:
            #final SCF calculation
            '''
            #Deprecated, see above deprecated section
            if "b3lyp_hfx" in psi4_config and psi4_config["b3lyp_hfx"] != 20:
                #if calculation was done with a HFX different from 20, use custom B3LYP
                print("customized b3lyp with different HFX: ", psi4_config["b3lyp_hfx"])
                e, wfn = psi4.energy("scf", self.get_hfx_functional('b3lyp', psi4_config["b3lyp_hfx"]),  molecule=mol, return_wfn=True)
            else:
                e, wfn = psi4.energy('b3lyp', molecule=mol, return_wfn=True)
            wfn.to_file(rundir + "b3lyp/wfn.180")
            '''
            if "_hfx_" in base_func:
                #if hfx in functional, expects in the form functional_hfx_XX, where XX is the percent of HFX included.
                func = base_func.split('_')[0]
                hfx = base_func.split('_')[-1]
                e, wfn = psi4.energy("scf", dft_functional=self.get_hfx_functional(func, hfx),  molecule=mol, return_wfn=True)
            else:
                e, wfn = psi4.energy(base_func, molecule=mol, return_wfn=True)
            wfn.to_file(rundir + base_func + "/wfn.180")
            success = True
        except:
            print("This calculation does not converge.")
        #If using def2-TZVP, did the first calculation in def2-SV(P), so do another calculation to get TZVP.
        if psi4_config["basis"] == "def2-tzvp" and success:
            #Project to def2-TZVP; allow for more iterations since doing a projection
            
            #copy the result of the def2-SV(P) calculation to Psi4's expected location
            #check if this is redundant, i.e., if Psi4 saves the wfn of the last calculation to the default location
            pid = str(os.getpid())
            targetfile = psi4_scr + filename + '.default.' + pid + '.180.npy'
            shutil.copyfile(rundir + base_func + "/wfn.180.npy", targetfile)

            psi4.set_options({"basis": "def2-tzvp", "maxiter": 200 if "maxiter" not in psi4_config else psi4_config["maxiter"]})
            try:
                '''
                #Deprecated, see above
                if "b3lyp_hfx" in psi4_config and psi4_config["b3lyp_hfx"] != 20:
                    print("customized b3lyp with different HFX: ", psi4_config["b3lyp_hfx"])
                    e, wfn = psi4.energy("scf", self.get_hfx_functional('b3lyp', psi4_config["b3lyp_hfx"]),  molecule=mol, return_wfn=True)
                else:
                    e, wfn = psi4.energy('b3lyp', molecule=mol, return_wfn=True)
                wfn.to_file(rundir + "b3lyp/wfn.180")
                '''
                if "_hfx_" in base_func:
                    #if hfx in functional, expects in the form functional_hfx_XX, where XX is the percent of HFX included.
                    func = base_func.split('_')[0]
                    hfx = base_func.split('_')[-1]
                    e, wfn = psi4.energy("scf", dft_functional=self.get_hfx_functional(func, hfx),  molecule=mol, return_wfn=True)
                else:
                    e, wfn = psi4.energy(base_func, molecule=mol, return_wfn=True)
                wfn.to_file(rundir + base_func + "/wfn.180")
            except:
                print("This calculation does not converge.")
        #Check if the calculation succeeded
        success = run_utils.check_success(path=rundir + base_func)
        #remove copied wfn file, psi4 log files
        for filename in os.listdir('./'):
            if ("psi." in filename) or ("default" in filename):
                print("removing: :", filename)
                os.remove(filename)
        return success

    def run_general(self, functional="b3lyp", rundir="./", return_wfn=False,
                    psi4_scr='./', filename='output', verbose=False):
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
            verbose: bool
                If true, will print wfn the calculation is started from.
        """
        #Make the subdirectory, load relevant information
        psi4_config = self.config
        rundir = rundir + functional.replace("(", "l-").replace(")", "-r")
        with open(psi4_config["charge-spin-info"], "r") as f:
            d = json.load(f)
        psi4_config.update(d)
        run_utils = RunUtils()
        run_utils.ensure_dir(rundir)
        if verbose:
            print('Wfnfile:', psi4_config['wfnfile'])

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
        success = run_utils.check_success(path=rundir)
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
        success = run_utils.check_success(path=rundir)
        #remove temporary files
        for filename in os.listdir("./"):
            if ("psi." in filename) or ("default" in filename):
                print("removing: :", filename)
                os.remove(filename)
        return success

    def run_with_check(self, method, functional='b3lyp', rundir="./", return_wfn=False,
                       psi4_scr='./', filename='output', verbose=False, retry_scf=False):
        """
        Will run a psi4 calculation using the method specified in method, using the
        functional specified in functional.

        Note: if using run_initial, will ignore the functional argument, and instead use
        whatever is specified in psi4_config['base_functional'], or b3lyp if not specified.

        Parameters:
            method: str
                Either run_initial or run_general, depending on what calculation is being run.
            functional: str
                Name of the functional one wants to use in the calculation.
            rundir: str
                Path where functional calculations are run from.
                Results will be stored in the path rundir+functional.
            return_wfn: bool
                Whether or not the wfn file should be written after the calculation.
            psi4_scr: str
                Path of the Psi4 scratch directory.
            filename: str
                Name out the output .dat file containing results.
            verbose: bool
                If true, will print wfn the calculation is started from (for run_general).
            retry_scf: bool
                If true, will retry a calculation if it failed due to SCF iterations reached.
                Should

        Returns:
            success: bool
                Whether or not the calculation was successful.
        """

        psi4_config = self.config
        base_func = psi4_config['base_functional']

        if method == 'run_initial':
            print(f'Running initial calculation with base functional {base_func}.')
            functional = base_func
            #raise NotImplementedError('Method run_initial can only be used with the base functional specified in the config file.')

        if method == 'run_initial':
            runfunc = self.run_initial
        elif method == 'run_general':
            runfunc = self.run_general
        else:
            raise NotImplementedError("Please specify a valid method.")

        print(f"==={functional}===")
        #Clean name to prevent forbidden characters
        functional = functional.replace("(", "l-").replace(")", "-r")
        #If the corresponding folder not present, run the calculation with the functions above
        if not os.path.isdir(functional):
            if functional == base_func:
                success = runfunc(rundir=rundir, return_wfn=return_wfn,
                                  psi4_scr=psi4_scr, filename=filename)
            else:
                success = runfunc(functional=functional, rundir=rundir,
                                  return_wfn=return_wfn, psi4_scr=psi4_scr,
                                  filename=filename, verbose=verbose)
            print(f'Success: {success}')
        else:
            #If calculation already attempted, check for convergence
            print("folder exists.")
            resubed = False
            #resubmit if no output file or if iterations not reached
            if not os.path.isfile(functional + "/output.dat"):
                resubed = True
            else:
                with open(functional + "/output.dat", "r") as fo:
                    txt = "".join(fo.readlines())
                if "==> Iterations <==" not in txt or (not (("@DF-UKS iter" in txt) or ("@DF-RKS iter" in txt) or ("@DF-UHF iter" in txt) or ("@DF-RHF iter" in txt))):
                    resubed = True
            if resubed:
                print("previous errored out. resubmitting...")
                if functional == base_func:
                    success = runfunc(rundir=rundir, return_wfn=return_wfn,
                                      psi4_scr=psi4_scr, filename=filename)
                else:
                    success = runfunc(functional=functional, rundir=rundir,
                                      return_wfn=return_wfn, psi4_scr=psi4_scr,
                                      filename=filename, verbose=verbose)
                print("success: ", success)
            else:
                #If checks above pass, ensure no SCF error in calculation and that .wfn written
                with open(functional + "/output.dat", "r") as fo:
                    txt = "".join(fo.readlines())
                #Checking for the wfn only for if return_wfn=True, but False by default for run_general
                if 'PsiException: Could not converge SCF iterations' not in txt: # and os.path.isfile(functional + "/wfn.180.npy"):
                    success = True
                else: #Did not converge in the given number of SCF iterations
                    if retry_scf:
                        if functional == base_func:
                            success = runfunc(rundir=rundir, return_wfn=return_wfn,
                                              psi4_scr=psi4_scr, filename=filename)
                        else:
                            success = runfunc(functional=functional, rundir=rundir,
                                              return_wfn=return_wfn, psi4_scr=psi4_scr,
                                              filename=filename, verbose=verbose)
                    else:
                        success = False
                print("success: ", success)


        return success

    def rescue_with_check(self, functional='b3lyp', wfnpath='b3lyp/wfn.180.npy', alpha=20):
        """
        Checks if a calculation is converged. If it is not converged,
        runs a calculation using functional and alpha HFX, starting
        from the wave function in wfnpath.


        Parameters:
            functional: str
                Name of the functional one wants to converge.
            wfnpath: str
                Path to the wave function one wants the calculation to be run from.
            alpha: int
                HFX percentage to use.

        Returns:
            success: bool
                Whether or not a calculation succeeded.
        """
        #If the calculation at this HFX has not been tried before, run calculation
        if not os.path.isdir(functional + "-%d" % alpha):
            os.makedirs(functional + "-%d" % alpha)
            success = self.run_general_hfx(functional, hfx=alpha, wfn=wfnpath)
            print("success: ", success)
        else:
            #If the calculation (at jj HFX) has been tried before, check for convergence
            print("attempted rescue: ", functional, 'HFX', str(alpha))
            resubed = False
            #Need resubmission if no output or if iterations not reached in output
            if not os.path.isfile(functional + "-%d" % alpha + "/output.dat"):
                resubed = True
            else:
                with open(functional + "-%d" % alpha + "/output.dat", "r") as fo:
                    txt = "".join(fo.readlines())
                if "==> Iterations <==" not in txt:
                    resubed = True
            #Resubmit calculation using the specified wavefunction
            if resubed and os.path.isfile(wfnpath):
                print("previously errored out. resubmitting...")
                success = self.run_general_hfx(functional, hfx=alpha, wfn=wfnpath)
                print("success: ", success)
            #If the error was not due to SCF error and no wfn file written
            elif ('PsiException: Could not converge SCF iterations') not in txt and (not os.path.isfile(functional + "-%d" % alpha + "/wfn.180.npy")):
                #Try running the calculation again, move the old output to -timeout so it can be checked later if desired
                _functional = functional + "-%d" % alpha
                if not os.path.isfile(_functional + "/output-timeout.dat"):
                    shutil.copy(_functional + "/output.dat",
                                _functional + "/output-timeout.dat")
                    print("Time out, direct resub...")
                    success = self.run_general_hfx(functional, hfx=alpha, wfn=wfnpath)
                    print("success: ", success)
                else:
                    #Only try resubmission due to SCF errors once, if already tried, give up
                    print("Already submit once for timeout.")
                    success = False
            else:
                print("give up resubmission.")
                success = False

        return success
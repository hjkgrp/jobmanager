# Guide to using the Psi4 Workflow within Jobmanager

The Psi4 workflow is a utility that allows for one to do single-point DFT calculations using a number of different functionals in an automated way. First, a user converges a calculation using TeraChem, with B3LYP/LACVP* (or B3LYP/def2-SV(P)). Then, using the setup detailed below, the workflow does a one-step calculation in Psi4 to get a Psi4 wavefunction with B3LYP and the corresponding basis set. This is then used as the starting point to get a converged Psi4 wavefunction with B3LYP. Basis set projection can be carried out from LACVP* to def2-SV(P) to def2-TZVP, if desired. The B3LYP wavefunction is then used to initialize single-points using the other functionals, allowing DFT calculations on what should be the same electronic state on a number of different functionals. The single-points all have relatively small maximum number of iterations allowed, in order to ensure the structures have the same electronic state.

In order to run the Psi4 workflow, one needs to have their files set up in the following way:
 - A parent directory from which jobmanager will be run, and which contains the configuration files
 - A subdirectory for each structure, which each contain that structure's geometry (as a `.xyz` file), a `.molden` file of the structure, which was generated with TeraChem using B3LYP/LACVP*, and a `.json` file detailing some of the calculation parameters.

After all of these files are set up, one can simply run `jobmanager` in the parent directory (with the `jobmanager` conda environment active), which will then launch the Psi4 calculations specified in the configuration file. The specific inputs and outputs are detailed more below:

In the parent directory, one needs two files:
 - A file named `configure`, which contains the following: `run_psi4: psi4_config.json`, where `psi4_config.json` is the name of the `.json` file that one has specified the calculations they want to run.
 - A `.json` file that contains the following arguments (named `psi4_config.json` in this example folder):
   - `functional`: a list of functionals that you want to run on each structure. B3LYP will be run automatically, but can also be specified (although it has no effect). For functionals that support different Hartree-Fock exchange percentages, one can specify the degree of exchange by appending `_hfx_xx` to the end of the functional, where `xx` is the percentage of Hartree-Fock exchange one wants to use. Currently, only integer percentages of HFX are supported.
   - `xyzfile`: the name of the `.xyz` file in each subfolder. Each `.xyz` file must have this same name across all subdirectories that you want the Psi4 workflow to work on. This geometry should be the same as the one used for a TeraChem single point, or the final frame of a TeraChem geometry optimization (found in the `/scr/` directory of a TeraChem optimization).
   - `moldenfile`: the name of the `.molden` file in each subfolder. Each `.molden` file must have this same name across all subdirectories that you want the Psi4 workflow to work on. Note that the Psi4 workflow is expecting a molden from a TeraChem calculation using the LACVP* or def2-SV(P) basis set, which can be found in the `/scr/` directory of a TeraChem calculation.
   - `memory`: The amount of memory that Psi4 will use for the calculations.
   - `num_threads`: The number of threads that Psi4 will use for the calculations.
   - `basis`: The basis set to use for the Psi4 calculations. `lacvps` will use LACVP*, and `def2-tzvp` will first project the LACVP* calculation to def2-SV(P), converge that, and then use that as an initial guess for a def2-TZVP calculation.
   - `charge-spin-info`: The name of the `.json` file in each subfolder that contains the charge, spin, and spin treatment for each structure. Note that this file has to have the same name across all subfolders.
   - `hfx_rescue`: If true, allows for the calculation to be converged at different HFX and then that used as an initial guess for the desired calculation. Should be left as `false` in order to keep the initial guess the same for all functionals.
   - `wfnfile`: Leave as `b3lyp/wfn.180.npy`. Gives the name of the Psi4 wavefunction file that all other calculations are initialized from. This should be converged in Psi4 using the specified basis set and B3LYP, and is generated automatically when the workflow is run.
   - `bashrc`: The user's `.bashrc` file, used to help initialize the conda environment.
   - `conda_env`: The path to the `jobmanager` conda environment.

Note that the workflow has been updated so that calculations can be initialized from an arbitrary functional, instead of just B3LYP. In order to use this functionality, add the keyword `base_functional` to the `.json` file outlined above and set it to the name of the functional to be used to initialize other calculations. Then, edit the path in `wfnfile` in the `.json` file to be whatever the name of your base functional is plus `/wfn.180.npy`. The workflow will work exactly the same as above, but now will initialize calculations from the set functional instead of b3lyp.

Another functionality that has been added to this workflow is the HFX workflow. To use this, one adds the keyword `hfx_levels` to the `.json` file above, where the value is a list of integers that represent percentages of Hartree-Fock exchange to include. Then, for each functional specified in the `functional` list, they will be run at each of the HFX levels specified, starting from the nearest converged calculation.

By default, the workflow will only allow for 50 iterations on the calculation with the base functional (using the same functional and basis set as the provided TeraChem molden). You can increase this with the `init_maxiter` keyword in the `psi4_config.json` file. If the def2-TZVP basis is specified, up to 200 iterations are allowed to project from def2-SV(P) to def2-TZVP. This can be increased with the `maxiter` keyword in the `psi4_config.json` file. Similarly, all subsequent calculations that use the base functional's `.wfn` file default to a maximum of 50 iterations. This can be increased with the `maxiter` keyword.

Some additional options that can be added to one's calculations (via new arguments in the `psi4_config.json` file):
 - Instead of just the Psi4 wfn object, a molden file can be written for the result of each calculation. The user should specify the `write_molden_base` argument as true or false based on if a molden file for the base functional should be written. If a molden file for all of the non-base functionals is desired, the user should specify the `write_molden` argument as true/false. The default for both is false. By default, only the wfn file for the base functional is written, with no molden files. The default name for the molden files written is `geo_psi4.molden`, if a different name is desired, the user should specify the `write_molden_name` argument, where the corresponding entry is a string (including the `.molden` file extension).
 - If partial charges are desired, the user should specify the `partial_charge_base` or `partial_charge` arguments based on if they want partial charges for the base functional or for all other functionals, respectively. Both of these arguments expect a true or false value. To specify the charge scheme, the user should use the `partial_charge_scheme` argument, where the value can be either a string that is a valid argument of Psi4's `oeprop` function (https://psicode.org/psi4manual/master/oeprop.html), specifically one of `MULLIKEN_CHARGES`, `LOWDIN_CHARGES`, or `MBIS_CHARGES` for partial charge schemes; or they can submit a list of strings containing multiple such arguments if properties with several schemes are desired.
 - The default convergence criteria are 3e-5 a.u. for both energy and density. The user can change these with the `D_CONVERGENCE` or `E_CONVERGENCE` keywords, for density and energy, respectively.

In each subdirectory, one needs three files:
 - A `.xyz` file, matching the naming convention established in the parent directory's `.json` file. This should give the geometry of the structure corresponding to the wavefunction in the `.molden` file, so either the geometry used in a single-point or the final frame of a geometry optimization.
 - A `.molden` file, matching the naming convention established in the parent directory's `.json` file. This should be the converged molden for the structure, calculated using B3LYP/LACVP* (or B3LYP/def2-SV(P)), which is found in the `/scr/` directory of a TeraChem calculation, by default.
 - A `.json` file, matching the naming convention established in the parent directory's `.json` file. This should contain three arguments:
    - `charge`: The overall charge of the structure.
    - `spin`: The spin of the structure.
    - `ref`: Either `uks` if one wants to run unrestricted Kohn-Sham for spin specialization, or `rks` if one wants to run restricted Kohn-Sham.

After these files are generated, one runs `jobmanager` in the parent directory (the one with the `configure` file), and it should perform the specified calculations on the structures in each subfolder.

The output structure is as follows (in each subfolder):
 - `nohup.out` files tell what actions were taken by jobmanager and which calculations had to be rerun.
 - `nohup.err` files show any error messages that appeared during the calculations.
 - For each functional, a folder will be created, which contains the file `output.dat`, containing the output from the Psi4 calculation using that functional.
 - For B3LYP specifically, a Psi4 wavefunction file `wfn.180.npy` will be generated as well, which is used as an initial guess for the other functionals. The `1step` wavefunctions are generated in the process of converting the TeraChem molden to a Psi4 wavefunction, and can be mostly ignored.

The given example folder (`example_psi4_run`) should run if the `jobmanager` command is run from within that folder, and gives a sense for the desired file structure and components needed.

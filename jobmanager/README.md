# Configure file options
## Dependent jobs
Invoked if present
### Lists
- `solvent`
- `functionalsSP`
### Boolean
- `vertEA`
- `vertIP`
- `thermo`
- `dissociation`
- `hfx_resample`
- `mbe`
### Multiple lines
- `dissociated_ligand_charge`
- `dissociated_ligand_spinmult`

## jobmanager settings
- `sleep`
- `max_jobs`
- `max_resub`
- `hard_job_limit`
## checks and recovery
- `geo_check`: (default = False) supports oct, bidentate_oct, tridentate_oct, and tetrahedral checks. Only used if 'bad_geo' provided to `job_recovery`
- `ss_cutoff`: (default = 1.0) maximum deviation of S<sup>2</sup>. Only used if 'spin_contaminated' provided to `job_recovery`
- `job_recovery`: (default = []) modes of recovery
  - `scf`
    - if oscillating (standard deviation of energy > 1) resubmit with precision double, dftgrid 5, and dynamicgrid no
    - otherwise resubmit calc with `levela` = 1.0 and `levelb` = 0.1
  - `bad_geo`: resubmits bad geometries from initial geometry with frozen metal-ligand distances
  - `spin_contaminated`: resubmits spin contaminated calculations with BLYP
  - `thermo_grad_error`: resubmits calcs with 'Maximum component of gradient is too large' message using tighter convergence criteria
- `use_molscontrol`: (default = 0) integer determining whether or not to use molscontrol

## Calculation settings
- `levela`
- `levelb`
- `method`
- `hfx`
- `dispersion`

## Other run types
Provide json file path after :
- `general_sp`
- `run_psi4`

# Psuedocode Procedure
`resub.main`
Read `configure` and
- If not Psi4 run, call `resub.resub` until no jobs submitted or active for 3 cycles.
- Else, call `resub.resub_psi4`

`resub.resub`
1. Read `configure`
2. Check on all jobs with `moltools.check_completeness` 
3. Kill running jobs with SCF errors
4. Check for and submit dependent (derivative) jobs of finished jobs
5. Resubmit/recover jobs based on errors and `configure`

`moltools.check_completeness`
1. Run `tools.check_completeness`
2. Run geometry checks

`tools.check_completeness`
1. Find all files in subdirectories ending with `.out`
2. Filter out filenames starting with `.` or containing `nohup`, `_v[0-9].`, `slurm-`, `smd.out`, and `atom`
3. For each file
    i. Read in all lines
    ii. Check all lines for keywords `TeraChem` and `ORCA`
    iii. Check all lines for other keywords and output information
4. Get active jobs using bash
5. Categorize jobs by completion and errors
6. Get active jobs by finding all `_jobscript` files
7. Sort out active/finished conflicts
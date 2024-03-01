# Classes


## `Manager`
- Carry out the steps between generating inputs & getting results
  - Manage `Job` instances and interfaces with queue system
  - Recovery/resubmission attempts


## `Job` 
- Organize calculation(s) to run (independent of queue system)
  - Run history
  - Dependent calculation structure
  - Batch jobs


## `Calc` 
- Interface with software/write inputs
- Contains information about a single calculation
  - `type` (single point, minimization, TS search)


## `System`
- Interface with molSimplify
  - Generate/check geometries
- Contains information about a chemical system
  - `spin`
  - `charge`
  - `props` - Properties/possible variables (e.g. ligands, geometry, metal)


## `Project`
- Generate inputs from desired goal
  - Create and operate on `System` instances
- Compile/summarize results



# Classes


## `Manager`
- Spawns and manages `Run` instances

## `Run`
- Runs a `Job` instance using queue system or local machine

## `Job` 
- Contains information about the calculation(s) to run (independent of queue system)
  - quantum chemistry software to use
  - dependent calculations
  - 

## `Calc` 
- Contains information about a single calculation (independent of software)
  - `type` (single point, minimization, TS search)


## `System`
- Contains information about a chemical system & methods to generate via molSimplify
  - `spin`
  - `charge`
  - `props` - ligands, geometry, metal


## `Project`
- Create `Calc` and `Job` instances
- Compile/summarize results

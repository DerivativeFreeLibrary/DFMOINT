# DFMOINT - a derivative-free algorithm for multiobjective mixed-integer (constrained) problems

DFMOINT is a new linesearch-based solution method with theoretical convergence guarantees. In particular, it has been proved that every limit point of certain subsequences 
defined by the algorithm are stationary for the problem.
The continuous variables, are managed by employing a strategy which takes advantage of dense sequences of search directions.
The subset of variables that must assume discrete values is dealt with using primitive directions. 

## Authors

Giampaolo Liuzzi, Stefano Lucidi<br>
Department of Computer, Control and Management Engineering<br>
"Sapienza" University of Rome

#### Copyright

(2023) G. Liuzzi, S. Lucidi

DFMOINT is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

## Usage

By typing `python main.py -h` at command prompt you will get a brief usage description of the code:

```
usage: main [-h] [-v] [--constrained] [--print_filters] [--alg {DFMOINT,NSGA2,MACO,NOMAD}] [--tol TOL]
            [--max_fun MAX_FUN] [--seed SEED] [--outlev OUTLEV] [--export] [--probid PROBID] [--nvar {10,15,20,25,30}]
            [--nint NINT] [--ctype {a,b,c,d,e,f,z}] [--prbsel PRBSEL] [--runall] [--runall_dms]

Algorithm DFMOINT

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         Print version number
  --constrained         Solve constrained problems
  --print_filters       Print all filters while iterating
  --alg {DFMOINT,NSGA2,MACO,NOMAD}
                        Name of algo to be run
  --tol TOL             Tolerance in the stopping condition
  --max_fun MAX_FUN     Maximum number of function avaluations
  --seed SEED           seed for random seq. initialization
  --outlev OUTLEV       Output level
  --export              export results
  --probid PROBID       id of problem to be solved
  --nvar {10,15,20,25,30}
                        Number of variables
  --nint NINT           Number of integer variables (must be < nvar)
  --ctype {a,b,c,d,e,f,z}
                        type of constraints (z for unconstrained)
  --prbsel PRBSEL       prb of DMS collection
  --runall              run code on all CEC09 problems
  --runall_dms          run code on all DMS problems
```

If it is your intention to run DFMOINT on a user-defined problem, the easyest thing is to edit/modify suitably the python file `problem_factory.py` 

# bqpsolvers

Solver interfaces for bqpjson data files.
See the header of each file for dependencies.


### Solvers

The current solver connectors are,

* mip_gurobi.py - a MIP formulation of B-QP using gurobipy
* mip_cplex.py - a MIP formulation of B-QP using cplex
* lns_hfs.py - an LNS formulation using the HFS solver
* bop_ortools.py - a BOP formulation using or-tools
* aqc_dwave.py - a D-Wave based QUBO formulation using dwave_sapi2


### Input and Output

These scripts have a standardized output line, which begings with `BQP_DATA`.
There is an attempt to standardize the various solver's parameters via consistent command line arguments in these scripts.


### Tests

```
./mip_gurobi.py -f test/data/ran1_b_1.json 
./lns_hfs.py -f test/data/ran1_b_1.json 
./bop_ortools.py -f test/data/ran1_b_1.json 
```

### License
bqpsolvers is developed at Los Alamos National Laboratory and is provided under a BSD-ish license with a "modifications must be indicated" clause.  See the `LICENSE.md` file for the full text.  This package is part of the Hybrid Quantum-Classical Computing suite, known internally as LA-CC-16-032.
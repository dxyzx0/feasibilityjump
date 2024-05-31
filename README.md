# Feasibility Jump with Big Integer Support

This repository contains the reference C++ implementation of the Feasibility
Jump, which is a primal heuristic for mixed-integer linear programming.

See the [open access paper](https://link.springer.com/article/10.1007/s12532-023-00234-8) for details about
the algorithm.

## Usage

To build, run:

```
make
```

This should produce the executable `pbo_fj`, which solves MIP problems given
in `.mps` format. Its command line interface is:

```
Usage: pbo_fj [--timeout|t TIMEOUT] [--save-solutions|-s OUTDIR] [--verbose|-v] [--heuristic-only|-h] INFILE
```
 
For example, to use the FJ heuristic alone, and save solutions to the current
directory, run:

```
./pbo_fj -h benchmark/instances/miplib2017benchmark/academictimetablesmall.mps --save-solutions .
```

# T. Dufaud

This directory contains the code corresponding to the solution
of Poisson 1D problem by direct method or iterative method.
It is organized in three directories:
```
src/ 
include/
bin/
benchmarks/
tests/
```

"src" contains the source codes, "include" contains the 
header files and "bin" contains the executables. 
The "benchmarks" directory contains the benchmarking codes
and the "tests" directory contains the unit tests.
The compilation and execution can be done using the Makefile.

Here are the principal targets: 
testenv: bin/tp_testenv
tp2poisson1D_direct: bin/tpPoisson1D_direct
tp2poisson1D_iter: bin/tpPoisson1D_iter

The command,
```bash
$ make target # Compile an executable bin/target 
```

```bash
$ make all # compile the executable corresponding to all targets
```

```bash
$ make run_target #Execute ./bin/target
```

```bash
$ make clean
rm *.o bin/*
```

Here:
```bash
$ make run_testenv
$ make run_tpPoisson1D_iter
$ make run_tpPoisson1D_direct
```
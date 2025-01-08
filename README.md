### Notes de moi :
- La documentation est dans le dossier `docs/` faite avec Doxygen.
- Les tests sont dans le dossier `tests/` et les benchmarks dans `benchmarks/`.
- Les réponses aux questions sont dans le fichier `rapport/Questions.md`.
- Le rapport de projet est dans le dossier `rapport/Compte_Rendu.md`.
- Pour les questions sur le format CSR/CSC, j'ai ajouté un fichier `src/lib_poisson1D_csr_csc.c` qui contient les définitions des fonctions pour les formats CSR et CSC. Les adaptations pour les algorithmes se trouvent à la suite des fonctions de la version classique dans `src/lib_poisson1D.c`.
- J'ai dû renommé les fonctions d'avant avec un `_tridiag` pour les différencier des implémentations avec les formats CSR/CSC.
- Désormais le binaire des méthodes itératives prend 2 arguments : méthode (0 pour Richardson, 1 pour Jacobi, 2 pour Gauss-Seidel) et format (0 pour classique, 1 pour CSR).

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
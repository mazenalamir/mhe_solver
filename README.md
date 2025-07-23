# A python module for nonlinear moving horizon estimation

This module implements the moving hozion estimation algorithm as described in the author's book:

---

M. Alamir, **Nonlinear Control of Uncertain Systems:** Conventional and Learning-Based alternatives with python, Springer-Nature, 2025.

---

## Moving horizon estimation and the module main features

Moving horizon estimation is based on on-line nonlinear programming solution over moving windows. The principle is to reconciliate the gathered measurements with some known although potentially partially known dynamics of the system being *observed*. Compared to some existing modules, `mhe_solver` shows the following appealing features: 

1. The simultaneous estimation of state and parameters (extended estimation problem)
2. The possibility to choose the updating rate through the parameter $q$ (see below), by so doing the estimator propagate the system's dynamics until the next updating instants. By comparing the $q$-associated computation time to the one really needed during the simulation, it is possible to adapt $q$ so that the resulting estimator is real-time implementable.
3. The module proposes a method to simulate the resulting mhe behavior given a previously simulated *true* system.

## Contents of the repository 

The repository contains the following files: 

- `requirements.txt`: Standard list of needed packages to be used through `pip install -r requirements.txt`
- `mhe_solver.py`: the main module implementing the algorithm and the simualtion facilities. This file contains also the implementation of the estimator for the user-defined modified Lorentz oscillator provided in the file `lorentz.py`.
- `user_defined_reactor.py`: An example of user-defined files that contained the only required information from the user.
- `notebook_reactor.ipynb`: A jupyter notebook containing the use-case on the previous user-defined example.




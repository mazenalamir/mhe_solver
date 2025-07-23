# mhe_solver

## A python module for nonlinear moving horizon estimation

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

## Hints & important features 

- The module considers that the input is a part of the measurement outputs. This is a non standard viewpoint which seemed to me quite obvious as the control input is always considered to be measured. Moreover this enables to add the appropirate noise if necessay. This is the reason why one needs the map that **extract the input from the output** since the input per se is necessary to inject in the dynamics in order to simualte the system.
  
- Notice that the module is defined for controlled system. However it can be used for autonomous systems provided that some fictitious input is defined as it is shown in the use-case defined in `user_defined_reactor.py` and `notebook_reactor.ipynb` files. The same holds for the presence of constraints for which the user-defined map should be defined. The following excerpts highlights the used *tricks*:

```python
# from the user_defined_reactor.py file

def u_reactor(y):
    return 0

def c_reactor(x, u, p):
    c = -1 
    return c
```

Notice how: 

- The map `u_reactor` function defining how the input is extracted from the output is obviously dummy but sufficient to execture the estimation on non-controlled system.
- The constraints map `c_reactor` defines an always satisfied scalar constraint. 

## Example of a user-defined file 

```python
import numpy as np
from casadi import vertcat

# The problem's parameters

nx, n_par, nu, nc, ny = 3, 4, 2, 2, 4
p_min = np.array([4, 0, 15, 0.0])
p_max = np.array([6, 4, 30, 0.2])
x_min = np.array([-10000]*nx)
x_max = np.array([+10000]*nx)

# The system's ODEs
def ode_lorentz(x,u,p):

    xdot = vertcat(
        p[0]*(x[1]-x[0])+u[0],
        x[0]*(p[1]-x[2])-x[1],
        x[0]*x[1]-p[2]*x[2]+p[3]*u[1]
    )
    return xdot

# The output map 
def h_lorentz(x,u,p):

    y = vertcat(x[0], x[2], u)
    return y

# How to extract the input from the measurement 
def u_lorentz(y):
    return y[-2:]

# The instant-wise constraints map 
def c_lorentz(x, u, p):

    c = vertcat(
        p[3]-p[1],
        p[1]-p[0]
    )
    return c
```

## How to cite the algorithms 

```
@Inbook{Alamir2025,
title="Nonlinear Moving-Horizon Extended Observers",
bookTitle="Nonlinear Control of Uncertain Systems: Conventional and Learning-Based Alternatives with python",
year="2025",
publisher="Springer-Nature",
isbn="ISBN-13-978-3031932861"
}

```




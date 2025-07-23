import numpy as np
from casadi import vertcat

nx, n_par, nu, nc, ny = 2, 1, 1, 1, 1
p_min = np.array([0.14])
p_max = np.array([0.18])
x_min = np.array([-10000]*nx)
x_max = np.array([+10000]*nx)

def ode_reactor(x,u,p):

    xdot = vertcat(
       -2*p[0]*x[0]**2,
       p[0]*x[0]**2
    )
    return xdot

def h_reactor(x,u,p):

    y = x[0]+x[1]
    return y

def u_reactor(y):
    return 0

def c_reactor(x, u, p):
    c = -1 
    return c


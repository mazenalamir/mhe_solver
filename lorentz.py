import numpy as np
from casadi import vertcat

# Definition of the parameters of the problem
nx, n_par, nu, nc, ny = 3, 4, 2, 2, 4
p_min = np.array([4, 0, 15, 0.0])
p_max = np.array([6, 4, 30, 0.2])
x_min = np.array([-10000]*nx)
x_max = np.array([+10000]*nx)

# The dynamics of the system 
def ode_lorentz(x,u,p):

    xdot = vertcat(
        p[0]*(x[1]-x[0])+u[0],
        x[0]*(p[1]-x[2])-x[1],
        x[0]*x[1]-p[2]*x[2]+p[3]*u[1]
    )
    return xdot

# The measured output 
def h_lorentz(x,u,p):

    y = vertcat(x[0], x[2], u)
    return y

# How to extract u from the measurement vector
def u_lorentz(y):
    return y[-2:]

# The constraints function 
# constraints are to be enforced instant-wise 
def c_lorentz(x, u, p):

    c = vertcat(
        p[3]-p[1],
        p[1]-p[0]
    )
    return c


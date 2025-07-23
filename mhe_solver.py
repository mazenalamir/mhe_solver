"""
This module implement the moving horizon estimation algorithm 
described the book: 

Nonlinear control of uncertain systems: 
Conventional and learning-based alternatives with python 
by Mazen Alamir
"""
import numpy as np
import pandas as pd
from time import time
from casadi import vertcat, Function, MX, mtimes
from casadi import nlpsol, reshape
from casadi import jacobian, if_else, power
import copy

class Container:
    # Dummy class to use dot-based attributes.
    pass

class MHE:

    def __init__(self,
                 nx, n_par, nu, nc, ny,
                 x_min, x_max,
                 p_min, p_max,
                 ode,
                 meas_func,
                 constraints,
                 u_from_y, # map extracting u from y
                 dt,    # maximum dt for runge kutta
                 N_O,   # observation horizon
                 q,     # observation delay
                 Q, R, Qf, Pf, rho_c, # Penalty matrices
                 max_iter=500
                 ):

        self.nx = nx
        self.nu = nu
        self.np = n_par
        self.nc = nc
        self.ny = ny
        self.x_min = x_min
        self.x_max = x_max
        self.p_min = p_min
        self.p_max = p_max
        self.dt = dt
        self.N_O = N_O
        self.q = q
        self.Q = Q
        self.R = R
        self.Qf = Qf,
        self.Pf = Pf
        self.rho_c = rho_c
        self.nz = (N_O+1)*nx+n_par
        self.max_iter = max_iter
        self.u_from_y = u_from_y
        # Create the casadi maps
        x = MX.sym('x', nx)
        u = MX.sym('u', nu)
        p = MX.sym('p', n_par)
        xdot = ode(x, u, p)
        # r.h.s of the ode (continuous time)
        self.f = Function('f', [x, u, p], [xdot])
        # One-step-ahead map (discrete-time)
        k1 = self.f(x, u, p)
        k2 = self.f(x+0.5*k1*dt, u, p)
        k3 = self.f(x + 0.5 * k2 * dt, u, p)
        k4 = self.f(x + k3 * dt, u, p)
        x_next = x + dt/6*(k1+2*(k2+k3)+k4)
        self.F = Function('F', [x, u, p], [x_next])
        # Measurement map
        y = meas_func(x, u, p)
        self.h = Function('h', [x,u,p], [y])
        # Constraints map
        c = constraints(x, u, p)
        self.c = Function('c', [x,u,p], [c])
        # The cost functions
        x_hat = MX.sym('xhat', nx)
        p_hat = MX.sym('phat', n_par)
        lam = MX.sym('lam', 1)
        z = MX.sym('z', self.nz)
        Y = MX.sym('Y', N_O*ny)
        J_eq, J_meas  = 0, 0
        J_constr, J_previous = 0, 0
        g, lbg, ubg = [], [], [] # for Casadi-Ipopt formulation
        p = z[-n_par:]
        #------ loop over the observation horizon
        for i in range(N_O):
            xi = z[i*nx:(i+1)*nx]
            xip1 = z[(i+1)*nx:(i+2)*nx]
            yi = Y[i*ny:(i+1)*ny]
            ui = u_from_y(yi)
            e_x = xip1-self.F(xi, ui, p)
            e_y = yi-meas_func(xi, ui, p)
            ci = self.c(xi, ui, p)
            g += [ci[j] for j in range(self.nc)]
            lbg += [-1e12] * self.nc
            ubg += [0] * self.nc
            J_eq += mtimes(e_x.T, mtimes(Q, e_x))
            J_meas += mtimes(e_y.T, mtimes(R, e_y))
            for j in range(nc):
                J_constr += if_else(ci[j]<0, 0, ci[j]*ci[j])
        e_xf = xi - x_hat
        e_p = p - p_hat
        J_previous += mtimes(e_xf.T, mtimes(Qf, e_xf))
        J_previous += mtimes(e_p.T, mtimes(Pf, e_p))
        #-- Filter and Normalize wrt time
        J_eq *= lam / (N_O * nx * self.dt)
        J_meas *= lam / (N_O * ny)
        J_previous *= (1-lam)
        J_constr *= rho_c/(nc * N_O)
        J_ipopt = J_eq+J_meas+J_previous
        J_grad = J_ipopt + J_constr
        #--- Define the exported functions
        self.J_eq = Function('J_eq', [z, Y, x_hat, p_hat, lam], [J_eq])
        self.J_meas = Function('J_meas', [z, Y, x_hat, p_hat, lam], [J_meas])
        self.J_previous = Function('J_previous', [z, Y, x_hat, p_hat, lam], [J_previous])
        self.J_constr = Function('J_constr', [z, Y, x_hat, p_hat, lam], [J_constr])
        self.J_ipopt = Function('J_ipopt', [z, Y, x_hat, p_hat, lam], [J_ipopt])
        self.J_grad = Function('J_grad', [z, Y, x_hat, p_hat, lam], [J_grad])

        # Define the Casadi-Ipopt solver
        self.problem = {'f':J_ipopt,
                        'x':z,
                        'g':vertcat(*g),
                        'p':vertcat(Y, x_hat, p_hat, lam)}

        self.solver = nlpsol('solver', 'ipopt',
                             self.problem, {'ipopt': {'max_iter': self.max_iter}})

        lb_state = np.array(list(x_min) * (self.N_O + 1)).flatten()
        ub_state = np.array(list(x_max) * (self.N_O + 1)).flatten()
        lb_p = np.array(p_min).flatten()
        ub_p = np.array(p_max).flatten()
        self.lbz = np.array(list(lb_state) + list(lb_p))
        self.ubz = np.array(list(ub_state) + list(ub_p))

        self.compute_mhe_ipopt = lambda z0, Y, x_hat, p_hat, lam: \
                            self.solver(x0=z0, lbx=self.lbz, ubx=self.ubz,
                                        lbg=lbg,ubg=ubg,
                                        p=vertcat(Y, x_hat, p_hat, lam))

        # Compute the gradient for the fast gradient algorithm

        G = jacobian(J_grad, z)
        self.grad = Function('grad', [z, Y, x_hat, p_hat, lam], [G])

    def compute_mhe_fg(self, z, Y, x_hat, p_hat, lam, gam, c, Niter):

        def project(z):

            # projection on the hypercube of admissible value
            z_proj = copy.copy(np.array(z).flatten())
            for i in range(self.nz):
                if z_proj[i] > self.ubz[i]:
                    z_proj[i] = self.ubz[i]
                elif z_proj[i] < self.lbz[i]:
                    z_proj[i] = self.lbz[i]

            return z_proj

        zstar_previous = copy.copy(z)
        eta_previous = copy.copy(z)

        for _ in range(Niter):

            G = self.grad(zstar_previous, Y, x_hat, p_hat, lam).full().flatten()
            eta = zstar_previous - gam * G
            zstar = project(eta + c * (eta - eta_previous))
            eta_previous = copy.copy(eta)
            zstar_previous = copy.copy(zstar)

        return zstar

    def generate_random_U(self,
                          nf=10,
                          T=5,
                          t_sim=10):

        # Return the profiles of u as a list
        # of arrays of shape (nu,)
        # so that it can be used by simulate

        N_sim = int(t_sim/self.dt)
        t = np.linspace(0, t_sim, N_sim)
        alpha = [np.random.randn(nf) for j in range(self.nu)]
        phi = 2*np.pi*np.random.rand(nf)
        U = np.array([np.array([
            alpha[j][i]*np.sin(2*i*np.pi*t+phi[i])
                  for i in range(nf)
        ]).sum(axis=0)
         for j in range(self.nu)]).T
        U = list(U)
        return U

    def simulate(self, x0, U, p):
        # utilities that simulates the system for
        # a given sequence of exogenous inputs
        X = [list(x0)]
        Y = []
        for i in range(len(U)):
            u = U[i]
            Y.append(list(self.h(X[-1],u,p).full().flatten()))
            X.append(list(self.F(X[-1],u, p).full().flatten()))
        X = np.array(X).reshape(-1, self.nx)
        Y = np.array(Y).reshape(-1, self.ny)
        t = np.array([i*self.dt for i in range(X.shape[0])])
        return t, X, Y

    def simulate_mhe(self, t, X, Y, x_hat_0, p_hat_0, lam,
                     option='ipopt', gam=0.05, c=0.8, Niter=20):

        # t, X and Y are outputs of self.sim
        # they are matrices.

        #----------- some reusable functions

        def compute_zstar(k, Y, X_hat, P_hat):

            # Determine the instants delimiting the
            # observation window.
            i1 = k - self.N_O - self.q
            i2 = i1 + self.N_O
            # collect the corresponding measurement and
            # previous estimation
            Y_minus = Y[i1:i2, :].ravel()
            x_hat = X_hat[i2]
            p_hat = P_hat[i2]
            # warm start for the initial guess
            lex = np.array([ell.flatten() for ell
                            in X_hat[i1:i2+1]]).flatten()
            z = np.array(list(lex) + list(p_hat)).flatten()
            # solve the optimization problem
            if option == 'ipopt':
                t1 = time()
                R = self.compute_mhe_ipopt(z, Y_minus,
                                           x_hat, p_hat, lam)
                cpu_ = time()-t1
                z_star = R['x'].full().flatten()
            elif option == 'fg':
                t1 = time()
                z_star = self.compute_mhe_fg(z, Y_minus,
                                             x_hat, p_hat,
                                             lam, gam, c, Niter)
                cpu_ = time()-t1
            # extract x_star and p_star from the computed solution
            x_star = z_star[0:-self.np].reshape(-1, self.nx)[-1, :]
            p_star = z_star[-self.np:]

            return x_star, p_star, cpu_

        def update_x_hat(k, Y, x_star, p_star):

            x = x_star
            for i in range(self.q):
                u = self.u_from_y(Y[k - self.q + i, :])
                x = self.F(x, u, p_star)

            return x.full()
        #--------------------------------------------------
        # Fill up the measurement buffer for the first
        # optimization to be possible to perform

        X_hat, P_hat = [x_hat_0], [p_hat_0]
        cpu, t_cpu = [], []
        k = 0 # the real time of the system
        # fill the Y_minus for the first time
        for i in range(self.N_O+self.q):
            x = X_hat[-1]
            u = self.u_from_y(Y[k,:])
            p = P_hat[-1]
            X_hat.append(self.F(x, u, p).full().flatten())
            P_hat.append(p)
            k += 1

        # Run the optimizer for the first time
        # and update x_hat accordingly
        x_star, p_star, cpu_ = compute_zstar(k, Y, X_hat, P_hat)
        cpu.append(cpu_)
        t_cpu.append(k * self.dt)
        X_hat[k] = update_x_hat(k, Y, x_star, p_star)
        P_hat[k] = p_star

        # Ready to start the  main loop
        # until the end of the simulated data
        while k < len(X) - self.q:

            for i in range(self.q):
                x = X_hat[-1]
                u = self.u_from_y(Y[k,:])
                p = P_hat[-1]
                X_hat.append(self.F(x, u, p).full().flatten())
                P_hat.append(p)
                k += 1

            x_star, p_star, cpu_ = compute_zstar(k, Y, X_hat, P_hat)
            cpu.append(cpu_)
            t_cpu.append(k * self.dt)

            # Update the value at instant k
            X_hat[k] = update_x_hat(k, Y, x_star, p_star)
            P_hat[k] = p_star

        X_hat.append(X_hat[-1])
        X_hat = np.array([ell.flatten() for ell in X_hat])
        P_hat.append(P_hat[-1])
        P_hat = np.array([ell.flatten() for ell in P_hat])
        t_cpu, cpu = np.array(t_cpu), np.array(cpu)

        nDelay = len(X)-len(X_hat)

        if nDelay > 0:
            X_hat = np.vstack([np.array([x_hat_0 for _ in range(nDelay)]), X_hat])
            P_hat = np.vstack([np.array([p_hat_0 for _ in range(nDelay)]), P_hat])

        return X_hat, P_hat, t_cpu, cpu

if __name__ == "__main__":

    from lorentz import ode_lorentz, h_lorentz, nx, nu, n_par, nc, ny
    from lorentz import x_min, x_max, p_min, p_max,c_lorentz, u_lorentz

    Q = np.eye(nx)
    R = np.eye(ny)
    Qf = np.eye(nx)
    Pf = np.eye(n_par)

    mhe = MHE(nx=nx, nu=nu, n_par=n_par,nc=nc, ny=ny,
              x_min= x_min, x_max=x_max,
              p_min=p_min, p_max=p_max,
              ode=ode_lorentz, meas_func=h_lorentz,
              constraints=c_lorentz,
              u_from_y=u_lorentz,
              dt=0.01, N_O=10, q=3,
              Q=Q, R=R, Qf=Qf, Pf=Pf, rho_c=1e2)


    p = np.array([5, 3, 18, 0.1])
    lam = 0.9
    U = 0.1*np.ones(nu * (mhe.N_O))
    x0 = np.ones(nx)
    Usim = list(U.reshape(-1, mhe.nu))
    t, X, Y = mhe.simulate(x0, Usim, p)
    X = X.ravel()
    Y = Y.ravel()
    # initialize x_hat and p_hat to the correct values
    x_hat = X[-mhe.nx:]
    p_hat = 1.0 * p
    # construct the corresponding z vector
    z = np.zeros(mhe.nz)
    z[0:nx*(mhe.N_O+1)] = X
    z[-mhe.np:] = p
    # detune the measurement to check cost functions
    Y_test = 1.0*Y
    for i in range(mhe.N_O):
        Y_test[i*ny:i*ny+ny-nu] = 1.0*Y[i*ny:i*ny+ny-nu]

    J_eq = mhe.J_eq(z, Y_test, x_hat, p_hat, lam)
    J_meas = mhe.J_meas(z, Y_test, x_hat, p_hat, lam)
    J_previous = mhe.J_previous(z, Y_test, x_hat, p_hat, lam)
    J_constr = mhe.J_constr(z, Y_test, x_hat, p_hat, lam)
    J_iopt = mhe.J_ipopt(z, Y_test, x_hat, p_hat, lam)
    J_grad = mhe.J_grad(z, Y_test, x_hat, p_hat, lam)
    print(J_eq, J_meas, J_previous, J_constr, J_iopt, J_grad)

    Usim = list(U.reshape(-1,mhe.nu))
    t, X, Y = mhe.simulate(x0, Usim, p)

    U_rand = mhe.generate_random_U(nf=6, T=5, t_sim=10)
    print(U_rand[0].shape)

    print(X.shape)
    print(Y.shape)
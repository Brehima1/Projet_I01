import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
from scipy.sparse import spmatrix
import copy
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import mesh
from darcy_1D import getmatrhs1D



def time_loop(mesh, K, bcfunc, rhsfun, t0, tf, dt, theta, uinit, lineaire=False):

    t = t0
    ts = [t0]
    u = [uinit]
    uold = uinit
    while(t < tf):
        if (1.1*dt >= tf -t):
            # mieux t = t + min(dt, tf-t)µ
            dt = tf -t
        
        bcleft, bcright = bcfunc
        lval = theta*bcleft['value'](t+dt) + (1-theta)*bcleft['value'](t)
        rval = theta*bcright['value'](t+dt) + (1-theta)*bcright['value'](t)
        bcl = {'type': bcleft['type'], 'value': lval}
        bcr = {'type': bcright['type'], 'value': rval}
        bcs = [bcl, bcr]
        if lineaire:
            unew = do_step(mesh, K, bcs, rhsfun(t+dt), dt, theta, uold)
        else:
            unew, critere = picard(mesh, K, bcs, rhsfun(t+dt), dt, theta, uold)
            #print("done ", t, "    critere= " , critere)

        
        u.append(unew)
        t += dt
        ts.append(t)
        uold = unew

    return ts, u

def do_step(mesh, K, bclist, rhsfun, dt, theta, uold):

    A, f = getmatrhs1D(mesh, K, bclist, rhsfun)
    if theta == 0:
        unew = uold + (dt*f -dt*A.dot(uold))/mesh.mes
    else:
        b = 1/theta * (f + mesh.mes*uold/dt - (1-theta)*A.dot(uold))
        B = A.copy()
        D =A.diagonal()
        D += mesh.mes/(theta*dt)
        B.setdiag(D)
        unew = spsolve(B, b)
    return unew

def picard(mesh, Kfun, bclist, rhsfun, dt, theta, uold):
    converged = False
    tol =1e-5

    ukold = uold

    while not(converged):
        K = Kfun(ukold)
        A, f = getmatrhs1D(mesh, K, bclist, rhsfun)
        b = 1/theta * (f + mesh.mes*uold/dt - (1-theta)*A.dot(uold))
        B = A.copy()
        D =A.diagonal()
        D += mesh.mes/(theta*dt)
        B.setdiag(D)
        uknew = spsolve(B, b)
        
        critere = np.max (np.abs (uknew-ukold )) / np.max(np.abs(ukold))
        #print('Critere   ', critere)
        converged = critere < tol
        ukold = uknew

    return uknew,critere
if __name__ == '__main__':
    test_homogene()
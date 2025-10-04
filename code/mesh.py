from collections import namedtuple
import numpy as np

Mesh = namedtuple('Mesh', ['domain', 'sommets', 'centres', 'mes', 'dist'])

def make_uniform_mesh_1(L, N):
    h =L/(N-1)
    domain = (-L/2, L/2)
    sommets = np.linspace(-L/2, L/2, N+1)
    centres = np.linspace(-L/2+h/2, L/2-h/2, N)
    mes = h * np.ones(N)
    dist = h * np.ones(N+1)
    dist[0] = h/2
    dist[N] = h/2

    return Mesh(domain, sommets, centres, mes, dist)

def make_uniform_mesh(L, N):
    #Cas sur la demie droite
    h =L/(N-1)
    domain = (0, L)
    sommets = np.linspace(0, L, N+1)
    centres = np.linspace(h/2, L-h/2, N)
    mes = h * np.ones(N)
    dist = h * np.ones(N+1)
    dist[0] = h/2
    dist[N] = h/2

    return Mesh(domain, sommets, centres, mes, dist)

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt

import mesh


class SingularMatrixError(Exception):

    def __init__(self, message):
        self.message = message

def compute_kbar(mesh, K, moyenne):
    """Calcule la perméabilité moyenne aux interfaces

    mesh: maillage
    K: tableau (une valeur par intervalle) ou nombre (si constante)
    moyenne: type de moyenne souhaité
    """

    Kbar = np.zeros_like(mesh.sommets)
    if moyenne == 'constante':
        if K.shape[0] == mesh.sommets.shape[0] - 1:  # Si K est défini sur les cellules
        
            Kbar[:-1] = K  
            Kbar[-1] = K[-1] 
        else:
            raise ValueError(f"Dimensions incompatibles : K a {K.shape[0]} éléments, mesh.sommets a {mesh.sommets.shape[0]}")
    
    elif moyenne == 'arithmetique':
        Kbar[1:-1] = (K[:-1] + K[1:])/2
        Kbar[0] = K[0]
        Kbar[-1] = K[-1]
    elif moyenne == 'harmonique':
        Kbar[1:-1] = (mesh.dist[1:-1]*K[:-1]*K[1:]) / ((mesh.sommets[1:-1]-mesh.centres[:-1])*K[1:] + (mesh.centres[1:] - mesh.sommets[1:-1])*K[:-1])
        Kbar[0] = K[0]
        Kbar[-1] = K[-1]
    return Kbar


def getmatrhs1D(mesh, K, bclist, rhs, moyenne='arithmetique'):
    """Calcule la matrice et le second membre du système

    Arguments: cf fonction solve

    Remarques:
    - tx contient les transmissivités. Notez qu'il y a 1 valeur de plus que le nombre d'intervalles
    - A est stockée sous forme de matrice creuse (scipy.sparse), en format csr. La fonction spdiags est           adaptée au cas 1D (3 diagonales)
    - Les conditions aux limites sont traitées à part, à la fin de la fonction
    """
    
    Kbar = compute_kbar(mesh, K, moyenne)
    tx = Kbar / mesh.dist
    # Pour utiliser dia_array, et traiter les CL à part
    tx[0] = 0
    tx[-1] = 0
    n = len(mesh.centres)
    DiagVecs = np.array([-tx[1:], tx[:-1]+tx[1:], -tx[:-1]])
    offsets = np.array((-1, 0, 1))
    A = sp.dia_array((DiagVecs, offsets), (n, n)).tocsr()

    b = mesh.mes*rhs(mesh.centres)

# Boundary conditions
    bcleft, bcright = bclist
    if bcleft['type'] == 'Dirichlet': 
        A[0,0] += Kbar[0]/mesh.dist[0]
        b[0] += Kbar[0]/mesh.dist[0]*bcleft['value']
    elif bclist[0]['type'] == 'Neumann':
        b[0] += bcleft['value']
        
    if bcright['type'] == 'Dirichlet':
        A[-1,-1] += Kbar[-1]/mesh.dist[-1]
        b[-1] += Kbar[-1]/mesh.dist[-1]*bcright['value']
    elif bcright['type'] == 'Neumann':
        b[-1] += bcright['value']

    return A, b


def solve(mesh, K, bclist, rhs, moyenne= 'arithmetique'):
    """Résolution de l'équation de Darcy 1D par volumes finis
    Cette fonction appelle getmatrhs1D pour calculer la matrice et lke seocnd membre
    puis résout en appelant scipy.sparse.linalg.solve (matrices creuses)

    Le cas d'un problème avec conditons aux limites de Neumann aux deux extrémités déclenche une erreur

    Arguments:
    mesh:    maillage
    K:       perméabilité, tableau (une valeur par intervalle) ou nombre (si constante)
    bclist:  conditions aux limites ([bcleft, bcright], chaque élément est un dict 
             contenant un champ 'type' et un chanp 'value'. 
             'type' peut valoir 'Dirichlet' ou 'Neumann´)
    rhs:     second membne. Fonction rhs(x), x comme mesh.centres, ou tableau de la même taille
    moyenne: choix 'harmonique' (défaut), 'artihmetique' (ne pas utiliser, 
             sauf pour raison pédagogique, 'constante'si c'est le cas
    """

    A, b = getmatrhs1D(mesh, K, bclist, rhs, moyenne)
    bcleft, bcright = bclist
    if bcleft['type'] == 'Neumann' and bcright['type'] == 'Neumann':
        raise SingularMatrixError('Singular matrix')
    u = spsolve(A, b)
    return u

def compute_flux(mesh, K, u, moyenne):
    """ Calcule les flux aux interfaces

    Doit être appelée après solve
    Arguemnts: mesh, K, moyenne: comme solve
               u: pression aux mailles
    """
    
    Kbar = compute_kbar(mesh, K, moyenne)
    tx = Kbar / mesh.dist
    tx[0] = 0
    tx[-1] = 0
    flux = tx[1:-1]*np.diff(u)
    return flux
if __name__ == '__main__':
    test()

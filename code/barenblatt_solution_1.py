import numpy as np
from scipy.integrate import trapz
from scipy.special import gamma

def barenblatt(x, t, M, m, n=1, num_points=10000):
    """
    Calcule la solution de Barenblatt dans le cas R^+.

    Paramètres :
    x : Position spatiale (doit être un tableau numpy de valeurs positives)
    t : Temps (doit être strictement positif)
    M : Masse totale
    m, n : Paramètres du modèle
    num_points : Nombre de points pour l'intégration numérique

    Retourne :
    U(x, t) : Valeur de la solution en chaque point x
    """
    if t <= 0:
        raise ValueError("Le temps t doit être strictement positif.")
    if np.any(x < 0):
        raise ValueError("La position x doit être dans R^+ (x >= 0).")

    # Paramètres du modèle
    alpha = n / (n * (m - 1) + 2)
    k = (m - 1) * alpha / (2 * n)
    s = (m-1)/m
    xi_max = 1/np.sqrt(k)  # Limite où (1 - k * xi^2) >= 0
    xi_values = np.linspace(-xi_max, xi_max, num_points)  # Intégration sur R^+
    integrand = np.maximum(1 - k * (xi_values**2), 0) ** (1 / (m - 1))
    # Calcul de l'intégrale avec la méthode des trapèzes
    d = n * np.pi ** (n / 2) / gamma(n / 2 + 1) * trapz(integrand, xi_values)
    # Calcul de C (ajusté pour R^+)
    C = ((2*M)/d)**s
    xi = x / t ** (alpha / n)
    F_xi = np.maximum(C - k * (xi**2), 0)**(1 / (m - 1))
    return t**(-alpha)*F_xi
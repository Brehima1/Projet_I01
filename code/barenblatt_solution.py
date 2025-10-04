from scipy.integrate import quad
import numpy as np

def barenblatt(m, n, M, x, t):
    """
    Barenblatt solution for the porous medium equation.

    Parameters:
        m : float
            Nonlinearity parameter (m > 1).
        n : int
            Spatial dimension.
        M : float
            Total mass of the solution.
        x : array_like
            Spatial coordinates (1D array).
        t : float
            Time.

    Returns:
        U : ndarray
            Solution values at given x and t.
    """
    if t == 0:
        return np.zeros_like(x)  # Return zero if t is zero

    if m <= 1:
        raise ValueError("Parameter m must be greater than 1 for the Barenblatt solution.")
    if n <= 0:
        raise ValueError("Spatial dimension n must be greater than 0.")

    # Compute alpha and k
    alpha = n / (n * (m - 1) + 2)
    k = (m - 1) * alpha / (2 * n)

    # Compute the integral to find normalization constant C
    def integrand(y):
        return np.maximum(0, 1 - k * y**2)**(1 / (m - 1)) * y**(n - 1)
    
    try:
        d, _ = quad(integrand, 0, np.inf)
    except Exception as e:
        raise RuntimeError("Integral computation failed.") from e

    gamma = n / (2 * (m - 1) * alpha)
    C = (M / (2 * d))**(1 / gamma)

    # Compute xi, F(xi), and U(x, t; M)
    xi = np.array(x) / (t**(alpha / n))
    
    # Apply the maximum condition before exponentiation
    arg_F = np.maximum(0, C - k * xi**2)
    F = arg_F**(1 / (m - 1))  # Compute only on non-negative values

    U = t**(-alpha) * F

    return U

import numpy as np
from scipy.special import lambertw

def true_tau_return(G, tau=0.311, delta=0.300):
    arg = - (G * delta / tau) * np.exp(delta / tau)
    # W_0 is the principal branch. For complex values, it returns the complex root.
    W = lambertw(arg, k=0)
    lam = (W / delta) - (1.0 / tau)
    alpha = np.real(lam)
    return -1.0 / alpha

for G in [0.1, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0]:
    print(f"G={G}, tau_env={true_tau_return(G):.3f}")

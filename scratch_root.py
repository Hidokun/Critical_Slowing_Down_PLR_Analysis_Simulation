import scipy.optimize as opt
import numpy as np
tau=0.311
delta=0.300
func = lambda x: np.tan(x) + (tau/delta)*x
x_c = opt.fsolve(func, 2.0)[0]
f_c = x_c/(2*np.pi*delta)
G_c = -1/np.cos(x_c)
print(f'x_c={x_c:.4f}, f_c={f_c:.4f}, G_c={G_c:.4f}')

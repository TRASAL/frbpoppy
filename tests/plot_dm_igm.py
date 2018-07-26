"""Check accuracy of using 1200z for the intergalactic dispersion measure."""
from frbpoppy.log import pprint
from scipy.integrate import quad as integrate
import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as cst

H_0 = 69.6
W_m = 0.286
W_b = 0.0486
W_v = 0.714
W_l = 0.6911
z_max = 8.0

# Initialize parameters
step = 0.001
W_k = 1.0 - W_m - W_v  # Omega curvature

if W_k != 0.0:
    pprint('Careful - Your cosmological parameters do not sum to 1.0')


def dm_igm(x):
    """Comoving distance (Hogg et al, 1999)."""
    return (1+z)/math.sqrt((W_m*(1+x)**3 + W_l))


dms = []
zs = []

f = 3*cst.c*H_0*W_b/(8*cst.pi*cst.G*cst.m_p)*1.05e-42

for z in np.arange(0, z_max+step, step):
    d = f*integrate(dm_igm, 0, z)[0]
    dms.append(d)
    zs.append(z)

plt.plot(zs, dms, label='Actual Ioka')
plt.plot(zs, [1200*z for z in zs], label='1200 Approx.')
plt.ylabel(r'DM$_{\text{IGM}}$ \si{\parsec\per\cm\cubed}$')
plt.xlabel('Redshift')
plt.yscale("log", nonposy='clip')
plt.tight_layout(pad=0)
plt.legend()
plt.savefig('./plots/dm_igm')

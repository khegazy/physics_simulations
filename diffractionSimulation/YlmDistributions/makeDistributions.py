import numpy as np
from scipy.special import sph_harm

theta = np.linspace(0, 2*np.pi, 200)
phi   = np.linspace(0, np.pi, 200)

phi, theta = np.meshgrid(phi, theta)
phi = phi.astype(np.double)
theta = theta.astype(np.double)

for j in np.arange(11)*2:
  sph = np.real(sph_harm(0, j, theta, phi))
  sph = sph.astype(np.double)
  sph /= np.sum(sph)
  sph.tofile("YlmDist_L-{}_M-{}".format(j, 0))

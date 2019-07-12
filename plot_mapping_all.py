import matplotlib.pyplot as plt
from math import *
from mpl_toolkits.mplot3d import Axes3D
import pylab
import numpy as np
from matplotlib import cm
import matplotlib.colors as colors
from scipy.interpolate import UnivariateSpline
np.set_printoptions(threshold=np.nan)


lowerlimit = 0
upperlimit = 25

f = open("flux_data.txt", "r")


numtheta = int(f.readline().split()[0])
numphi = int(f.readline().split()[0])

arraytheta = []

arrayphi = []

muon_flux = [[] for i in range(numtheta)]

i = 0

while i < numtheta:

    temp = f.readline().split()

    if len(temp) < 1:
        break

    for j in range(0, numphi):
        muon_flux[i].append(float(temp[j]))

    muon_flux[i] = np.array(muon_flux[i])

    i = i + 1

muon_flux = np.array(muon_flux)

temp = f.readline().split()

for i in range(0, numtheta):
    arraytheta.append(float(temp[i]))

temp = f.readline().split()

for i in range(0, numphi):
    arrayphi.append(float(temp[i]))

arraytheta = np.array(arraytheta)
arrayphi = np.array(arrayphi)

arraytheta2 = arraytheta
arrayphi2 = arrayphi

arrayphi, arraytheta = np.meshgrid(arrayphi, arraytheta)


fig = plt.figure()
ax = fig.add_subplot(111)
#print(array_phi_theta_flux)
f = plt.imshow(muon_flux[lowerlimit:upperlimit], extent=(np.amin(arrayphi), np.amax(arrayphi), np.amin(arraytheta[lowerlimit:upperlimit]), np.amax(arraytheta[lowerlimit:upperlimit])), aspect = 'auto')
plt.xlabel('phi')
plt.ylabel('cos theta')
plt.colorbar(f)
ax.set_title("total flux, normal B, 3 depth, pos and neg")
plt.clim(0.9995, 1.0005)
#print(array_phi_theta_flux[lowerlimit:upperlimit])

plt.show()
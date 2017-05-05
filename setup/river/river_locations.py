'''
Figure out where to input river water.
'''

import matplotlib.pyplot as plt
import cmocean.cm as cmo
import netCDF4 as netCDF
import numpy as np


d = netCDF.Dataset('../bathy/grid-smoothed-Positive-r0.4.nc')

mask_rho = d['mask_rho'][1:-1, 1:-1]  # no ghost cells
# to outline the rho mask without ghost cells but
# with correct indexing so that pcolormesh and lines
# are consistent
xrhog = np.arange(0.5, mask_rho.shape[1]+1.5)
yrhog = np.arange(0.5, mask_rho.shape[0]+1.5)
X, Y = np.meshgrid(xrhog, yrhog)
# for using with contour
xrhogp1 = np.arange(1, mask_rho.shape[1]+1)
yrhogp1 = np.arange(1, mask_rho.shape[0]+1)
Xp1, Yp1 = np.meshgrid(xrhogp1, yrhogp1)

# plot overall bathymetry
plt.figure()
plt.pcolormesh(X, Y, d['h'][1:-1, 1:-1], cmap=cmo.deep)
plt.contour(Xp1, Yp1, d['h'][1:-1, 1:-1], cmap=cmo.deep)
plt.contour(Xp1, Yp1, mask_rho, [0], colors='k')
plt.axis('tight')
plt.title('bathymetry')
plt.savefig('bathymetry.png', bbox_inches='tight')

# plot overall mask_rho
plt.figure()
plt.pcolormesh(X, Y, mask_rho, cmap='spring')
plt.axis('tight')
plt.title('mask_rho')
plt.savefig('mask_rho.png', bbox_inches='tight')

# plot zoomed in mask_rho with rho labels
plt.figure(figsize=(12, 8))
plt.pcolormesh(X, Y, mask_rho, cmap='spring')
plt.plot(X, Y, '0.2', lw=0.1)
plt.plot(X.T, Y.T, '0.2', lw=0.1)
plt.xlim(153, 175)
plt.ylim(90, 102)
for j in range(91, 102):
    for i in range(154, 175):
        plt.text(i, j, str(i) + ',\n' + str(j), fontsize=10,
                 horizontalalignment='center',
                 verticalalignment='center')
plt.title('mask_rho labeled')
plt.savefig('mask_rho_zoomed.png', bbox_inches='tight')

# plot zoomed in mask_rho with u labels
plt.figure(figsize=(12, 8))
plt.pcolormesh(X, Y, mask_rho, cmap='spring')
plt.plot(X, Y, '0.2', lw=0.1)
plt.plot(X.T, Y.T, '0.2', lw=0.1)
plt.xlim(153, 175)
plt.ylim(90, 102)
for j in range(90, 102):
    for i in range(154, 175):
        plt.text(i-0.5, j, str(i) + ',\n' + str(j), fontsize=10,
                 horizontalalignment='center',
                 verticalalignment='center')
plt.title('mask_rho with mask_u labeled')
plt.savefig('mask_rho_zoomed_ulabels.png', bbox_inches='tight')

# plot zoomed in mask_rho with v labels
plt.figure(figsize=(12, 8))
plt.pcolormesh(X, Y, mask_rho, cmap='spring')
plt.plot(X, Y, '0.2', lw=0.1)
plt.plot(X.T, Y.T, '0.2', lw=0.1)
plt.xlim(153, 175)
plt.ylim(90, 102)
for j in range(90, 102):
    for i in range(154, 175):
        plt.text(i, j-0.5, str(i) + ',\n' + str(j), fontsize=10,
                 horizontalalignment='center',
                 verticalalignment='center')
plt.title('mask_rho with mask_v labeled')
plt.savefig('mask_rho_zoomed_vlabels.png', bbox_inches='tight')

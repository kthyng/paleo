'''
Check smoothing from grid-shifted.nc to grid-smoothed.nc
'''

import matplotlib.pyplot as plt
import netCDF4 as netCDF
import cmocean.cm as cmo
import bathy_smoother.bathy_tools as bt
import numpy as np

smoothings = ['Positive', 'Negative', 'PlusMinus', 'Laplacian', 'lp']
rx0max = 0.4

gorig = netCDF.Dataset('grid-original.nc')
gshift = netCDF.Dataset('grid-shifted.nc')
gsm = [netCDF.Dataset('grid-smoothed-' + smoothing + '-r' + str(rx0max) + '.nc') for smoothing in smoothings]

# plt.pcolormesh(gshift['h'][:] - gsmooth['h'][:], cmap=cmo.cdom)
# plt.colorbar()

# plot roughness matrices
rm = []
for i, smoothing in enumerate(smoothings):
    # rm.append(bt.RoughnessMatrix(gsm[i]['h'][:], np.ones_like(gsm[i]['h'][:])))
    rm.append(bt.RoughnessMatrix(gsm[i]['h'][:], gsm[i]['mask_rho'][:]))

fig, axes = plt.subplots(1, len(smoothings)+1, figsize=(18, 4))
mappable = axes[0].pcolormesh(gshift['h'][:]*gshift['mask_rho'][:], cmap=cmo.bathy)
axes[0].contour(gshift['mask_rho'][:], [0], colors='k')
axes[0].set_title('shifted bathy')
plt.colorbar(mappable, ax=axes[0])

for i in range(len(smoothings)):
    hdiff = gshift['h'][:] - gsm[i]['h'][:]
    mappable = axes[i+1].pcolormesh(rm[i], cmap=cmo.speed)
    axes[i+1].contour(gsm[i]['mask_rho'][:], [0], colors='k')
    axes[i+1].set_title(smoothings[i] + '\n' + 
                        str(abs(rm[i]).max()))
    plt.colorbar(mappable, ax=axes[i+1])
fig.tight_layout()
fig.savefig('roughness_comp.png', bbox_inches='tight')

# plot grid files
fig, axes = plt.subplots(1, len(smoothings)+1, figsize=(18, 4))
mappable = axes[0].pcolormesh(gshift['h'][:]*gshift['mask_rho'][:], cmap=cmo.bathy)
axes[0].contour(gshift['mask_rho'][:], [0], colors='k')
axes[0].set_title('shifted bathy')
plt.colorbar(mappable, ax=axes[0])

for i in range(len(smoothings)):
    hdiff = gshift['h'][:] - gsm[i]['h'][:]
    mappable = axes[i+1].pcolormesh(hdiff[:]*gsm[i]['mask_rho'][:],
                                    cmap=cmo.vel, vmin=-50, vmax=50)
    axes[i+1].contour(gsm[i]['mask_rho'][:], [0], colors='k')
    axes[i+1].set_title(smoothings[i] + '-shifted bathy\n' + 
                        str(abs(hdiff).sum()))
    plt.colorbar(mappable, ax=axes[i+1])
fig.tight_layout()
fig.savefig('smoothing_comp.png', bbox_inches='tight')

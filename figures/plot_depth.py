'''
Plot depth near mouth
'''

import matplotlib.pyplot as plt
import cmocean.cm as cmo
import netCDF4 as netCDF
import tracpy
import numpy as np
import os


proj = tracpy.tools.make_proj('nwgom', usebasemap=True, **{'llcrnrlon': -98,
                              'llcrnrlat': 26})
grid = tracpy.inout.readgrid('../grid.nc', proj)

var = grid.h

# Find largest salinity value in water column
vmax = 100
vmin = 50

fig = plt.figure(figsize=(10, 6))

fname = 'depth.png'

ax = fig.add_subplot(111)
proj.drawcoastlines(ax=ax)
proj.fillcontinents('0.8', ax=ax)
proj.drawparallels(np.arange(18, 35), dashes=(1, 1), linewidth=0.15,
                   labels=[1, 0, 0, 0], ax=ax)
proj.drawmeridians(np.arange(-100, -80), dashes=(1, 1), linewidth=0.15,
                   labels=[0, 0, 0, 1], ax=ax)
ax.contour(grid.x_rho, grid.y_rho, grid.mask_rho, [0],
           colors='0.3', linewidths=0.5, alpha=0.7)
mapp = ax.pcolormesh(grid.x_psi, grid.y_psi, var[1:-1, 1:-1], cmap=cmo.deep,
                     vmin=vmin, vmax=vmax)
fig.colorbar(mapp, shrink=0.77, ax=ax)
fig.tight_layout()
fig.savefig(fname, bbox_inches='tight')

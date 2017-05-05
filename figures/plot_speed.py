'''
Plot speed near mouth
'''

import matplotlib.pyplot as plt
import cmocean.cm as cmo
import netCDF4 as netCDF
import tracpy
import numpy as np
import os


restart = True  # use to show last result from restart file for blowup case

proj = tracpy.tools.make_proj('nwgom', usebasemap=True, **{'llcrnrlon': -98,
                              'llcrnrlat': 18, 'urcrnrlon': -78, 'urcrnrlat': 30})
grid = tracpy.inout.readgrid('../grid.nc', proj)

if restart:
    d = netCDF.Dataset('../out/ocean_rst_old25.nc')
    # d = netCDF.Dataset('../ocean_rst.nc')
else:
  # d = netCDF.Dataset('../out/ocean_his_0001_old21.nc')
    d = netCDF.Dataset('../ocean_his_0001.nc')
u = d['u']
v = d['v']

vmax = 1.5
vmin = 0
cmap = cmo.speed

t = d['ocean_time']
dates = netCDF.num2date(t[:], t.units)

if not os.path.exists('speedmax'):
    os.mkdir('speedmax')

fig = plt.figure(figsize=(10, 6))

for it in range(u.shape[0]):
    fname = 'speedmax/' + str(it).zfill(2) + '.png'
    if os.path.exists(fname):
        continue

    # on rho grid
    var = np.sqrt(tracpy.op.resize(u[it, :, 1:-1, :], 2)**2 +
                  tracpy.op.resize(v[it, :, :, 1:-1], 1)**2)
    # Find largest value in water column
    varmax = var.max(axis=0)

    ax = fig.add_subplot(111)
    proj.drawcoastlines(ax=ax)
    proj.fillcontinents('0.8', ax=ax)
    proj.drawparallels(np.arange(18, 35), dashes=(1, 1), linewidth=0.15,
                       labels=[1, 0, 0, 0], ax=ax)
    proj.drawmeridians(np.arange(-100, -78, 2), dashes=(1, 1), linewidth=0.15,
                       labels=[0, 0, 0, 1], ax=ax)
    ax.contour(grid.x_rho, grid.y_rho, grid.h, np.arange(500, 4000, 500),
               colors='0.3', linewidths=0.5, alpha=0.7)
    mapp = ax.pcolormesh(grid.x_psi, grid.y_psi, varmax, cmap=cmap,
                         vmin=vmin, vmax=vmax)
    ax.set_title(dates[it].isoformat())
    fig.colorbar(mapp, shrink=0.77, ax=ax, extend='max')
    fig.tight_layout()
    fig.savefig(fname, bbox_inches='tight')
    fig.clear()

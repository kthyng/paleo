'''
Plot salinity near mouth
'''

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import netCDF4 as netCDF
import tracpy
import numpy as np
import os


restart = False  # use to show last result from restart file for blowup case

proj = tracpy.tools.make_proj('nwgom', usebasemap=True, **{'llcrnrlon': -98,
                              'llcrnrlat': 18, 'urcrnrlon': -78, 'urcrnrlat': 30})
# proj = tracpy.tools.make_proj('nwgom', usebasemap=True, **{'llcrnrlon': -98,
#                               'llcrnrlat': 26})
grid = tracpy.inout.readgrid('../grid.nc', proj)

if restart:
    loc = '../runs/new1/ocean_rst.nc'
else:
    loc = '../runs/new3/ocean_his_0001.nc'
d = netCDF.Dataset(loc)
salt = d['salt']

if ('new1' in loc) or ('new2' in loc):
    # Find largest salinity value in water column
    saltlow = salt[:, :, 1:-1, 1:-1].min(axis=1)
    cmap = cmo.haline
else:
    # Find largest salinity value in water column
    saltlow = salt[:, :, 1:-1, 1:-1].max(axis=1)
    cmap = cmo.haline_r

if 'new1' in loc:
    smax = 36.5; smin = 33.9
elif 'new2' in loc:
    smax = 36.225; smin = 35.2
elif 'new3' in loc:
    smax = 36.43; smin = 36.225
elif 'new4' in loc:
    smax = 37.6; smin = 36.225

t = d['ocean_time']
dates = netCDF.num2date(t[:], t.units)

fig = plt.figure(figsize=(10, 6))

for it in range(salt.shape[0]):
    fname = 'saltmax/' + str(it).zfill(3) + '.png'
    # if os.path.exists(fname):
    #     continue
    ax = fig.add_subplot(111)
    proj.drawcoastlines(ax=ax)
    proj.fillcontinents('0.8', ax=ax)
    proj.drawparallels(np.arange(18, 35), dashes=(1, 1), linewidth=0.15,
                       labels=[1, 0, 0, 0], ax=ax)
    proj.drawmeridians(np.arange(-100, -80, 2), dashes=(1, 1), linewidth=0.15,
                       labels=[0, 0, 0, 1], ax=ax)
    ax.contour(grid.x_rho, grid.y_rho, grid.h, np.arange(500, 4000, 500),
               colors='0.3', linewidths=0.5, alpha=0.7)
    mapp = ax.pcolormesh(grid.x_psi, grid.y_psi, saltlow[it], cmap=cmap,
                         vmin=smin, vmax=smax)
    ax.set_title(dates[it].isoformat())
    fig.colorbar(mapp, shrink=0.77, ax=ax)
    fig.tight_layout()
    plt.show()
    fig.savefig(fname, bbox_inches='tight')
    fig.clear()

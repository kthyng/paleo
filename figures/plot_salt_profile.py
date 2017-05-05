'''
Plot salinity profile near mouth
'''

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import cmocean.cm as cmo
import netCDF4 as netCDF
import tracpy
import numpy as np
import os
import octant.depths


proj = tracpy.tools.make_proj('nwgom', usebasemap=True, **{'llcrnrlon': -98,
                              'llcrnrlat': 18, 'urcrnrlon': -78, 'urcrnrlat': 30})
grid = tracpy.inout.readgrid('../grid.nc', proj)

loc = '../runs/new3/ocean_his_0001.nc'
d = netCDF.Dataset(loc)
salt = d['salt']

if ('new1' in loc) or ('new2' in loc):
    cmap = cmo.haline
else:
    cmap = cmo.haline_r

if 'new1' in loc:
    smax = 36.5; smin = 33.9
elif 'new2' in loc:
    smax = 36.225; smin = 35.2
elif 'new3' in loc:
    smax = 36.5; smin = 36.225
elif 'new4' in loc:
    smax = 37.6; smin = 36.225

grid.Vtransform = 2
grid.Vstretching = 4
grid.km = 20
grid.theta_s = 0.0001
grid.theta_b = 0.0
grid.hc = 0

grid.zrt0 = octant.depths.get_zrho(grid.Vtransform, grid.Vstretching,
                                   grid.km, grid.theta_s,
                                   grid.theta_b, grid.h, grid.hc,
                                   zeta=0, Hscale=3)

# Find largest salinity value in water column
saltlow = salt[:, :, 1:-1, 161]

# xlimits
xind = [28,100]

t = d['ocean_time']
dates = netCDF.num2date(t[:], t.units)

# flip left and right so river mouth at left of plot
xvar = grid.y_rho[1:-1,161]; xvar = xvar[::-1]
saltlow = saltlow[:,:,::-1]
zvar = grid.zrt0[:,1:-1,161]; zvar = zvar[:,::-1]

dist = abs(xvar - xvar[xind[0]])/1000.

fig = plt.figure(figsize=(10, 5))

for it in range(salt.shape[0]):
    # it = 100
    fname = 'saltprof/' + str(it).zfill(3) + '.png'
    # if os.path.exists(fname):
    #     continue
    ax = fig.add_subplot(111)
    mapp = ax.pcolormesh(dist, zvar, saltlow[it], cmap=cmap,
                         vmin=smin, vmax=smax)
    ax.invert_yaxis()
    ax.set_title(dates[it].isoformat())
    ax.set_xlim(dist[xind[0]], dist[xind[1]])
    ax.set_xlabel('distance along axis [km]')
    ax.set_ylabel('depth [m]')
    fig.colorbar(mapp, ax=ax)#, shrink=0.77)
    fig.tight_layout()
    plt.show()
    fig.savefig(fname, bbox_inches='tight')
    fig.clear()

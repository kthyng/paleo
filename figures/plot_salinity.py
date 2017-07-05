'''
Plot salinity near mouth
'''

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import cmocean.cm as cmo
# import netCDF4 as netCDF
# import tracpy
import numpy as np
import os
import cartopy.crs as ccrs
import xarray as xr
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker
import pandas as pd


restart = False  # use to show last result from restart file for blowup case
which = 'new3'
density = '1028'

# proj = tracpy.tools.make_proj('nwgom-pyproj', **{'llcrnrlon': -98,
#                               'llcrnrlat': 18, 'urcrnrlon': -78, 'urcrnrlat': 30})
# proj = tracpy.tools.make_proj('nwgom', usebasemap=True, **{'llcrnrlon': -98,
#                               'llcrnrlat': 18, 'urcrnrlon': -78, 'urcrnrlat': 30})
# proj = tracpy.tools.make_proj('nwgom', usebasemap=True, **{'llcrnrlon': -98,
#                               'llcrnrlat': 26})
# grid = tracpy.inout.readgrid('../grid.nc', proj)
# grid = octant.grid.CGrid_geo(lonverts, latverts, proj)


if restart:
    loc = '../runs/' + which + '/' + density + '/ocean_rst.nc'
else:
    loc = '../runs/' + which + '/' + density + '/ocean_his_0001.nc'
# d = netCDF.Dataset(loc)
# salt = d['salt']
d = xr.open_dataset(loc)
salt = d.salt#.data

if ('new1' in loc) or ('new2' in loc):
    # Find largest salinity value in water column
    saltlow = salt[:, :, 1:-1, 1:-1].min(axis=1)
    cmap = cmo.haline
else:
    # Find largest salinity value in water column
    saltlow = salt.isel(eta_rho=slice(1,-1), xi_rho=slice(1,-1)).max(axis=1)
    # saltlow = salt[:, :, 1:-1, 1:-1].max(axis=1)
    cmap = cmo.haline_r

if 'new1' in loc:
    smax = 36.5; smin = 33.9
elif 'new2' in loc:
    smax = 36.225; smin = 35.2
elif 'new3' in loc:
    # smax = 40; smin = 36.225
    smax = 36.35; smin = 36.225  # old
elif 'new4' in loc:
    smax = 37.6; smin = 36.225
elif 'new5' in loc or 'new6' in loc:
    # smax = 37.5; smin = 36.225
    smax = 40; smin = 36.225

dates = d.ocean_time
# t = d['ocean_time']
# dates = netCDF.num2date(t[:], t.units)

if not os.path.exists('saltmax/'):
    os.mkdir('saltmax/')
if not os.path.exists('saltmax/' + which):
    os.mkdir('saltmax/' + which)
if not os.path.exists('saltmax/' + which + '/' + density):
    os.mkdir('saltmax/' + which + '/' + density)

fig = plt.figure(figsize=(10, 6))
        # fig = plt.figure(figsize=(9.4, 7.7), dpi=100)

# for it in range(salt.shape[0]):
for date in dates:
    fname = 'saltmax/' + which + '/' + density + '/' + pd.to_datetime(date.data).isoformat()[0:13] + '.png'
    # fname = 'saltmax/' + str(it).zfill(3) + '.png'
    # if os.path.exists(fname):
    #     continue
    ax = fig.add_axes([0.06, 0.01, 0.93, 0.95], projection=ccrs.Mercator(central_longitude=-85.0))
    ax.set_frame_on(False) # kind of like it without the box
    ax.set_extent([-98, -87.5, 22.8, 30.5], ccrs.PlateCarree())
    gl = ax.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
    # the following two make the labels look like lat/lon format
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabels_top = False  # turn off labels where you don't want them
    gl.ylabels_right = False

    # plot isobaths
    ax.contour(d.lon_rho, d.lat_rho, d.h, np.arange(500, 4000, 500), colors='0.6', transform=ccrs.PlateCarree(), linewidths=0.5)
    # ax = fig.add_subplot(111)
    # proj.drawcoastlines(ax=ax)
    # proj.fillcontinents('0.8', ax=ax)
    # proj.drawparallels(np.arange(18, 35), dashes=(1, 1), linewidth=0.15,
    #                    labels=[1, 0, 0, 0], ax=ax)
    # proj.drawmeridians(np.arange(-100, -80, 2), dashes=(1, 1), linewidth=0.15,
    #                    labels=[0, 0, 0, 1], ax=ax)
    # ax.contour(grid.x_rho, grid.y_rho, grid.h, np.arange(500, 4000, 500),
    #            colors='0.3', linewidths=0.5, alpha=0.7)
    # import pdb; pdb.set_trace()
    mapp = ax.pcolormesh(d.lon_psi, d.lat_psi, saltlow.sel(ocean_time=date), cmap=cmap,
                         vmin=smin, vmax=smax, transform=ccrs.PlateCarree())
    datestr = pd.to_datetime(date.data).strftime('%Y %b %02d %H:%M')
    ax.set_title(datestr)
    fig.colorbar(mapp, shrink=0.77, ax=ax)
    # fig.tight_layout()
    plt.show()
    fig.savefig(fname, bbox_inches='tight')
    fig.clear()

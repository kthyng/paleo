#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py
Created on 2008-04-24

Copyright (c) 2008 Robert D. Hetland
"""
__docformat__ = "restructuredtext en"

from datetime import datetime, timedelta

import numpy as np
import netCDF4 as netCDF
import os


which = 'new6'
river_density = 1029

time_units = 'days since 1970-01-01 00:00:00'
startdate = datetime(2010, 1, 1, 0, 0)
nmonths = 12
runtime = 30*nmonths*86400.  # about nmonths months in seconds
driver = 3*3600.  # 3 hour steps in seconds
# every 3 hours
dates = [startdate + timedelta(seconds=driver)*i for i in range(int(runtime/driver))]
river_jd = netCDF.date2num(dates, units=time_units)

# first east-west, then north-south
river_Xposition = \
    np.array([160., 160., 160., 159., 160., 161.,
              155., 156., 157., 158., 159., 162., 163., 164., 165., 166., 167., 168., 169., 170., 171., 172., 173., 174.])

river_Eposition = \
    np.array([92., 93., 94., 95., 96., 97.,
              92., 92., 92., 92., 92., 99., 99., 99., 98., 98., 98., 98., 97., 97., 97., 97., 97., 97.])

river_direction = \
    np.array([0., 0., 0., 0., 0., 0.,
              1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.])

ntimes = len(river_jd)  # number of river times
nriver = len(river_Xposition)  # number of river inputs
river = np.arange(1, nriver+1)
s_rho = 20

# calculate river properties: depth_river, temp_river, salt_river
os.chdir('../initialization')
import make_initialization
_, river_temp, river_salt = make_initialization.calcs(river_density, doprint=False)
# import runpy
# file_globals = runpy.run_module("make_initialization.py 0 0 --calc_river " + str(river_density))
# os.system('/opt/anaconda2/bin/python2 make_initialization.py 0 0 --calc_river ' + str(river_density))
os.chdir('../river')

# Set river temp
# river_temp = 18.769

# salinity to higher than background salinity of 36.225
# river_salt = 38.909

# vertical profile
# Want river input to be weighted by the thickness of each layer, which is being
# calculated as (due to Vtransform=2 and Vstretching=4):
theta_b=0; theta_s=0.0001
s = np.linspace(-1.0, 0.0, s_rho+1)
C = (1.0 - np.cosh(theta_s * s)) / (np.cosh(theta_s) - 1.0)

# vertical profile for river water input in
# Weight according to layer thickness to be Uniform
# num vertical levels x number of rivers
river_Vshape = (np.diff(C))[:, np.newaxis].repeat(nriver, axis=1)
# # do linear. bottom to surface ordering probably.
# river_Vshape = np.linspace(0, 2.0/s_rho, s_rho)[:, np.newaxis].repeat(nriver, axis=1)

# number for vertically-integrated transport into each cell, m^3/s
# # Changing this to be total flux over time, so units are m^3, still divided into # river inputs
# transport = 1e14/nriver  # see evernote for calculation
# Changing this to be total flux over time, so units are m^3, NOT divided into # river inputs
if which == 'new5':
    transport = 3e14
elif which == 'new6':
    transport = 1.2e15

# sign for transport — gives direction of input into rho cells
direction = np.array([1., 1., 1., 1., 1., 1.,
                      -1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1., -1.])

# Change in river transport over time
t = np.linspace(0, runtime, ntimes)
mid = 0.5*runtime  # half way
width = 0.15*runtime
# can solve for volume amplitude
# amp should be in m^3/s
amp = transport/((np.exp(-(t - mid)**2/(2*width**2))).sum()*driver)  # peak volume
gauss = amp*np.exp(-(t - mid)**2/(2*width**2))
# plt.plot(gauss/1e6); plt.show()

# time vs number of rivers
river_transport = np.ones((ntimes, nriver))*direction
# transport changes in time as a gaussian
river_transport *= gauss[:,np.newaxis]

# now weight across rivers according to cell wall width
xfac = 5000.  # m: average cell wall width in this direction
yfac = 7000.  # m: average cell wall width in this direction
iy = river_direction == 0  # cell wall river enters is delta y
ix = river_direction == 1  # cell wall river enters is delta x
walltot = iy.sum()*yfac + ix.sum()*xfac  # total of cell wall width weights, to divide by
# size is number of river inputs
# weights = iy.astype(float)*(1-xfac/yfac+1) + ix.astype(float)*(xfac/yfac)
weights = iy.astype(float)*yfac/walltot + ix.astype(float)*xfac/walltot  # unit norm
# weight across river inputs by cell wall face
river_transport *= weights

fname = which + '_' + str(river_density) + '.nc'
nc = netCDF.Dataset(fname, 'w', format='NETCDF4')
# nc = netCDF.Dataset('tglo_river_frc_2010-01-01.nc', 'w', format='NETCDF4')
nc.createDimension('x_rho', 256)
nc.createDimension('eta_rho', 128)
nc.createDimension('s_rho', s_rho)
nc.createDimension('river', nriver)
nc.createDimension('river_time', ntimes)


def write_nc_var(name, dimensions, var, units=None):
    nc.createVariable(name, 'f8', dimensions)
    if units is not None:
        nc.variables[name].units = units
    nc.variables[name][:] = var

write_nc_var('river', ('river', ), river, 'River index')
write_nc_var('river_Xposition', ('river', ), river_Xposition, 'River xi index')
write_nc_var('river_Eposition', ('river', ), river_Eposition, 'River eta index')
write_nc_var('river_direction', ('river', ), river_direction, 'River direction')

write_nc_var('river_Vshape', ('s_rho', 'river'), river_Vshape, 'River vertical shape')
write_nc_var('river_time', ('river_time', ), river_jd, time_units)
write_nc_var('river_transport', ('river_time', 'river'), river_transport,
             'River transport [m3/s]')
write_nc_var('river_temp', ('river_time', 's_rho', 'river'),
             river_temp,
             'River temperature [deg C]')
write_nc_var('river_salt', ('river_time', 's_rho', 'river'), river_salt,
             'River salinity [psu]')

nc.close()

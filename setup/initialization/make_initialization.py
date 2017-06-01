 # -*- coding: utf-8 -*-
'''
Make initial tracer netCDF files for simulation.
Using data from ODV, World Ocean Atlas 2013

run setup/make_initialization.py True False
'''

import netCDF4 as netCDF
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import octant
import xarray as xr
import gsw
import argparse


mpl.rcParams.update({'font.size': 14})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'


# parse the input arguments
parser = argparse.ArgumentParser()
parser.add_argument('doplot', type=int, default=1, help='1 to plot')
parser.add_argument('makeinitialfile', type=int, default=0, help='1 to make file for initialization')
args = parser.parse_args()

# doplot = True  # plot profile data
# makeinitialfile = True  # make netCDF file for initialization
doplot = bool(args.doplot)
makeinitialfile = bool(args.makeinitialfile)
print(doplot, makeinitialfile)

data = netCDF.Dataset('data_from_WOA13_1.00deg_1955-2012_Annual_centerGOM.nc')

# already chose station from sight in ODV
lon = -93.5+360
lat = 24.5

# find index of station
ilon = np.where(lon == data['longitude'][:])[0]
ilat = np.where(lat == data['latitude'][:])[0]
ind = list(set(ilat).intersection(ilon))[0]  # station index in netCDF file

depth = data['var1'][ind]
temp_present = data['var2'][ind]
salt_present = data['var3'][ind]

# change temp and salt to be inline with paleo values
temp = temp_present - 2  # two degrees colder
salt = salt_present + salt_present*0.035  # 3.5% more saline

# linear equation of state used in ROMS
TCOEF = 1.7e-4  # 1/Celsius
SCOEF = 7.6e-4  # 1/nondimensional
rho0_present = gsw.rho(35, 4.2, 0)  # kg/m^3
T0_present = 4.2  # deg C
S0_present = 35  # ppt. This is that same as in ROMS simulation.
rho_present = rho0_present*(1 + SCOEF*(salt_present-S0_present) - TCOEF*(temp_present-T0_present))
# calculate paleo background numbers from present numbers
T0 = T0_present - 2  # deg C
S0 = S0_present + S0_present*0.035  # ppt. This is that same as in ROMS simulation.
# read paleo density off T-S diagram since not sure it makes sense to mix
# present and paleo in linear EOS calculation
rho0 = gsw.rho(36.225, 2.2, 0)  # kg/m^3
# density of initial data profile
rho = rho0*(1 + SCOEF*(salt-S0) - TCOEF*(temp-T0))

# Calculate equivalent temperature values which with no salinity change
# will give the density profile in rho
tempeq = (1/TCOEF)*(1-rho/rho0) + T0

if doplot:

    fig, axes = plt.subplots(1, 2, figsize=(13.7,  8.15), sharey=True)
    color = '#934977'  # temp
    # paleo
    axes[0].plot(temp, depth, '-', color=color, lw=2)
    axes[0].plot(temp[::4], depth[::4], '*', color=color, ms=8)
    # present day
    axes[0].plot(temp_present, depth, '-', color=color, lw=1.0)
    axes[0].plot(temp_present[::4], depth[::4], '*', color=color, ms=5)
    axes[0].set_ylim(depth.max(), 0)
    axes[0].set_xlabel('($\star$) Temperature [C]', color=color)
    axes[0].tick_params('x', labelcolor=color, color=color)
    axes[0].set_ylabel('Depth [m]')
    axes[0].text(5, 3500, 'present day', fontsize=10, color='0.4',
                 rotation=90, verticalalignment='bottom')
    axes[0].text(1, 3500, 'paleo', fontsize=10, color='0.1',
                 rotation=90, verticalalignment='bottom')
    axes[0].text(15, 3500, 'paleo', fontsize=10, color='0.1',
                 rotation=90, verticalalignment='bottom')
    lims = axes[0].get_xlim()

    ax0twin = axes[0].twiny()
    color = '#497793'  # salt
    ax0twin.plot(salt, depth, ls='-', color=color, lw=2)
    ax0twin.plot(salt_present, depth, ls='-', color=color, lw=1.0)
    ax0twin.set_xlabel('Salinity [psu]', color=color)
    ax0twin.tick_params('x', labelcolor=color, color=color)

    color = '#934977'  # temp
    axes[1].plot(tempeq, depth, '-', color=color, lw=2)
    axes[1].plot(tempeq[::4], depth[::4], '*', color=color, ms=8, zorder=7)
    axes[1].set_ylim(depth.max(), 0)
    axes[1].set_xlabel(r'($\star$) Equivalent temperature [C]', color=color)
    axes[1].tick_params('x', labelcolor=color, color=color)
    axes[1].set_xlim(lims[0], lims[1])

    ax1twin = axes[1].twiny()
    color = '#779349'  # density
    ax1twin.plot(rho, depth, ls='-', color=color, lw=2, alpha=0.8)
    # ax1twin.plot(rho_present, depth, ls='-.', color=color, lw=0.5, alpha=0.7)
    ax1twin.set_xlabel(r'Paleo density [kg/m$^3$]', color=color)
    ax1twin.tick_params('x', labelcolor=color, color=color)
    ax1twin.get_xaxis().get_major_formatter().set_useOffset(False)

    plt.tight_layout()

    fig.savefig('WOA_profiles.pdf', bbox_inches='tight')
    fig.savefig('WOA_profiles.png', bbox_inches='tight')

if makeinitialfile:

    # following this for format:
    # https://www.myroms.org/index.php?page=initial

    grid = netCDF.Dataset('../../grid.nc')

    # coordinates
    lon_rho = grid['lon_rho']
    lat_rho = grid['lat_rho']
    ocean_time = 86400*14610  # from ocean_tglo.in

    # dimensions
    N = 20
    s_rho = octant.depths.get_srho(N)
    eta_rho = lon_rho.shape[0]
    xi_rho = lon_rho.shape[1]

    # find depths
    zr = octant.depths.get_zrho(Vtransform=2, Vstretching=4, N=N, theta_s=0.0001,
                                theta_b=0, h=grid['h'][:], hc=0, zeta=0,
                                Hscale=3)

    # use these depths to interpolate value from data profile
    T = np.interp(zr, depth, tempeq)

    # salinity initial file
    salt = np.ones_like(T)*S0
    u = np.zeros_like(salt)
    v = np.zeros_like(salt)

    # zeta, ubar, vbar, u, v
    zeta = np.zeros_like(lon_rho)
    ubar = np.zeros_like(lon_rho)
    vbar = np.zeros_like(lon_rho)

    # export to netcdf
    ds = xr.Dataset({'zeta': (['eta_rho', 'xi_rho'], zeta),
                     'ubar': (['eta_rho', 'xi_rho'], ubar),
                     'vbar': (['eta_rho', 'xi_rho'], vbar),
                     'u': (['s_rho', 'eta_rho', 'xi_rho'], u),
                     'v': (['s_rho', 'eta_rho', 'xi_rho'], v),
                     'temp': (['s_rho', 'eta_rho', 'xi_rho'], T),
                     'salt': (['s_rho', 'eta_rho', 'xi_rho'], salt)},
                    coords={'lon_rho': (['eta_rho', 'xi_rho'], lon_rho),
                            'lat_rho': (['eta_rho', 'xi_rho'], lat_rho),
                            'z_r': (['s_rho', 'eta_rho', 'xi_rho'], zr),
                            's_rho': (s_rho),
                            'ocean_time': (ocean_time)})  # ,
                            # 'reference_time': pd.Timestamp('2010-01-01')})
    ds.to_netcdf('ocean_ini.nc')

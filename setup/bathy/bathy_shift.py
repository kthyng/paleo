'''
Shift bathymetry up by 120 meters.

Need to change the grid mask also.
'''

import xarray as xr


ds = xr.open_dataset('grid-original.nc')
h = ds['h'].data

# remember min bathy
# dsmin = h.min()
dsmin = 25  # Set this deeper for more input

# shift bathymetry
h -= 120

# mask out anything that went negative since this is now land
ind = h < 0
mask_rho = ds['mask_rho'].data
mask_rho[ind] = 0

# make anything between 0 and dsmin deepened to dsmin
# also want all newly masked area to be dsmin to be consistent across mask
# in case that matters
ind = h < dsmin
h[ind] = dsmin
ds['h'].data = h  # update data structure

# fix rho mask where there is a hole
mask_rho[20, 176] = 0

# eliminate bottom right corner where salinity is weird
mask_rho[0, 255] = 0

# update other masks
ds['mask_u'].data = mask_rho[:, 1:]*mask_rho[:, :-1]
ds['mask_v'].data = mask_rho[1:, :]*mask_rho[:-1, :]
ds['mask_psi'].data = mask_rho[1:, 1:] * mask_rho[:-1, 1:] * \
    mask_rho[1:, :-1] * mask_rho[:-1, :-1]
ds['mask_rho'].data = mask_rho

# overwrite mask with new version
ds.drop('bath').to_netcdf('grid-shifted.nc', mode='w', format='NETCDF4', engine='netcdf4')

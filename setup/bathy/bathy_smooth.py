'''
Smooth bathymetry after original grid is shifted 120 meters in
bathy_shift.py.
'''

from bathy_smoother import LP_bathy_smoothing, bathy_smoothing
import xarray as xr
import numpy as np


smoothing = 'Positive'
rx0max = 0.4  # between 0 and 0.4

ds = xr.open_dataset('grid-shifted.nc')

Hobs = ds['h'].data
MSK = ds['mask_rho'].data
AreaMatrix = (1/ds['pm'].data)*(1/ds['pn'].data)

newgridfilename = 'grid-smoothed-' + smoothing + '-r' + str(rx0max) + '.nc'

if smoothing == 'lp':
    hnew = LP_bathy_smoothing.LP_smoothing_rx0(MSK, Hobs, rx0max,
                                               np.zeros_like(Hobs),
                                               np.ones_like(Hobs)*10000)

elif smoothing == 'Positive':
    hnew = bathy_smoothing.smoothing_Positive_rx0(MSK, Hobs, rx0max)

elif smoothing == 'Negative':
    hnew = bathy_smoothing.smoothing_Negative_rx0(MSK, Hobs, rx0max)

elif smoothing == 'PositiveVolume':
    hnew = bathy_smoothing.smoothing_PositiveVolume_rx0(MSK, Hobs, rx0max,
                                                        AreaMatrix)

elif smoothing == 'NegativeVolume':
    hnew = bathy_smoothing.smoothing_NegativeVolume_rx0(MSK, Hobs, rx0max,
                                                        AreaMatrix)

elif smoothing == 'PlusMinus':
    hnew, _, _ = bathy_smoothing.smoothing_PlusMinus_rx0(MSK, Hobs, rx0max,
                                                         AreaMatrix)

elif smoothing == 'Laplacian':
    hnew = bathy_smoothing.smoothing_Laplacian_rx0(MSK, Hobs, rx0max)

# match new bathymetry min to the masking bathy min
ind1 = ds['mask_rho'].data == 1  # water
ind0 = ds['mask_rho'].data == 0  # land
# Set masked depths to be equal to active min depth
hnew[ind0] = hnew[ind1].min()

# add in vertical parameters
ds['Vtransform'] = 2
ds['Vstretching'] = 4
ds['km'] = 20
ds['theta_s'] = 0.0001
ds['theta_b'] = 0.0
ds['hc'] = 0


ds['h'].data = hnew
ds.to_netcdf(newgridfilename, mode='w')

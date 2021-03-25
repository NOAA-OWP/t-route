import pathlib
import xarray as xr

domain_dir_old = pathlib.Path(
    '/glade/work/jamesmcc/domains/private/florence_933020089/NWM/DOMAIN')
domain_dir_new = pathlib.Path(
    '/glade/work/jamesmcc/domains/private/florence_cutout_v21/NWM/DOMAIN')

#-----------------------------------------------------------------------------
# Routelink

file_old = domain_dir_old / 'Route_Link.nc'
file_new = domain_dir_new / 'Route_Link.nc'

ds_old = xr.open_dataset(file_old)
ds_new = xr.open_dataset(file_new)

# for var in ds_new.variables:
#     print(var)

for var in ds_new.variables:
    if not ds_old[var].equals(ds_new[var]):
        print(var)

# ChSlp
# So
# ascendingIndex


#-----------------------------------------------------------------------------
# geo_em.d0?.nc

file_old = domain_dir_old / 'geo_em.d01.nc'
file_new = domain_dir_new / 'geo_em.d0x.nc'

ds_old = xr.open_dataset(file_old)
ds_new = xr.open_dataset(file_new)


for var in ds_new.variables:
    if not var in ds_old.variables:
        print(var)

# XLAT_C
# XLONG_C

# for var in ds_new.variables:
for var in ds_old.variables:
    if not ds_old[var].equals(ds_new[var]):
        print(var)

# LANDMASK
# LU_INDEX
# SCB_DOM
# SCT_DOM

import os
import xarray as xr

path = '/home/amin/Amin_fork/t-route/test/LowerColorado_TX/output/202108231300.flowveldepth.nc'


dataset = xr.open_dataset(path, engine='netcdf4')
print(dataset)
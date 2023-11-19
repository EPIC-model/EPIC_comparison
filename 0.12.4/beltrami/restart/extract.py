from tools.nc_reader import nc_reader
from tools.nc_fields import nc_fields
import numpy as np

ncr = nc_reader()
ncr.open('../ps3d_beltrami_64_fields.nc')

step = 60

time = ncr.get_all('t')[step]

ncst = nc_reader()
ncst.open('../ps3d_beltrami_64_field_stats.nc')
ke0 = ncst.get_all('ke')[0]
ke60 = ncst.get_all('ke')[step*10]
#print(ncst.get_all('t')[step*10])

print("Time:", time)
print("KE decay [percent]:", (ke0 - ke60) / ke0 * 100.0)

xi = ncr.get_dataset(step=step, name='x_vorticity', copy_periodic=False)
eta = ncr.get_dataset(step=step, name='y_vorticity', copy_periodic=False)
zeta = ncr.get_dataset(step=step, name='z_vorticity', copy_periodic=False)

origin = ncr.get_box_origin()
extent = ncr.get_box_extent()

xi = np.transpose(xi, axes=[2, 1, 0])
eta = np.transpose(eta, axes=[2, 1, 0])
zeta = np.transpose(zeta, axes=[2, 1, 0]) 

ncf = nc_fields()

ncf.open('restart_step_' + str(step) + '.nc')

ncf.set_time(time)

ncf.add_field('x_vorticity', xi, unit='1/s', long_name='x-vorticity')
ncf.add_field('y_vorticity', eta, unit='1/s', long_name='y-vorticity')
ncf.add_field('z_vorticity', zeta, unit='1/s', long_name='z-vorticity')

ncf.add_parameter('l_lower_boundry_zeta_zero', 'true')
ncf.add_parameter('l_upper_boundry_zeta_zero', 'true')

nz, ny, nx = xi.shape

ncf.add_box(origin, extent, [nx, ny, nz-1])


ncf.close()

from linear_interpl import find_bounds
from tools.nc_reader import nc_reader

#times = [5, 5.5, 6, 7.0, 8.0, 10.0]

ncr = nc_reader()

#ncr.open('rayleigh-taylor/epic_rt_384_fields.nc')

times = [600]

ncr.open('../data/moist/moist_fields.nc')

ts = ncr.get_all('t')

steps = []

for t in times:
    i, j = find_bounds(t, ts)

    steps = steps + [i, j]

    print(ts[i], t, ts[j])

ncr.close()

print(steps)

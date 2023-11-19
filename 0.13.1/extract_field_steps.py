# python extract_field_step.py --filename moist_fields.nc --step 12 --dataset humidity
import argparse
from tools.nc_reader import nc_reader
from tools.nc_fields import nc_fields
import numpy as np
import os

try:
    parser = argparse.ArgumentParser(
        description="Extract selected datasets from a field file for specific steps."
    )

    parser.add_argument('--filename',
                        type=str,
                        help='file path and file name')

    parser.add_argument('--datasets',
                        type=str,
                        nargs='+',
                        help='what kind of data to extract?')

    parser.add_argument('--steps',
                        type=int,
                        nargs='+',
                        help='steps to extract')

    args = parser.parse_args()

    if not os.path.exists(args.filename):
        raise IOError("File '" + args.filename + "' does not exist.")

    ncr = nc_reader()
    ncr.open(args.filename)


    ncf = nc_fields()

    dirname = os.path.dirname(args.filename)
    basename = os.path.basename(args.filename)

    ncf.open(os.path.join(dirname, 'extracted_steps_from_' + basename))

    t = ncr.get_all('t')

    origin = ncr.get_box_origin()
    extent = ncr.get_box_extent()
    ncells = ncr.get_box_ncells()
    ncf.add_box(origin, extent, ncells)

    steps = np.asarray(args.steps)

    ts = np.zeros(steps.size)

    for i, step in enumerate(args.steps):
        print("Extracting step", step)

        ts[i] = t[step]

        for name in args.datasets:

            dset = ncr.get_dataset(step=step, name=name, copy_periodic=False)

            # must get (z, y, x) ordering
            dset = np.transpose(dset, axes=[2, 1, 0])

            ncf.add_field(name, dset, unit='-', long_name=name, time_index=i)

    ncf.add_axis(axis='t', values=ts)

    x = ncr.get_all('x')
    ncf.add_axis(axis='x', values=x)
    y = ncr.get_all('y')
    ncf.add_axis(axis='y', values=y)
    z = ncr.get_all('z')
    ncf.add_axis(axis='z', values=z)


    ncf.close()
    ncr.close()

except Exception as ex:
    print(ex)


#4.990608288627607 5 5.000265192094116
#5.490053170888121 5.5 5.500048223439671
#5.990017875484795 6 6.000106183246031
#6.990143737620282 7.0 7.000070195093754
#7.990066435716057 8.0 8.000145284010056
#9.990192954625106 10.0 10.000108335991012
#[499, 500, 549, 550, 599, 600, 699, 700, 799, 800, 999, 1000]
# 499 500 549 550 599 600 699 700 799 800 999 1000

#python extract_field_steps.py --filename rayleigh-taylor/epic_rt_384_fields.nc --steps 499 500 549 550 599 600 699 700 799 800 999 1000 --datasets z_vorticity buoyancy


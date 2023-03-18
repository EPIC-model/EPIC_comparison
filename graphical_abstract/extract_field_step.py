# python extract_field_step.py --filename moist_fields.nc --step 12 --dataset humidity
import argparse
from tools.nc_reader import nc_reader
from tools.nc_fields import nc_fields
import numpy as np
import os

try:
    parser = argparse.ArgumentParser(
        description="Extract a single dataset from a field file for a specific step."
    )

    parser.add_argument('--filename',
                        type=str,
                        help='file path and file name')

    parser.add_argument('--dataset',
                        type=str,
                        help='what kind of data to extract?')
    
    parser.add_argument('--step',
                        type=int,
                        help='step to extract')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.filename):
        raise IOError("File '" + args.filename + "' does not exist.")

    ncr = nc_reader()
    ncr.open(args.filename)

    dset = ncr.get_dataset(step=args.step-1, name=args.dataset, copy_periodic=False)
        
    # must get (z, y, x) ordering
    dset = np.transpose(dset, axes=[2, 1, 0])
    
    ncf = nc_fields()

    dirname = os.path.dirname(args.filename)
    basename = os.path.basename(args.filename)
    
    ncf.open(os.path.join(dirname, 'extracted_step_' + str(args.step) + '_from_' + basename))

    ncf.add_field(args.dataset, dset, unit='-', long_name=args.dataset) 
    
    x = ncr.get_all('x')
    ncf.add_axis(axis='x', values=x)
    y = ncr.get_all('y')
    ncf.add_axis(axis='y', values=y)
    z = ncr.get_all('z')
    ncf.add_axis(axis='z', values=z)

    nz, ny, nx = dset.shape
    origin = ncr.get_box_origin()
    extent = ncr.get_box_extent()
    ncf.add_box(origin, extent, [nx, ny, nz-1])
    
    ncf.close()
    ncr.close()
    
except Exception as ex:
    print(ex)

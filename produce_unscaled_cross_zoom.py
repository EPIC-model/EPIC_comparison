from epic_utils import (
    produce_limited_cross_section,
    get_limited_coordinates
)
from utils import write_field_to_file
import xarray as xr

MAIN_DIR = "/work/e710/e710/shared/epic_comparison"
xz_loc = 3140
resolution = 256

# Set up file paths
epic_input_file_name = (
    MAIN_DIR + "/epic_rev_" + str(resolution) + "/moist_0000000012_parcels.nc"
)

boxes=[
        [1250,   2000,    5250,     5000],
        [3700,   3400,    4485,     4185],
        [3900,   3600,    4292.5,   3992.5],
        [4000,   3700,    4196.25,  3896.25],
        [4050,   3750,    4148.125, 3848.125],
        ]

for boxnr,box in enumerate(boxes):
    xp, yp, zp = get_limited_coordinates(epic_input_file_name,cross_type="xz",box=box,ndir=2048)
    xz_cross = produce_limited_cross_section(
        epic_input_file_name,
        xp,
        yp,
        zp,
        xz_loc,
        cross_type="xz",
        kernel="third",
        r_limit_fac=3.0,
    )
    write_field_to_file(
        xz_cross,
        xp,
        zp,
        "humidity",
        "epic_zooms_" + str(boxnr) + "_lim3r_t6.nc",
        cross_type="xz",
        mode="w",
    )
    xp, yp, zp = get_limited_coordinates(epic_input_file_name,cross_type="xz",box=box,ndir=4096)
    xz_cross = produce_limited_cross_section(
        epic_input_file_name,
        xp,
        yp,
        zp,
        xz_loc,
        cross_type="xz",
        kernel="sharp",
        r_limit_fac=1.0,
    )
    write_field_to_file(
        xz_cross,
        xp,
        zp,
        "humidity",
        "epic_zooms_" + str(boxnr) + "_sharp_t6.nc",
        cross_type="xz",
        mode="w",
    )




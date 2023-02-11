from epic_utils import (
    produce_epic_cross_section,
    get_epic_projection_coordinates,
    get_idealised_moist_parameters,
)
from mpic_utils import (
    produce_mpic_cross_section,
    get_mpic_projection_coordinates,
)
from utils import write_field_to_file
import xarray as xr

MAIN_DIR = "/home/stevenboeing/epic_comparison_klad"
RESOLUTIONS = [32]
xz_loc = 3.140
MPIC_RES_TSTEP_DICT = {32: "0099"}

# Combine these for a loop
METHODS = ["sharp", "lim3r", "lim5r"]
KERNELS = ["sharp", "third", "third"]
R_LIMIT_FACS = [1, 3, 5]
REFINEMENT = [64, 16, 16]

for resolution in RESOLUTIONS:
    # Set up file paths
    monc_input_file_name = (
        MAIN_DIR + "/monc_mpic_smooth_" + str(resolution) + "/mpic_diagnostic_3d_3.nc"
    )
    epic_input_file_name = (
        MAIN_DIR + "/epic_rev_" + str(resolution) + "/moist_0000000012_parcels.nc"
    )
    mpic_parcel_glob = (
        MAIN_DIR
        + "/pmpic_smooth_bubble_"
        + str(resolution)
        + "/parcels_?????_0"
        + MPIC_RES_TSTEP_DICT[resolution]
        + ".nc"
    )
    mpic_grid_file = (
        MAIN_DIR
        + "/pmpic_smooth_bubble_"
        + str(resolution)
        + "/grids_"
        + MPIC_RES_TSTEP_DICT[resolution]
        + ".nc"
    )
    epic_grid_file_name = (
        MAIN_DIR + "/epic_rev_" + str(resolution) + "/moist_fields.nc"
    )
    for mindex, method in enumerate(METHODS):
        # EPIC work
        xp, yp, zp = get_epic_projection_coordinates(
            epic_input_file_name,
            refinement_factor=REFINEMENT[mindex],
        )
        xz_cross = produce_epic_cross_section(
            epic_input_file_name,
            loc=xz_loc*1000.,
            cross_type="xz",
            kernel=KERNELS[mindex],
            refinement_factor=REFINEMENT[mindex],
            r_limit_fac=R_LIMIT_FACS[mindex],
        )
        len_condense, q_scale = get_idealised_moist_parameters(epic_input_file_name)
        write_field_to_file(
            xz_cross/q_scale,
            xp/len_condense,
            zp/len_condense,
            "humidity",
            "epic_" + str(resolution) + "_" + method + "_t6.nc",
            cross_type="xz",
            mode="w",
        )
        # MPIC work
        xp, yp, zp = get_mpic_projection_coordinates(
            mpic_grid_file, refinement_factor=REFINEMENT[mindex]
        )
        xz_cross = produce_mpic_cross_section(
            mpic_parcel_glob,
            mpic_grid_file,
            loc=xz_loc,
            cross_type="xz",
            kernel=KERNELS[mindex],
            refinement_factor=REFINEMENT[mindex],
            r_limit_fac=R_LIMIT_FACS[mindex],
        )
        write_field_to_file(
            xz_cross,
            xp,
            zp,
            "humidity",
            "mpic_" + str(resolution) + "_" + method + "_t6.nc",
            cross_type="xz",
            mode="w",
        )
    ds_nc = xr.open_dataset(monc_input_file_name)
    qt_monc= (ds_nc["q_cloud_liquid_mass"][0,:,:,1:]+ds_nc["q_vapour"][0,:,:,1:]).interp(y=resolution*xz_loc/6.28, method="nearest").data
    x_monc=(ds_nc["x"].data+0.5)*(6.28/resolution)
    z_monc=ds_nc["zn"].data[1:]
    write_field_to_file(
        qt_monc,
        x_monc,
        z_monc,
        "humidity",
        "monc_" + str(resolution) + "_t6.nc",
        cross_type="xz",
        mode="w",
    )

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

MAIN_DIR = "/work/e710/e710/shared/epic_comparison"
xz_loc = 3.140
MPIC_RES_TSTEP_DICT = {32: "0099", 64: "0193", 128: "0475", 256: "0941"}


def make_panels(RESOLUTIONS, METHODS, KERNELS, R_LIMIT_FACS, REFINEMENT):
    for resolution in RESOLUTIONS:
        # Set up file paths
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
        for mindex, method in enumerate(METHODS):
            # EPIC work
            xp, yp, zp = get_epic_projection_coordinates(
                epic_input_file_name,
                refinement_factor=REFINEMENT[mindex],
            )
            xz_cross = produce_epic_cross_section(
                epic_input_file_name,
                loc=xz_loc * 1000.0,
                cross_type="xz",
                kernel=KERNELS[mindex],
                refinement_factor=REFINEMENT[mindex],
                r_limit_fac=R_LIMIT_FACS[mindex],
            )
            len_condense, q_scale = get_idealised_moist_parameters(epic_input_file_name)
            write_field_to_file(
                xz_cross / q_scale,
                xp / len_condense,
                zp / len_condense,
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
            )  # Combine these for a loop


# RESOLUTIONS = [32, 64, 128, 256]
# METHODS = ["sharp", "lim3r"]
# KERNELS = ["sharp", "third"]
# R_LIMIT_FACS = [1, 3]
# REFINEMENT = [64, 16]

# make_panels(RESOLUTIONS,METHODS,KERNELS,R_LIMIT_FACS,REFINEMENT)

RESOLUTIONS = [32]
METHODS = ["lim5r"]
KERNELS = ["third"]
R_LIMIT_FACS = [5]
REFINEMENT = [16]

make_panels(RESOLUTIONS, METHODS, KERNELS, R_LIMIT_FACS, REFINEMENT)

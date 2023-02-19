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

MAIN_DIR = "/work/e710/e710/shared/epic_comparison"
xz_loc = 3140


def make_panels(RESOLUTIONS, METHODS, KERNELS, R_LIMIT_FACS, REFINEMENT):
    for resolution in RESOLUTIONS:
        # Set up file paths
        epic_input_file_name = (
            MAIN_DIR + "/epic_rev_" + str(resolution) + "/moist_0000000012_parcels.nc"
        )
        for mindex, method in enumerate(METHODS):
            # EPIC work
            xp, yp, zp = get_epic_projection_coordinates(
                epic_input_file_name,
                refinement_factor=REFINEMENT[mindex],
            )
            xz_cross = produce_epic_cross_section(
                epic_input_file_name,
                loc=xz_loc,
                cross_type="xz",
                kernel=KERNELS[mindex],
                refinement_factor=REFINEMENT[mindex],
                r_limit_fac=R_LIMIT_FACS[mindex],
            )
            len_condense, q_scale = get_idealised_moist_parameters(epic_input_file_name)
            write_field_to_file(
                xz_cross,
                xp,
                zp,
                "humidity",
                "epic_unsc_" + str(resolution) + "_" + method + "_t6.nc",
                cross_type="xz",
                mode="w",
            )


RESOLUTIONS = [256]
METHODS = ["lim3r"]
KERNELS = ["third"]
R_LIMIT_FACS = [3]
REFINEMENT = [16]

make_panels(RESOLUTIONS, METHODS, KERNELS, R_LIMIT_FACS, REFINEMENT)

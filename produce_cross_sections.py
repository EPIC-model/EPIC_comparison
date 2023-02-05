from utils import produce_epic_cross_section, write_field_to_file, get_projection_coordinates
from pylab import *

MAIN_DIR = "/home/stevenboeing/epic_comparison_klad"
RESOLUTIONS = [32]
xz_loc = 3140

for resolution in RESOLUTIONS:
    input_file_name=MAIN_DIR+"/epic_rev_" + str(resolution) + "/moist_0000000012_parcels.nc"
    # normal cross-section
    xp, yp, zp = get_projection_coordinates(input_file_name, refinement_factor=16)
    xz_cross = produce_epic_cross_section(
        input_file_name,
        loc=xz_loc,
        cross_type="xz",
        kernel="third",
        refinement_factor=16,
    )
    write_field_to_file(xz_cross, xp, zp, "humidity", "epic_"+str(resolution)+"_t600.nc", cross_type="xz", mode="w")
    # sharp cross-section
    # ~ xp, yp, zp = get_projection_coordinates(input_file_name, refinement_factor=64)
    # ~ xz_cross = produce_epic_cross_section(
        # ~ input_file_name,
        # ~ loc=xz_loc,
        # ~ cross_type="xz",
        # ~ kernel="sharp",
        # ~ refinement_factor=64,
    # ~ )
    # ~ write_field_to_file(xz_cross, xp, zp, "humidity", "epic_"+str(resolution)+"_t600_sharp.nc", cross_type="xz", mode="w")

import numpy as np
import netCDF4 as nc
import xarray as xr
from glob import glob

MAIN_DIR = "/home/stevenboeing/epic_comparison_klad"
RESOLUTIONS = [32]
xz_loc = 3.140
MPIC_RES_TSTEP_DICT = {32: "0099"}

nbins=100
bin_edges=np.linspace(0.,0.10,nbins)
bin_centres=0.5*(bin_edges[1:]+bin_edges[:-1])

bin_edges_coord = {
    "bin_edges": ("bin_edges", bin_edges, {"long_name": "Bin edges"})
}
bin_centres_coord = {
    "bin_centres": ("bin_centres", bin_centres, {"long_name": "Bin centres"})
}

ds = xr.Dataset(
    coords={
        **bin_edges_coord,
        **bin_centres_coord,
    }
)

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
    ds_nc = nc.Dataset(epic_input_file_name)
    vol = ds_nc.variables["volume"][0,:]
    h = ds_nc.variables["humidity"][0,:]
    z = ds_nc.variables["z_position"][0,:]
    len_condense = ds_nc.groups["physical_quantities"].scale_height
    q_scale = ds_nc.groups[
        "physical_quantities"
    ].saturation_specific_humidity_at_ground_level
    hl=(h/q_scale-np.exp(-z/len_condense))
    hl=hl*(hl>0.0)
    hist_epic,bins_epic=np.histogram(hl,weights=vol,density=True,bins=bin_edges)
    ds['hist_epic_'+str(resolution)]=(
                ("bin_centres"),
                hist_epic,
            )
    ds['hist_epic_'+str(resolution)].attrs={"long_name": "EPIC histogram for $"+str(resolution)+"^3$ points", "units": "-"}
    hl_all=np.array([])
    vol_all=np.array([])
    for input_file_name in glob(mpic_parcel_glob):
        ds_nc = nc.Dataset(input_file_name)
        vol = ds_nc.variables["vol"][:]
        h = ds_nc.variables["h"][:]
        z = ds_nc.variables["z"][:]
        hl=h-np.exp(-z)
        hl=hl*(hl>0.0)
        hl_all=np.hstack((hl_all,hl))
        vol_all=np.hstack((vol_all,vol))
    hist_mpic,bins_mpic=np.histogram(hl_all,weights=vol_all,density=True,bins=bin_edges)
    ds['hist_mpic_'+str(resolution)]=(
                ("bin_centres"),
                hist_mpic,
            )
    ds['hist_mpic_'+str(resolution)].attrs={"long_name": "MPIC histogram for $"+str(resolution)+"$^3$ points", "units": "-"}
    ds_nc = nc.Dataset(monc_input_file_name)
    hl = ds_nc.variables["q_cloud_liquid_mass"][:,:,:,1:]
    hist_monc,bins_monc=np.histogram(hl,density=True,bins=bin_edges)
    ds['hist_monc_'+str(resolution)]=(
                ("bin_centres"),
                hist_monc,
            )
    ds['hist_monc_'+str(resolution)].attrs={"long_name": "MONC histogram for $"+str(resolution)+"$^3$ points", "units": "-"}
ds.to_netcdf('pdfs.nc')



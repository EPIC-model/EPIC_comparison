from numba import njit
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import xarray as xr
import glob
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib.lines import Line2D

# Calculates equilibrium size distribution for "random splitting and merging" distibution

plt.rcParams.update(
    {
        "figure.figsize": (8, 4),
        "figure.dpi": 200,
        "font.family": "serif",
        "font.size": 12,
        "text.usetex": True,
        'text.latex.preamble': "\n".join([
        r"\usepackage{amsmath}",
        r"\usepackage[utf8]{inputenc}",
        r"\usepackage[T1]{fontenc}",
        r"\usepackage{siunitx}",
        ])
})


@njit
# Calculate exponential bin index, for a given volume
# Exponentially spaced volumes between V_min and V_max are used
# Uses the sum of a geometric sequence to determine the index
def calc_bin(V, V_min, V_bin_min, bin_scaling):
    """
    V_min: minimum parcel volume
    V_bin_min: minimum size spacing
    bin_scaling: size increase factor per bin
    """
    return int(
        np.log(1.0 - ((V - V_min) / V_bin_min) * (1.0 - bin_scaling))
        / np.log(bin_scaling)
    )


@njit
def calc_distribution(
    V_min,
    V_max,
    V_init,
    n_tot=100000,
    t_start_scale=10000,
    t_end_scale=20000,
    bin_factor=20,
    bin_scaling=1.05,
):
    """
    V_min, V_max: minimal and maximal parcel size
    V_init: mean initial volume (we may want to code up a non-random version as well, for better analogy to EPIC)
    n_tot: total number of parcels
    t_start_scale: number of "superiterations" when sampling starts (a superiteration consists of n_tot random splitting events)
    t_end_scale: number of "superiterations" when sampling ends
    bin_factor: fractional size (in terms of V_min) of smallest bin
    bin_scaling: size increase factor per bin
    """
    # Sanity_checks
    if 2.0 * V_init - V_min > V_max:
        raise Exception(
            "Error: 2.0*V_init-V_min>V_max (meaning V may be initialised larger than V_max)"
        )
    if V_max < 2.0 * V_min:
        raise Exception("Error: V_max<2.0*V_min (insufficient dynamic range)")
    # Set up default for V_max if e.g. zero provided and volume table
    if V_max < V_min:
        V_max = V_min * n_tot
    # Set up bins for pdf
    V_bin_min = V_min / bin_factor
    max_bin_index = calc_bin(V_max, V_min, V_bin_min, bin_scaling)
    counts = np.zeros(
        max_bin_index + 1
    )  # Add 1, so there is an index corresponding to V_max
    # Initialise volumes
    V = np.zeros(n_tot + 1)  # Include spare index for while loop below
    V_sum = 0.0
    n_init = 0
    while V_sum < V_min * n_tot:
        V_temp = V_min + 2.0 * (V_init - V_min) * np.random.rand()
        V[n_init] = V_temp
        V_sum = V_sum + V_temp
        n_init = n_init + 1
    # Reset last parcel, which took the volume over V_min*n_tot
    V[n_init] = 0.0
    n_parcels = n_init - 1
    # Set up the problem
    for n_cycle in range(t_end_scale):
        for operation_in_cycle in range(n_tot):
            # Randomly selects a parcels, split it, and merge if needed
            this_parcel = int(n_parcels * np.random.rand())
            V_temp = 0.5 * V[this_parcel]
            if V_temp > V_min:
                # Just split the parcel, and add a parcel
                V[this_parcel] = V_temp
                n_parcels = n_parcels + 1  # First increase n_parcels
                V[n_parcels - 1] = V_temp  # Then add at last index
            else:
                # Backfill from last active parcel
                V[this_parcel] = V[n_parcels - 1]
                n_parcels = n_parcels - 1
                # Merge the split parcels into random other parcels
                rand1 = int(np.random.rand() * n_parcels)
                V[rand1] = V[rand1] + V_temp
                rand2 = int(np.random.rand() * n_parcels)
                V[rand2] = V[rand2] + V_temp
                # Split if these are larger than Vmax
                if V[rand1] > V_max:
                    V_temp = 0.5 * V[rand1]
                    V[rand1] = V_temp
                    n_parcels = n_parcels + 1
                    V[n_parcels - 1] = V_temp
                if V[rand2] > V_max:
                    V_temp = 0.5 * V[rand2]
                    V[rand2] = V_temp
                    n_parcels = n_parcels + 1
                    V[n_parcels - 1] = V_temp
        # Collect statistics
        if n_cycle > t_start_scale:
            for n_parcel in range(n_parcels):
                bin_index = calc_bin(V[n_parcel], V_min, V_bin_min, bin_scaling)
                counts[bin_index] += 1
    # Construct pdf information to return
    lower_bin_edges = np.zeros(len(counts) + 1)
    # Restrict edges to Vmax, for proper normalisation
    for n_bin in range(len(lower_bin_edges)):
        lower_bin_edges[n_bin] = min(
            V_min + V_bin_min * (1 - bin_scaling**n_bin) / (1 - bin_scaling), V_max
        )
    bin_centres = 0.5 * (lower_bin_edges[1:] + lower_bin_edges[:-1])
    bin_widths = lower_bin_edges[1:] - lower_bin_edges[:-1]
    total_counts = sum(counts)
    densities = (counts / total_counts) / bin_widths
    densities = densities / (densities > 0.0)  # Select valid indices
    return bin_centres, densities


V_data = np.array([])
for this_file in ["epic_beltrami_32_0000000006_parcels.nc"]:
    this_data = xr.open_dataset(this_file)
    V_data = np.append(V_data, this_data["volume"].data[0], axis=0)
V_min = (np.pi**3 / 32.0**3) / 40.0
V_data = np.array(V_data) / V_min
print(len(V_data), np.nanmax(V_data))

fig = plt.figure()
ax = fig.add_subplot(111)

ax.hist(V_data, bins=np.linspace(0, 10, 201), density=True, histtype="step", zorder=-1, color='b')

ax.hist(
    V_data[V_data >= 1.0],
    bins=np.linspace(0, 10, 201),
    density=True,
    histtype="step",
    zorder=-1,
    color='r',
)

colors=['k','m','c']
V_min = 1.0
V_init = 2.0
for nr,V_max in enumerate([32.0, 8.0, 4.0]):
    bin_centres, densities = calc_distribution(
        V_min,
        V_max,
        V_init,
        n_tot=100000,
        t_start_scale=4000,
        t_end_scale=8000,
        bin_factor=20,
        bin_scaling=1.05,
    )
    ax.scatter(bin_centres, densities, s=2.0, zorder=0, color=colors[nr])

ax.set_xlabel(r'$V/V_{\min}$')
ax.set_ylabel(r'probability density')
ax.set_xlim(0.,10.)
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(MultipleLocator(0.5))

legend=ax.legend(
    (
        r'EPIC Beltrami',
        r'EPIC Beltrami ($V \geq V_{\min}$ only)',
        r'Conceptual $V_{\max}/V_{\min} = 32$',
        r'Conceptual $V_{\max}/V_{\min} = 8$',
        r'Conceptual $V_{\max}/V_{\min} = 4$',
    ),
)
handles = legend.legendHandles
handles[0]= Line2D([0], [0], color='b')
handles[1]= Line2D([0], [0], color='r')
legend=ax.legend(
    handles,
    (
        r'EPIC Beltrami',
        r'EPIC Beltrami ($V \geq V_{\min}$ only)',
        r'Conceptual $V_{\max}/V_{\min} = 32$',
        r'Conceptual $V_{\max}/V_{\min} = 8$',
        r'Conceptual $V_{\max}/V_{\min} = 4$',
    ),
)
plt.savefig("epic_distro.pdf",bbox_inches='tight')

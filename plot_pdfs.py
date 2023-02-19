import matplotlib.pyplot as plt
import xarray as xr
from cycler import cycler

plt.rcParams.update(
    {
        "figure.figsize": (8, 5),
        "figure.dpi": 200,
        "font.family": "serif",
        "font.size": 11,
        "text.usetex": False,
    }
)
#    "text.usetex": True,
#    'text.latex.preamble': "\n".join([
#        r"\usepackage{amsmath}",
#        r"\usepackage[utf8]{inputenc}",
#        r"\usepackage[T1]{fontenc}",
#        r"\usepackage{siunitx}",
#        ])
# })

RESOLUTIONS = [32, 64, 128, 256]
default_cycler = cycler(
    color=[
        "#ff8080",
        "#ff80ff",
        "#9999ff",
        "#ff3333",
        "#ff1aff",
        "#3333ff",
        "#cc0000",
        "#b300b3",
        "#0000cc",
        "#660000",
        "#4d004d",
        "#000066",
    ]
)
plt.rc("axes", prop_cycle=default_cycler)
ds = xr.open_dataset("pdfs.nc")

fig = plt.Figure(figsize=(7, 15))
for resolution in RESOLUTIONS:
    plt.stairs(ds["hist_epic_" + str(resolution)], ds["bin_edges"])
    plt.stairs(ds["hist_mpic_" + str(resolution)], ds["bin_edges"])
    plt.stairs(ds["hist_monc_" + str(resolution)], ds["bin_edges"])
plt.ylim(0, 2.5)
plt.xlim(0, 0.1)
plt.xlabel("$q_l$ [kg/kg]")
plt.ylabel("probability density")
plt.legend(["EPIC", "MPIC", "MONC"])
fig.tight_layout()
plt.savefig("compare_hl_pdf.png")
plt.savefig("compare_hl_pdf.jpeg")

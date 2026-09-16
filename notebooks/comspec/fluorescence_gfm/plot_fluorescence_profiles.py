"""Overview figure of the reconstructed fluorescence profiles: fig/comspec/fluorescence_profiles.png."""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))   # -> project
sys.path.insert(0, os.path.join(ROOT, "notebooks"))
import rcparams  # noqa: F401,E402  (project figure style)

OUT = os.path.join(ROOT, "data", "fluorescence")
R_DISPLAY = 300         # smoothing resolving power for the overview panels


def smooth(prof, col, R_grid=5000):
    sigma = R_grid / R_DISPLAY / 2.355
    return gaussian_filter1d(prof[col].to_numpy(), sigma)


def main():
    P = pd.read_csv(os.path.join(OUT, "profiles", "all_species_T070K.csv"))
    lam = P["lam_um"].to_numpy()
    groups = [("Main species, T_rot = 70 K", ["H2O", "CO2", "CO"]),
              ("Organics, T_rot = 70 K", ["CH4", "C2H6", "CH3OH", "C2H2", "H2CO", "C2H4"]),
              ("N- and S-bearing, T_rot = 70 K", ["NH3", "HCN", "HC3N", "OCS", "H2S"])]
    fig, axes = plt.subplots(4, 1, figsize=(20, 26))
    for ax, (title, sps) in zip(axes[:3], groups):
        for s in sps:
            y = smooth(P, s)
            ax.plot(lam, np.where(y > 0, y, np.nan), lw=2, label=s)
        ax.set_yscale("log"); ax.set_ylim(1e-7, 3e-2); ax.set_xlim(0.7, 5.0)
        ax.set_title(title); ax.set_ylabel(r"$g_\lambda$ [photons s$^{-1}$ molec$^{-1}$ µm$^{-1}$]")
        ax.legend(ncol=3, loc="upper left")
    ax = axes[3]
    H = pd.read_csv(os.path.join(OUT, "profiles", "H2O_gprofile.csv"))
    for T, c in zip((30, 70, 130), ("tab:blue", "tab:green", "tab:red")):
        y = smooth(H, f"g_lam_T{T:03d}K")
        ax.plot(H["lam_um"], np.where(y > 0, y, np.nan), lw=2, color=c, label=f"H$_2$O, {T} K")
    C = pd.read_csv(os.path.join(OUT, "profiles", "CO_gprofile.csv"))
    y = smooth(C, "g_lam_T070K"); ax.plot(C["lam_um"], np.where(y > 0, y, np.nan), lw=2, color="k", ls="--", label="CO, 70 K")
    ax.set_yscale("log"); ax.set_ylim(1e-6, 1e-2); ax.set_xlim(2.4, 5.0)
    ax.set_title("H$_2$O rotational-temperature dependence and the CO / H$_2$O hot-band blend")
    ax.set_ylabel(r"$g_\lambda$ [photons s$^{-1}$ molec$^{-1}$ µm$^{-1}$]"); ax.set_xlabel("wavelength [µm]")
    ax.legend(ncol=2, loc="upper left")
    for ax in axes[:3]:
        ax.set_xlabel("wavelength [µm]")
    fig.suptitle("Reconstructed GSFC-style fluorescence g-factors at 1 au (smoothed to R = 300)", y=0.995)
    fig.tight_layout()
    os.makedirs(os.path.join(ROOT, "fig", "comspec"), exist_ok=True)
    fig.savefig(os.path.join(ROOT, "fig", "comspec", "fluorescence_profiles.png"))
    print("saved fig/comspec/fluorescence_profiles.png")


if __name__ == "__main__":
    main()

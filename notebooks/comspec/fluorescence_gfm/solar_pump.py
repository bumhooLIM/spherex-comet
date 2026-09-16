"""Solar pumping spectrum at 1 au (Villanueva et al. 2011, Appendix B recipe).

F_nu(nu) = [Kurucz ATLAS9 line-free solar continuum, scaled to 1 au] x [disk-integrated solar
pseudo-transmittance of Toon (JPL, 2024 merged ATMOS/ACE/Kitt Peak atlas), 600-33300 cm-1 at
0.01 cm-1].  The product reproduces the TSIS-1 Hybrid Solar Reference Spectrum to +1.5 % over
0.7-2.7 um (see data/fluorescence/README.md).
"""
import os
import numpy as np
import pandas as pd

C_LIGHT = 2.99792458e8
KURUCZ_TO_IRR = 2.720e-4        # 4 pi (Rsun/AU)^2: erg cm-2 s-1 sr-1 nm-1 -> W m-2 um-1 (00asun.readme)
KURUCZ_FILE = "fsunallp.10000resam25"
TOON_FILE = "toon_solar_merged_20240731_600_33300_100.out.gz"
TSIS_FILE = "tsis1_hsrs_1nm.csv"


def load_kurucz(inputs):
    a = np.loadtxt(os.path.join(inputs, "solar", KURUCZ_FILE))
    wl_nm, flux, cont = a[:, 0], a[:, 1], a[:, 2]
    return wl_nm, flux * KURUCZ_TO_IRR, cont * KURUCZ_TO_IRR       # W m-2 um-1 at 1 au, vacuum nm


def load_toon(inputs):
    npz = os.path.join(inputs, "solar", "toon_cache.npz")
    if os.path.exists(npz):
        d = np.load(npz)
        return d["nu"], d["tr"].astype(float)
    arr = np.loadtxt(os.path.join(inputs, "solar", TOON_FILE), skiprows=3)
    nu, tr = arr[:, 0], arr[:, 1].astype(np.float32)
    np.savez_compressed(npz, nu=nu, tr=tr)
    return nu, tr.astype(float)


def load_tsis_1nm(inputs):
    d = pd.read_csv(os.path.join(inputs, "solar", TSIS_FILE))
    return d.iloc[:, 0].to_numpy(), d.iloc[:, 1].to_numpy() * 1e3    # nm, W m-2 um-1


class SolarPump:
    """Callable F_nu(nu_cm) [W m-2 Hz-1] at 1 au."""

    def __init__(self, inputs, scale=1.0, use_lines=True):
        wl_nm, flux, cont = load_kurucz(inputs)
        self.nu_k = 1e7 / wl_nm[::-1]
        self.cont_flam = cont[::-1] * scale
        self.flux_flam = flux[::-1] * scale
        self.nu_t, self.tr = load_toon(inputs)
        self.use_lines = use_lines

    def continuum_flam(self, nu_cm):
        return np.interp(nu_cm, self.nu_k, self.cont_flam)

    def transmittance(self, nu_cm):
        if not self.use_lines:
            return np.ones_like(np.asarray(nu_cm, float))
        return np.interp(nu_cm, self.nu_t, self.tr, left=1.0, right=1.0)

    def flam(self, nu_cm):
        return self.continuum_flam(nu_cm) * self.transmittance(nu_cm)      # W m-2 um-1

    def fnu(self, nu_cm):
        nu_cm = np.asarray(nu_cm, float)
        lam_m = 1e-2 / nu_cm
        return self.flam(nu_cm) * 1e6 * lam_m ** 2 / C_LIGHT

    __call__ = fnu

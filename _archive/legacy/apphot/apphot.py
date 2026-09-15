import re
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.stats import SigmaClip
from photutils.aperture import (
    CircularAperture, 
    CircularAnnulus, 
    ApertureStats, 
    aperture_photometry
)
from skimage.draw import disk  # Critical missing import
from typing import Any, List, Tuple, Union, Dict, Optional
import argparse
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.gridspec import GridSpec
from astropy.time import Time
from astropy.io import fits
from astropy.visualization import ZScaleInterval
import directory  # type: ignore
import rcparams  # type: ignore
import time  # Added for overall execution timin
import warnings

def reconstruct_cutout_wcs(row: Any) -> WCS:
    """
    Reconstructs a full Astropy WCS object for a SPHEREx cutout.
    Optimized to accept NamedTuples (from df.itertuples()) for maximum speed.
    """
    hdr = fits.Header()
    
    # 2. Base Coordinate System (using getattr for NamedTuple compatibility)
    hdr['WCSAXES'] = getattr(row, 'WCSAXES', 2)
    hdr['RADESYS'] = getattr(row, 'RADESYS', 'ICRS')
    hdr['CTYPE1']  = getattr(row, 'CTYPE1')
    hdr['CTYPE2']  = getattr(row, 'CTYPE2')
    hdr['CUNIT1']  = getattr(row, 'CUNIT1', 'deg')
    hdr['CUNIT2']  = getattr(row, 'CUNIT2', 'deg')
    hdr['LONPOLE'] = getattr(row, 'LONPOLE')
    hdr['LATPOLE'] = getattr(row, 'LATPOLE')
    
    # 3. Coordinate Reference Values
    hdr['CRVAL1'] = getattr(row, 'CRVAL1')
    hdr['CRVAL2'] = getattr(row, 'CRVAL2')
    hdr['CDELT1'] = getattr(row, 'CDELT1', 1.0)
    hdr['CDELT2'] = getattr(row, 'CDELT2', 1.0)
    
    # 4. Shift the Reference Pixel (CRPIX) using LTV offsets
    hdr['CRPIX1'] = getattr(row, 'CRPIX1') + getattr(row, 'ltv1')
    hdr['CRPIX2'] = getattr(row, 'CRPIX2') + getattr(row, 'ltv2')
    
    # 5. Rotation/Pixel Scale Matrix (PC Matrix)
    hdr['PC1_1'] = getattr(row, 'PC1_1')
    hdr['PC1_2'] = getattr(row, 'PC1_2')
    hdr['PC2_1'] = getattr(row, 'PC2_1')
    hdr['PC2_2'] = getattr(row, 'PC2_2')
    
    # 6. Inject the Physical Image Bounds (NAXIS)
    cutout_size = int(getattr(row, 'cutout_size'))
    hdr['NAXIS']  = 2
    hdr['NAXIS1'] = cutout_size
    hdr['NAXIS2'] = cutout_size
    
    # 7. Dynamically Scoop all SIP Distortion Coefficients
    sip_pattern = re.compile(r"^[AB]P?(_ORDER|_\d_\d)$")
    
    # Check if we are dealing with a NamedTuple (_fields) or a Series (index)
    fields = getattr(row, '_fields', row.index if hasattr(row, 'index') else [])
    
    for col_name in fields:
        val = getattr(row, col_name)
        if sip_pattern.match(col_name) and pd.notna(val):
            hdr[col_name] = val
            
    return WCS(hdr, relax=True)


def spherex_phase_grouping(fits_summary, threshold_days=28, key_dateobs='DATE-OBS'):
    """
    Groups FITS file records into discrete observation phases based on the time gap 
    between consecutive observations.

    This function sorts the input DataFrame by the specified observation date column. 
    It then calculates the time difference between consecutive rows and assigns a new 
    phase number whenever the gap exceeds the defined `threshold_days`.

    Args:
        fits_summary (pandas.DataFrame): A DataFrame containing summary information 
            of FITS files, which must include a column with observation dates.
        threshold_days (int or float, optional): The maximum allowable time gap (in days) 
            between consecutive observations to remain in the same phase. A gap strictly 
            greater than this value triggers a new phase. Defaults to 28.
        key_dateobs (str, optional): The column name in `fits_summary` that contains 
            the observation datetimes. Defaults to 'DATE-OBS'.

    Returns:
        pandas.DataFrame: A copy of the input DataFrame, sorted chronologically by 
            `key_dateobs`, with an additional integer column named 'PHASE' (starting 
            from 1) indicating the grouped phase of each observation.
    """
    
    fits_summary_updated = fits_summary.copy()
    dateobs = pd.to_datetime(fits_summary_updated[key_dateobs])
    
    sort_idx = dateobs.argsort()
    fits_summary_updated = fits_summary_updated.iloc[sort_idx].reset_index(drop=True)
    dateobs = dateobs.iloc[sort_idx].reset_index(drop=True)
    
    time_diff = dateobs.diff()
    is_new_phase = time_diff > pd.Timedelta(days=threshold_days)
    fits_summary_updated['phase'] = is_new_phase.cumsum() + 1
    
    return fits_summary_updated

def create_gaia_subset(
    gaia_all: Any, 
    fits_summary: pd.DataFrame, 
    del_ra_deg: float, 
    del_dec_deg: float, 
    gmag_limit: float, 
    ra_col: str = 'ra', 
    dec_col: str = 'dec'
) -> Any:
    """
    Filters a massive Gaia catalog for multiple overlapping target regions.
    Optimized via magnitude cuts, global bounding boxes, unique pointings, 
    and dynamic spherical RA adjustments.
    """
    global_mask = gaia_all['phot_g_mean_mag'] < gmag_limit
    valid_coords = fits_summary[[ra_col, dec_col]].dropna()
    if valid_coords.empty:
        return gaia_all[0:0] 
        
    dec_min_global = valid_coords[dec_col].min() - del_dec_deg
    dec_max_global = valid_coords[dec_col].max() + del_dec_deg
    global_mask &= (gaia_all['dec'] >= dec_min_global) & (gaia_all['dec'] <= dec_max_global)
    
    gaia_filtered = gaia_all[global_mask]
    ra_gaia = gaia_filtered['ra']
    dec_gaia = gaia_filtered['dec']
    
    exact_spatial_mask = np.zeros(len(gaia_filtered), dtype=bool)
    grid_res = min(del_ra_deg, del_dec_deg) * 0.5
    
    unique_pointings = (valid_coords / grid_res).round() * grid_res
    unique_pointings = unique_pointings.drop_duplicates()
    
    eff_del_ra_base = del_ra_deg + (grid_res / 2.0)
    eff_del_dec = del_dec_deg + (grid_res / 2.0)
    
    for ra_obj, dec_obj in unique_pointings.itertuples(index=False):
        local_dec_mask = (dec_gaia >= dec_obj - eff_del_dec) & (dec_gaia <= dec_obj + eff_del_dec)
        
        if abs(dec_obj) > 89.0:
            eff_del_ra_local = 180.0 
        else:
            cos_dec = np.cos(np.radians(dec_obj))
            eff_del_ra_local = min(180.0, eff_del_ra_base / cos_dec)
            
        ra_min = (ra_obj - eff_del_ra_local) % 360
        ra_max = (ra_obj + eff_del_ra_local) % 360
        
        if ra_min < ra_max:
            local_ra_mask = (ra_gaia >= ra_min) & (ra_gaia <= ra_max)
        else:
            local_ra_mask = (ra_gaia >= ra_min) | (ra_gaia <= ra_max)
            
        exact_spatial_mask |= (local_ra_mask & local_dec_mask)
        
    return gaia_filtered[exact_spatial_mask]


def flag_to_mask(flag, flag_number):
    mask = np.zeros_like(flag, dtype=bool)
    for f in flag_number:
        mask |= (flag & (1 << f)) != 0
    return mask.astype(bool)


def gaia_to_mask(gaia_subset, sci, gmag_limit, mag_col='phot_g_mean_mag', radius_scale=1.0):
    """
    Generates a boolean mask array where circular regions around Gaia sources are masked (True).
    The radius of each mask depends on the star's magnitude and the image's PSF FWHM.

    Args:
        gaia_subset (numpy.ndarray or pandas.DataFrame): The filtered catalog of Gaia sources.
        sci (astropy.io.fits.ImageHDU): The science image HDU.
        gmag_limit (float): The limiting magnitude. Fainter stars will have smaller 
            (or zero) mask radii.
        mag_col (str, optional): The column name for G-band magnitude. Defaults to 'phot_g_mean_mag'.
        radius_scale (float, optional): The scaling factor for the mask radius. Defaults to 1.0.
        dec_col (str, optional): The column name for Declination. Defaults to 'dec'.
        mag_col (str, optional): The column name for G-band magnitude. Defaults to 'phot_g_mean_mag'.
        radius_scale (float, optional): The scaling factor for the mask radius. Defaults to 0.3.

    Returns:
        numpy.ndarray: A 2D boolean array of shape (ny, nx) where True indicates a masked star pixel.
    """
    
    hdr = sci.header
    ny, nx = hdr['NAXIS2'], hdr['NAXIS1']
    wcs = WCS(hdr)
    
    psf_fwhm_arcsec = hdr.get('PSF_FWHM', 6.0)
    psf_fwhm_pixel = psf_fwhm_arcsec / hdr.get('PIX-SCL', 6.0) 
    
    mask_source = np.zeros((ny, nx), dtype=bool)
    
    ra_gaia   = np.array(gaia_subset['ra'], dtype=float)
    dec_gaia  = np.array(gaia_subset['dec'], dtype=float)
    gmag_gaia = np.array(gaia_subset[mag_col], dtype=float)
    
    skycoords_gaia = SkyCoord(ra=ra_gaia*u.deg, dec=dec_gaia*u.deg)
    x_gaia, y_gaia = wcs.world_to_pixel(skycoords_gaia)
    
    # Vectorized radius calculation to skip the Python loop for faint stars
    radii = radius_scale * psf_fwhm_pixel * (gmag_limit - gmag_gaia)
    valid_idx = (radii > 0) & ~np.isnan(x_gaia) & ~np.isnan(y_gaia)
    
    for xc, yc, radius in zip(x_gaia[valid_idx], y_gaia[valid_idx], radii[valid_idx]):
        rr, cc = disk((yc, xc), radius, shape=(ny, nx))
        mask_source[rr, cc] = True
        
    return mask_source


def perform_aperture_photometry(sci_data, err_data, master_mask, xcen, ycen, r_ap_list, ap_in_out):
    """
    Performs robust, sub-pixel accurate aperture photometry for multiple aperture sizes simultaneously.

    Parameters
    ----------
    sci_data : np.ndarray
        2D array of science data. Values should be in milliJanskys (mJy) per pixel.
    err_data : np.ndarray or None
        2D array of 1-sigma photometric uncertainties (standard deviation). 
    master_mask : np.ndarray
        2D boolean mask. A value of `True` indicates a bad/masked pixel.
    xcen : float
        The X pixel coordinate of the target's centroid.
    ycen : float
        The Y pixel coordinate of the target's centroid.
    r_ap_list : list of float or float
        A single radius or a list of aperture radii to evaluate simultaneously.
    ap_in_out : tuple of float
        A tuple of two floats: `(r_in, r_out)` for the background annulus.

    Returns
    -------
    pd.DataFrame
        A Pandas DataFrame with one row per aperture radius, containing explicit unit suffixes.
    """
    # 1. Unpack the annulus radii
    r_in, r_out = ap_in_out
    
    # Ensure r_ap_list is treated as a list even if a single float is passed
    if not isinstance(r_ap_list, (list, tuple, np.ndarray)):
        r_ap_list = [r_ap_list]
        
    columns = [
        'xcenter_pixel', 'ycenter_pixel', 
        'r_ap_pixel', 'r_in_pixel', 'r_out_pixel', 
        'aperture_area_pixel2', 'annulus_median_mjy_per_pix', 'bkg_std_mjy_per_pix', 
        'nsky_pixel2', 'nbadpix', 'aperture_sum_mjy', 'source_sum_mjy', 
        'source_sum_err_mjy', 'snr', 'abmag', 'abmag_err', 'badphot'
    ]
               
    # Catch invalid coordinates early
    if np.isnan(xcen) or np.isnan(ycen):
        return pd.DataFrame([[np.nan]*len(columns)] * len(r_ap_list), columns=columns)
        
    # Strip any hidden Astropy Units to prevent UnitConversionErrors
    sci_data = np.asarray(sci_data)
    if err_data is not None:
        err_data = np.asarray(err_data)
        
    positions = [(xcen, ycen)]
    
    try:
        # 2. Estimate Sky Background (COMPUTED ONLY ONCE)
        annulus = CircularAnnulus(positions, r_in=r_in, r_out=r_out)
        sigclip = SigmaClip(sigma=3.0, maxiters=5)
        sky_stats = ApertureStats(sci_data, annulus, mask=master_mask, sigma_clip=sigclip)
        
        msky = float(getattr(sky_stats.median[0], 'value', sky_stats.median[0]))
        ssky = float(getattr(sky_stats.std[0], 'value', sky_stats.std[0]))
        nsky = float(getattr(sky_stats.sum_aper_area[0], 'value', sky_stats.sum_aper_area[0]))
        
        # 3. Define all Apertures and Execute Photometry in One Pass
        apertures = [CircularAperture(positions, r=r) for r in r_ap_list]
        phot_table = aperture_photometry(sci_data, apertures, error=err_data, mask=master_mask)
        df_phot_raw = phot_table.to_pandas()
        
        # 4. Extract and Calculate Stats per Aperture
        records = []
        for i, r_ap in enumerate(r_ap_list):
            ap = apertures[i]
            
            # Exact fractional area and bad pixel count specifically for this aperture size
            ap_stats = ApertureStats(sci_data, ap, mask=master_mask)
            ap_area = float(getattr(ap_stats.sum_aper_area[0], 'value', ap_stats.sum_aper_area[0]))
            
            bad_stats = ApertureStats(master_mask.astype(float), ap)
            nbadpix = float(getattr(bad_stats.sum[0], 'value', bad_stats.sum[0]))
            
            # Photutils appends _0, _1, _2 to columns when processing a list of apertures
            suffix = f'_{i}'
            ap_sum = df_phot_raw.loc[0, f'aperture_sum{suffix}']
            
            if f'aperture_sum_err{suffix}' in df_phot_raw.columns:
                ap_sum_err_sq = df_phot_raw.loc[0, f'aperture_sum_err{suffix}']**2
            else:
                # Fallback if no error array is provided
                ap_sum_err_sq = np.abs(ap_sum - (ap_area * msky)) 
                
            # 5. Flux Correction
            source_sum = ap_sum - (ap_area * msky)
            
            # 6. DAOPHOT Errors
            sky_mean_err_term = (ap_area**2 * ssky**2) / nsky if nsky > 0 else 0.0
            source_sum_err = np.sqrt(ap_sum_err_sq + (ap_area * ssky**2) + sky_mean_err_term)
            snr = source_sum / source_sum_err if source_sum_err > 0 else np.nan
            
            # 7. AB Magnitude Conversion
            abmag, abmag_err = np.nan, np.nan
            if source_sum > 0:
                flux_jy = source_sum * 1e-3 
                abmag = -2.5 * np.log10(flux_jy) + 8.90
                abmag_err = (2.5 / np.log(10.0)) * (source_sum_err / source_sum)
                
            # 8. Flag Bad Photometry
            badphot = (ap_area == 0) or (source_sum <= 0) or (nsky < 10)
            
            records.append({
                'xcenter_pixel': xcen,
                'ycenter_pixel': ycen,
                'r_ap_pixel': r_ap,
                'r_in_pixel': r_in,
                'r_out_pixel': r_out,
                'aperture_area_pixel2': ap_area,
                'annulus_median_mjy_per_pix': msky,
                'bkg_std_mjy_per_pix': ssky,
                'nsky_pixel2': nsky,
                'nbadpix': nbadpix,
                'aperture_sum_mjy': ap_sum,
                'source_sum_mjy': source_sum,
                'source_sum_err_mjy': source_sum_err,
                'snr': snr,
                'abmag': abmag,
                'abmag_err': abmag_err,
                'badphot': badphot
            })
            
        return pd.DataFrame(records, columns=columns)
        
    except Exception as e:
        # Failsafe: Return a safely sized empty DataFrame if the math crashes
        return pd.DataFrame([[np.nan]*len(columns)] * len(r_ap_list), columns=columns)


def plot_coverage(data_summary: pd.DataFrame, gaia_subset_phase: dict, objdesig: str, gmag_limit: float, FIG_DIR: Path):
    num_phases = len(gaia_subset_phase)
    fig = plt.figure(figsize=(8 * num_phases, 8))
    gs = GridSpec(1, num_phases, figure=fig)

    for i, (phase, gaia_subset) in enumerate(gaia_subset_phase.items()):
        ax = fig.add_subplot(gs[i])
        
        # Plot Gaia Background
        ax.scatter(
            gaia_subset['ra'], gaia_subset['dec'], 
            s=10, alpha=0.5, color='green', 
            label=f'Gaia Source ($G<{gmag_limit}$)' if i == 0 else None
        )
        
        data_summary_gp_phase = data_summary[data_summary['phase'] == phase]

        # SPEED OPTIMIZATION: Now properly using itertuples()
        for idx, row in enumerate(data_summary_gp_phase.itertuples(index=False)):
            wcs = reconstruct_cutout_wcs(row)
            
            ny = int(getattr(row, 'cutout_size'))
            nx = ny
            ra_obj = getattr(row, 'ra')
            dec_obj = getattr(row, 'dec')

            # Plot Target Center
            ax.scatter(
                ra_obj, dec_obj, s=50, color='yellow', edgecolor='black', marker='*',
                label='Target' if idx == 0 else None
            )
            
            # Calculate Footprint Corners
            x_corners = [0, nx, nx,  0, 0]
            y_corners = [0,  0, ny, ny, 0]
            sky_corners = wcs.pixel_to_world(x_corners, y_corners)
            
            ra_poly = sky_corners.ra.deg
            dec_poly = sky_corners.dec.deg
            
            # GEOMETRY FIX: Prevent RA Wrap-Around Streaks
            if np.max(ra_poly) - np.min(ra_poly) > 180.0:
                ra_poly = np.where(ra_poly > 180.0, ra_poly - 360.0, ra_poly)
                
            label = 'SPHEREx Cutout FoV' if idx == 0 else None
            ax.plot(ra_poly, dec_poly, color='red', linewidth=1.5, alpha=0.2, label=label)
            
        # Formatting Time Safely using Pandas Datetime
        dates = pd.to_datetime(data_summary_gp_phase['DATE-OBS'])
        date_obs_min = dates.min().strftime('%Y-%m-%d')
        date_obs_max = dates.max().strftime('%Y-%m-%d')
        
        ax.set_title(f"Phase {phase}\n({date_obs_min} to {date_obs_max})")
        ax.grid(True, linestyle='--', alpha=0.3)
        ax.set_xlabel(r'$\alpha~(\degree)$')
        ax.set_ylabel(r'$\delta~(\degree)$')
        
        if i == 0:
            ax.legend()

    plt.tight_layout() 
    plt.suptitle(f"{objdesig}", y=1.02) 
    plt.savefig(FIG_DIR / f'coverage_{objdesig}.png', bbox_inches='tight')
    # plt.show()
    plt.close(fig) # Prevent memory leaks
    
def plot_spec(phot_summary: pd.DataFrame, objdesig: str, phase: int, FIG_DIR: Path):
    fig = plt.figure(figsize=(15, 6))
    ax = fig.add_subplot()

    phot_summary = phot_summary.sort_values(by='wl').reset_index(drop=True)
    x = phot_summary["wl"]
    apsum_mjy = phot_summary["source_sum_mjy"]
    apsum_err_mjy = phot_summary["source_sum_err_mjy"]

    sc = ax.scatter(x, apsum_mjy, c=phot_summary['r_hel'], cmap='RdBu', marker='o', zorder=5)

    ax.errorbar(x, apsum_mjy, yerr=apsum_err_mjy, fmt='none', ecolor='gray', alpha=0.6, zorder=4)
    ax.plot(x, apsum_mjy, lw=0.5, marker="none", alpha=0.7)

    # Clearly labeled Volatile Bands
    ax.axvspan(2.6, 2.8, color='skyblue', alpha=0.3, label=r'H$_2$O (2.7 $\mu$m)') 
    ax.axvspan(4.2, 4.4, color='skyblue', alpha=0.3, label=r'CO$_2$ (4.3 $\mu$m)') 
    ax.axvspan(4.6, 4.8, color='skyblue', alpha=0.3, label=r'CO (4.7 $\mu$m)')

    cbar = plt.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label(r"$R_{HEL}$ (au)")
    
    ax.set_xlabel(r"Wavelength ($\mu$m)")
    ax.set_ylabel("Aperture Flux (mJy)")
    
    # JD safely converted
    date_obs_min = Time(phot_summary['jd_utc'].min(), format='jd').to_datetime().strftime('%Y-%m-%d')
    date_obs_max = Time(phot_summary['jd_utc'].max(), format='jd').to_datetime().strftime('%Y-%m-%d')
    ax.set_title(f"{objdesig} | Phase {phase} ({date_obs_min} to {date_obs_max})")

    # Safe extraction of aperture radius
    r_ap_val = phot_summary['r_ap_km'].iloc[0] if 'r_ap_km' in phot_summary.columns else "N/A"
    
    ax.annotate(f"{len(phot_summary)} data points\n$\\rho = {r_ap_val}$ km",
                xy=(0.05, 0.95), xycoords='axes fraction',
                ha='left', va='top', 
                bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="gray", alpha=0.9))

    ax.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(FIG_DIR / f"spec_{objdesig}_phase{phase}.png", bbox_inches='tight')
    # plt.show()
    plt.close(fig)


def plot_cutout(phot_summary_filtered: pd.DataFrame, gaia_all: Any, objdesig: str, FITS_DIR: Path, FIG_DIR: Path):
    num_targets = len(phot_summary_filtered)
    if num_targets == 0:
        print("Warning: phot_summary_filtered is empty. Skipping plot.")
        return
        
    ncols = 5
    nrows = math.ceil(num_targets / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
    axes = np.atleast_1d(axes).flatten()

    for idx, row in enumerate(phot_summary_filtered.itertuples(index=False)):
        ax = axes[idx]
        
        filename = getattr(row, 'filename')
        xcen = getattr(row, 'xcen')
        ycen = getattr(row, 'ycen')
        phase = getattr(row, 'phase')
        wl = getattr(row, 'wl')
        r_hel = getattr(row, 'r_hel')
        vmag = getattr(row, 'vmag')
        r_ap_pixel = getattr(row, 'r_ap_pixel') 
        r_out_pixel = getattr(row, 'r_out_pixel') 
        
        fpath_fits = FITS_DIR / filename
        
        # FAILSAFE: If a FITS file failed to extract, mark it and skip cleanly
        if not fpath_fits.exists():
            ax.set_title(f"File Missing:\n{filename}", fontsize=8, color='red')
            ax.axis('off')
            continue

        with fits.open(fpath_fits) as hdul:
            sci_data = hdul[1].data.astype(np.float32)
            flag_data = hdul[3].data.astype(np.uint32)
            
            mask_flag = flag_to_mask(flag_data, flag_number=[2, 6, 7, 9, 10, 11, 12, 14, 15, 17, 22, 24, 26, 27, 28, 29])
            mask_source = gaia_to_mask(gaia_all, hdul[1], gmag_limit=18.0)
            
            master_mask = mask_flag | mask_source | np.isnan(sci_data)
            sci_masked = sci_data.copy()
            sci_masked[master_mask] = np.nan

        x_int, y_int = int(xcen), int(ycen)
        x_min, x_max = max(0, x_int - int(5*r_ap_pixel)), min(sci_masked.shape[1], x_int + int(5*r_ap_pixel))
        y_min, y_max = max(0, y_int - int(5*r_ap_pixel)), min(sci_masked.shape[0], y_int + int(5*r_ap_pixel))
        
        sci_data_cutout = sci_masked[y_min:y_max, x_min:x_max]
        
        try:
            vmin, vmax = ZScaleInterval().get_limits(sci_data_cutout) 
        except (ValueError, IndexError):
            vmin, vmax = 0, 1

        ax.imshow(
            sci_data_cutout, 
            origin='lower', 
            cmap='magma', 
            vmin=vmin, 
            vmax=vmax,
            extent=[x_min, x_max, y_min, y_max] 
        )

        # Plot Aperture Geometry
        circle_ap = plt.Circle((xcen, ycen), r_ap_pixel, color='red', fill=False, linestyle='-', linewidth=1.5, alpha=1)
        ax.add_patch(circle_ap)

        circle_out = plt.Circle((xcen, ycen), r_out_pixel, color='cyan', fill=False, linestyle='-', linewidth=1.5, alpha=1)
        ax.add_patch(circle_out)
        
        ax.set_title(f"{wl:.2f} $\mu$m | $r_{{hel}}$: {r_hel:.2f} au | V: {vmag:.2f} | r_ap: {r_ap_pixel:.2f} pix", fontsize=10)

    # Cleanup Empty Subplots
    for i in range(num_targets, len(axes)):
        axes[i].axis('off')

    plt.tight_layout()
    
    # Dynamically grab the phase from the first row for the filename, otherwise fallback to 'all'
    phase_val = getattr(phot_summary_filtered.iloc[0], 'phase') if hasattr(phot_summary_filtered, 'iloc') else 'all'
    
    plt.savefig(FIG_DIR / f"cutouts_{objdesig}_phase{phase_val}.png", bbox_inches='tight')
    # plt.close()
    plt.close(fig)

def combine_fits(
    df_phase: pd.DataFrame, 
    gaia_subset: Any, 
    objdesig: str, 
    phase: int, 
    chosen_ap_pixel: float, 
    FITS_DIR: Path, 
    COMBFITS_DIR: Path, 
    FIG_DIR: Path,
    bands: Optional[Dict[str, tuple]] = None,
    gmag_limit: float = 18.0
):
    """
    Extracts 2D cutouts for specified wavelength bands, masks background sources and flags, 
    computes a median-combined image, saves the result to a FITS file, and plots a summary figure.

    Parameters
    ----------
    df_phase : pd.DataFrame
        The photometry summary DataFrame filtered to a single phase and a single aperture size.
    gaia_subset : numpy.ndarray or pd.DataFrame
        The catalog of Gaia sources strictly filtered for the spatial footprint of this phase.
    objdesig : str
        The target's designation (e.g., '2P').
    phase : int
        The integer identifying the current observation phase.
    chosen_ap_pixel : float
        The base aperture radius in pixels used to determine the cutout size (5 * chosen_ap_pixel).
    FITS_DIR : pathlib.Path
        Directory containing the original individual FITS cutouts.
    COMBFITS_DIR : pathlib.Path
        Directory where the newly combined median FITS files will be saved.
    FIG_DIR : pathlib.Path
        Directory where the 1x5 band summary plots will be saved.
    bands : dict, optional
        A dictionary mapping band names to (wl_min, wl_max) tuples. 
        Defaults to predefined volatile bands if None.
    gmag_limit : float, optional
        The Gaia G-band magnitude limit used for star masking. Defaults to 18.0.
    """
    
    if bands is None:
        bands = {
            'dust_cont_1': (1.3, 2.0),
            'h2o': (2.6, 2.8),
            'dust_cont_2': (3.0, 4.0),
            'co2': (4.2, 4.4),
            'co': (4.6, 4.8)
        }

    # Calculate the universal 5*r_ap cutout radius
    cutout_radius = int(np.ceil(5 * chosen_ap_pixel))

    def get_padded_cutout(data_array, x_center, y_center, radius):
        """Extracts a (2*radius+1) square cutout, padding with NaNs if hitting detector edges."""
        h, w = data_array.shape
        out = np.full((2*radius+1, 2*radius+1), np.nan, dtype=np.float32)
        
        x_int, y_int = int(round(x_center)), int(round(y_center))
        x_min, x_max = x_int - radius, x_int + radius + 1
        y_min, y_max = y_int - radius, y_int + radius + 1
        
        valid_x_min, valid_x_max = max(0, x_min), min(w, x_max)
        valid_y_min, valid_y_max = max(0, y_min), min(h, y_max)
        
        out_x_min = valid_x_min - x_min
        out_x_max = out_x_min + (valid_x_max - valid_x_min)
        out_y_min = valid_y_min - y_min
        out_y_max = out_y_min + (valid_y_max - valid_y_min)
        
        if valid_x_min < valid_x_max and valid_y_min < valid_y_max:
            out[out_y_min:out_y_max, out_x_min:out_x_max] = data_array[valid_y_min:valid_y_max, valid_x_min:valid_x_max]
            
        return out

    # Ensure output directory exists
    COMBFITS_DIR.mkdir(parents=True, exist_ok=True)

    # Setup the Plot
    fig, axes = plt.subplots(1, len(bands), figsize=(4 * len(bands), 4))
    axes = np.atleast_1d(axes).flatten()

    for idx, (band_name, (wl_min, wl_max)) in enumerate(bands.items()):
        ax = axes[idx]
        
        # Filter dataframe for the current band
        band_df = df_phase[(df_phase['wl'] >= wl_min) & (df_phase['wl'] <= wl_max)]
        
        if band_df.empty:
            ax.set_title(f"{band_name}\n({wl_min}-{wl_max} $\\mu$m)\nNo Data", fontsize=10)
            ax.axis('off')
            continue
            
        cutout_stack = []
        
        # Extract cutouts
        for row in band_df.itertuples(index=False):
            fpath_fits = FITS_DIR / getattr(row, 'filename')
            
            # Failsafe: Skip if file is missing
            if not fpath_fits.exists():
                continue
            
            with fits.open(fpath_fits) as hdul:
                sci_data = hdul[1].data.astype(np.float32)
                flag_data = hdul[3].data.astype(np.uint32)
                
                mask_flag = flag_to_mask(flag_data, flag_number=[2, 6, 7, 9, 10, 11, 12, 14, 15, 17, 22, 24, 26, 27, 28, 29])
                # Call the pre-defined masking function
                mask_source = gaia_to_mask(gaia_subset, hdul[1], gmag_limit=gmag_limit)
                
                master_mask = mask_flag | mask_source | np.isnan(sci_data)
                sci_masked = sci_data.copy()
                sci_masked[master_mask] = np.nan
                
                cutout = get_padded_cutout(sci_masked, getattr(row, 'xcen'), getattr(row, 'ycen'), cutout_radius)
                cutout_stack.append(cutout)

        # Failsafe: If no valid cutouts were extracted
        if len(cutout_stack) == 0:
            ax.set_title(f"{band_name}\n({wl_min}-{wl_max} $\\mu$m)\nAll Files Missing", fontsize=10)
            ax.axis('off')
            continue

        # Perform the NaN-Median Stack
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning) 
            combined_image = np.nanmedian(np.array(cutout_stack), axis=0)

        # Save the Combined FITS
        hdu = fits.PrimaryHDU(combined_image)
        hdr = hdu.header
        hdr['OBJDESIG'] = (objdesig, 'Target designation')
        hdr['PHASE']    = (phase, 'Observation phase')
        hdr['BAND']     = (band_name, 'Custom SPHEREx spectral band')
        hdr['WL_MIN']   = (wl_min, 'Minimum wavelength (um)')
        hdr['WL_MAX']   = (wl_max, 'Maximum wavelength (um)')
        hdr['NCOMBINE'] = (len(cutout_stack), 'Number of exposures combined')
        hdr['BUNIT']    = ('mJy/pixel', 'Unit of array data')
        hdr['RAD_PIX']  = (cutout_radius, 'Cutout radius in pixels')
        
        out_filename = COMBFITS_DIR / f"{objdesig}_phase{phase}_{band_name}.fits"
        hdu.writeto(out_filename, overwrite=True)

        # Plot the Combined Image
        try:
            vmin, vmax = ZScaleInterval().get_limits(combined_image)
        except (ValueError, IndexError):
            vmin, vmax = 0, 1
            
        extent = [-cutout_radius, cutout_radius, -cutout_radius, cutout_radius]
        
        ax.imshow(combined_image, origin='lower', cmap='magma', vmin=vmin, vmax=vmax, extent=extent)
        
        # Draw the target aperture geometry
        circle_ap = plt.Circle((0, 0), chosen_ap_pixel, color='red', fill=False, linestyle='-', linewidth=1.5, alpha=1)
        ax.add_patch(circle_ap)
        
        ax.set_title(f"{band_name}\n({wl_min}-{wl_max} $\\mu$m) | N={len(cutout_stack)}", fontsize=11)
        ax.set_xlabel('$\Delta$ X [pix]')
        if idx == 0:
            ax.set_ylabel('$\Delta$ Y [pix]')

    # CRITICAL FIX: Plot formatting placed OUTSIDE the for loop
    plt.tight_layout()
    plt.savefig(FIG_DIR / f"combined_{objdesig}_phase{phase}.png", bbox_inches='tight')
    plt.close(fig)

if __name__ == "__main__":
    start_time = time.time()
    print(f"========== SPHEREx APPHOT PIPELINE START ==========")
    
    # ==========================================
    # 1. CONFIGURATION & PARAMETERS
    # ==========================================
    FITS_DIR     = directory.FITS_DIR
    DB_DIR       = directory.DB_DIR
    FIG_DIR      = directory.FIG_DIR
    REFCAT_DIR   = directory.REFCAT_DIR
    OUT_DIR      = directory.APPHOT_DIR
    COMBFITS_DIR = directory.COMBFITS_DIR
    
    # get objdesig from command line arguments
    parser = argparse.ArgumentParser(description="Run photometry pipeline for a specific target.")
    parser.add_argument('--objdesig', type=str, required=True, help='Target designation (e.g., 2P)')
    args = parser.parse_args()
    objdesig = args.objdesig
    
    fpath_db = DB_DIR / "db_filtered.parq"
    fpath_gaia = REFCAT_DIR / "gaiadr3_all.npy"
    
    # Target Parameters
    del_ra_deg, del_dec_deg = 0.2, 0.2
    gmag_limit = 18.0
    aperture_radii_km = [10000, 15000, 20000, 25000, 30000, 40000, 60000, 80000, 100000]

    print(f"[1/7] Configuration loaded for target: {objdesig}")

    # ==========================================
    # 2. DATA LOADING & PHASE GROUPING
    # ==========================================
    print(f"[2/7] Loading database and grouping phases...")
    # Use pushdown filtering to only load memory for obj
    data_summary = pd.read_parquet(
        fpath_db, 
        filters=[("objdesig", "==", objdesig)]
    )
    
    objdesig = objdesig.replace(" ", "")  # Remove spaces for file naming
    
    # Phase grouping based on observation date
    data_summary = spherex_phase_grouping(data_summary, threshold_days=28, key_dateobs='DATE-OBS')
    
    num_phases = data_summary['phase'].nunique()
    print(f"  -> Extracted {len(data_summary)} total observations.")
    print(f"  -> Grouped into {num_phases} distinct observation phases.")

    # ==========================================
    # 3. GAIA SUBSETTING & MEMORY MANAGEMENT
    # ==========================================
    print(f"[3/7] Filtering Gaia DR3 Catalog (Limit G < {gmag_limit})...")
    # 1. Load the memory map
    gaia_all_mmap = np.load(fpath_gaia, mmap_mode='r')
    
    # 2. PERFORM A SINGLE GLOBAL FILTER (The Speed Hack)
    # We evaluate the magnitude limit directly from the disk once.
    print("  -> Loading bright stars into active RAM...")
    valid_mag_mask = gaia_all_mmap['phot_g_mean_mag'] < gmag_limit
    
    # Extract only the passing stars into a fast, in-memory numpy array.
    # This drops the 11GB array down to a few hundred megabytes.
    gaia_light_ram = gaia_all_mmap[valid_mag_mask] 
    
    # We no longer need the 11GB disk reference for this run
    del gaia_all_mmap 

    gaia_subset_phase = dict()
    
    for phase in data_summary['phase'].unique():
        phase_mask = data_summary['phase'] == phase
        data_summary_phase = data_summary[phase_mask]
        
        # 3. Pass the lightweight RAM array to your subsetter
        gaia_subset = create_gaia_subset(
            gaia_light_ram, data_summary_phase, 
            del_ra_deg=del_ra_deg, del_dec_deg=del_dec_deg, gmag_limit=gmag_limit,
            ra_col='ra', dec_col='dec'
        )
        gaia_subset_phase[phase] = gaia_subset
        print(f"  -> Phase {phase}: {len(gaia_subset)} background stars isolated.")

    # Generate the Coverage Plot
    print("  -> Generating global coverage plot...")
    plot_coverage(data_summary, gaia_subset_phase, objdesig, gmag_limit, FIG_DIR)

    # ==========================================
    # 4. CROSS-MATCHING (OPTIMIZED)
    # ==========================================
    print("[4/7] Performing KD-Tree Cross-Matching with Gaia...")
    # Combine phase subsets to create a manageable catalog for KD-Tree cross-matching
    gaia_combined_subset = np.concatenate(list(gaia_subset_phase.values()))
    
    # Extract column names safely
    gaia_cols = gaia_combined_subset.dtype.names if not isinstance(gaia_combined_subset, pd.DataFrame) else gaia_combined_subset.columns
    gmag_col = 'phot_g_mean_mag' if 'phot_g_mean_mag' in gaia_cols else 'gmag'

    gaia_coords = SkyCoord(ra=gaia_combined_subset['ra']*u.deg, dec=gaia_combined_subset['dec']*u.deg)
    target_coords = SkyCoord(ra=data_summary['ra'].values*u.deg, dec=data_summary['dec'].values*u.deg)

    idx, sep2d, _ = target_coords.match_to_catalog_sky(gaia_coords)
    matched_gaia_sources = gaia_combined_subset[idx]

    # Convert angular separation to pixel distance
    data_summary['neargaia_gmag'] = matched_gaia_sources[gmag_col]
    data_summary['neargaia_dist_pixel'] = sep2d.arcsec / data_summary['pix_scale'].values
    print("  -> Cross-matching complete. Separations converted to pixels.")

    # ==========================================
    # 5. DYNAMIC APERTURE PARAMETERIZATION
    # ==========================================
    print("[5/7] Calculating dynamic aperture pixel radii...")
    pixel_scale_km = (data_summary["pix_scale"] * ((1*u.arcsec).to(u.rad).value) * data_summary["r_obs"] * (1*u.au).to(u.km).value)
    
    for i, r_km in enumerate(aperture_radii_km, start=1):
        data_summary[f"r_ap_{i:02d}_km"] = r_km
        data_summary[f"r_ap_{i:02d}_pix"] = r_km / pixel_scale_km

    # Set background annulus size based on the maximum aperture
    max_ap_col = f"r_ap_{len(aperture_radii_km):02d}_pix"
    data_summary["r_in_pix" ] = 1.5 * data_summary[max_ap_col]
    data_summary["r_out_pix"] = 3.0 * data_summary[max_ap_col]

    # ==========================================
    # 6. PHOTOMETRY BATCH PROCESSING
    # ==========================================
    phot_results = []
    
    list_columns = [
        'filename', 'objdesig', 'obsid', 'phase', 'detector', 'wl', 'wlwidth', 'sun_jy', 
        'xcen', 'ycen', 'ltv1', 'ltv2', 'cutout_size', 'pix_scale', 
        'ra', 'dec', 'r_hel', 'r_obs', 'alpha', 'hel_ecl_lon', 
        'hel_ecl_lat', 'obs_ecl_lon', 'obs_ecl_lat', 'racosdec_rate', 
        'dec_rate', 'sky_motion', 'sky_motion_pa', 'vmag', 'jd_utc', 'jd_tdb',
        'neargaia_gmag', 'neargaia_dist_pixel'
    ]

    total_files = len(data_summary)
    print(f"[6/7] Beginning aperture photometry extraction for {total_files} cutouts...")

    for row_idx, row in enumerate(data_summary.itertuples(index=False), start=1):
        fpath_fits = FITS_DIR / getattr(row, 'filename')
        phase = getattr(row, 'phase')
        
        with fits.open(fpath_fits) as hdul:
            sci_data = hdul[1].data.astype(np.float32)
            err_data = np.sqrt(hdul[2].data.astype(np.float32))
            
            flag_data = hdul[3].data.astype(np.uint32)
            mask_flag = flag_to_mask(flag_data, flag_number=[2, 6, 7, 9, 10, 11, 12, 14, 15, 17, 22, 24, 26, 27, 28, 29])
            
            mask_source = gaia_to_mask(gaia_subset_phase[phase], hdul[1], gmag_limit=gmag_limit, radius_scale=0.6)
            master_mask = mask_flag | mask_source | np.isnan(sci_data)

        # Build dynamic list of pixel radii for this row
        r_ap_pix_list = [getattr(row, f"r_ap_{i:02d}_pix") for i in range(1, len(aperture_radii_km) + 1)]

        phot = perform_aperture_photometry(
            sci_data=sci_data,
            err_data=err_data,
            master_mask=master_mask,
            xcen=getattr(row, 'xcen'),
            ycen=getattr(row, 'ycen'),
            r_ap_list=r_ap_pix_list,
            ap_in_out=(getattr(row, 'r_in_pix'), getattr(row, 'r_out_pix'))
        )
        
        # Inject metadata into the result dataframe
        phot['r_ap_km'] = aperture_radii_km
        for col in list_columns:
            phot[col] = getattr(row, col)

        phot_results.append(phot)
        
        # Status Tracker
        if row_idx % 500 == 0 or row_idx == total_files:
            print(f"  -> Processed {row_idx}/{total_files} cutouts...")

    print("  -> Photometry complete. Compiling results...")
    # Concatenate everything exactly once
    phot_summary = pd.concat(phot_results, ignore_index=True)
    out_file = OUT_DIR / f"{objdesig}.csv"
    phot_summary.to_csv(out_file, index=False)
    print(f"  -> Saved master summary to: {out_file}")

    # ==========================================
    # 7. PHASE-GROUPED PLOTTING & STACKING
    # ==========================================
    print(f"[7/7] Generating Phase Summary Plots and Band Stacks...")
    for phase in phot_summary['phase'].unique():
        # Isolate the base data for this phase
        df_phase = phot_summary[phot_summary['phase'] == phase]
        
        # A. Dynamically Determine Consistent Aperture Size
        ap_min_pixels = df_phase.groupby('r_ap_km')['r_ap_pixel'].min()
        valid_ap_kms = ap_min_pixels[ap_min_pixels > 2.0].index.tolist()
        
        if valid_ap_kms:
            chosen_ap_km = min(valid_ap_kms)
        else:
            chosen_ap_km = df_phase['r_ap_km'].max()
            
        print(f"  -> Phase {phase}: Selected {chosen_ap_km} km aperture (min pixel radius > 2).")

        # B. Filter Dataframes for plotting/stacking
        mask_phot = (df_phase['badphot'] == False) & (df_phase['snr'] > 1)
        mask_ap = (df_phase['r_ap_km'] == chosen_ap_km)
        
        df_spec = df_phase[mask_phot & mask_ap].copy()
        df_cutout = df_phase[mask_phot & mask_ap].copy() # Cutouts don't strictly need the SNR filter
        
        if not df_spec.empty:
            print(f"  -> Phase {phase}: Plotting spectrum with {len(df_spec)} valid exposures...")
            plot_spec(df_spec, objdesig, phase, FIG_DIR)
        else:
            print(f"  -> Phase {phase}: [WARNING] No valid data points passed the SNR/badphot filter.")
        
        if not df_cutout.empty:
            print(f"  -> Phase {phase}: Plotting {len(df_cutout)} individual FITS cutouts...")
            plot_cutout(df_cutout, gaia_subset_phase[phase], objdesig, FITS_DIR, FIG_DIR)
            
            # C. Combine FITS into Median Stacks
            print(f"  -> Phase {phase}: Stacking FITS by volatile bands...")
            # Calculate a representative pixel radius for the cutout bounds
            median_ap_pixel = df_cutout['r_ap_pixel'].median()
            
            combine_fits(
                df_phase=df_cutout, 
                gaia_subset=gaia_subset_phase[phase], 
                objdesig=objdesig, 
                phase=phase, 
                chosen_ap_pixel=median_ap_pixel, 
                FITS_DIR=FITS_DIR, 
                COMBFITS_DIR=COMBFITS_DIR, 
                FIG_DIR=FIG_DIR,
                gmag_limit=gmag_limit
            )
            
    # Calculate execution time
    elapsed = time.time() - start_time
    minutes, seconds = divmod(int(elapsed), 60)
    print(f"========== PIPELINE EXECUTION COMPLETE ==========")
    print(f"Total time elapsed: {minutes} minutes, {seconds} seconds.")
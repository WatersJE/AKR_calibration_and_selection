#!usr/bin/env python3
"""
11/02/2021

Produce calibrated AKR fluxes from a daily Wind/Waves L2 file


"""
import os
import argparse
import numpy as np
import pandas as pd

import read_wind_waves_l2
import waves_rad1_l2_analysis as l2
import background
import snr
import calibration as cal


def default_freqs():

    freqs = np.loadtxt('default_freq_list.txt')

    return freqs


def process_l2(date, data_dir):
    """
    Read L2 data file and concatenate observations into dataframe for whole day

    """

    if isinstance(date, pd.Timestamp):

        date = date.strftime('%Y%m%d')

    l2_fp = os.path.join(data_dir, filename.format(date))

    l2_obj, n_sw = read_wind_waves_l2.read_l2_hres(l2_fp, False)

    if isinstance(date, pd.Timestamp):
        date_str = date.strftime('%Y%m%d')
    elif isinstance(date, str):
        date_str = date

    # CHECK VALIDATION
    # create and store empty dataframe in event that all sweeps are None
    # (NotImplementedError raised in this case - use try except)
    try:
        l2_df = l2.concat_l2_sweeps(l2_obj, n_sw,
            log_filename=os.path.join(
                data_dir, 'L2_validation_{}.log'.format(date_str)))
    except NotImplementedError:
        # create empty data
        return None

    # calculate spin normalised standard deviation of spin-axis antenna
    orig_power = l2_df.loc[:, 'AMPL_Z']

    metric = lambda p_z: np.std(p_z) / np.mean(p_z)

    sigma_z = orig_power.groupby(orig_power.index // 8).apply(metric)

    # calculate background
    bg_df = background.new_background_file(l2_df)

    snr_avgs = np.zeros(l2_df.shape[0] // 8)

    # calculate SNR
    for i, spin_df in l2_df.groupby(l2_df.index // 8):

        # assign info vars
        power = spin_df['AMPL_Z']
        freq = spin_df['FREQ'].unique()[0]
        
        snr_avgs[i] = snr.snr_average(bg_df, power, freq, 'Z')

    # 9th March 2021
    # now subtract background
    for i, f in enumerate(
        np.sort(l2_df['FREQ'].dropna().unique())):

        subdf = l2_df.loc[l2_df['FREQ'] == f, :]

        subdf['AMPL_Z'] = background.subtract_background(
            bg_df, subdf['AMPL_Z'], f, 'Z'
        )

        l2_df.loc[l2_df['FREQ'] == f, 'AMPL_Z'] = subdf['AMPL_Z']

    # Select cols to retain - make dataframe of spin averages and sigma_z
    keep_cols = ['FREQ', 'AMPL_Z', 'DATETIME_Z', 'SWEEP', 'sweep_flag']

    spin_df = l2_df.loc[l2_df.index[::8], keep_cols]

    spin_df = spin_df.reset_index()

    spin_df['SNR_dB'] = snr_avgs

    # Remove RFI-contaminated channel (52 kHz)
        # (See later - for the time being just set 52 kHz 
        # to NaN at end)
    # spin_df = remove_rfi(spin_df, 52)

    z_power = l2_df.loc[:, 'AMPL_Z']

    # average over spins
    ampl_z = z_power.groupby(z_power.index // 8).apply(np.mean)

    spin_df['AMPL_Z'] = ampl_z

    spin_df['sigma_z'] = sigma_z

    return spin_df


def normalise_flux(flux, rad_dist, target='Earth'):
    """
    Normalise flux densities to 1 AU, first correcting by the radial distance
    from Earth

    If target is 'Sun', `rad_dist` is corrected appropriately to give 
    distance from Sun in Solar radii, assuming spacecraft is on nightside and close
    to GSE X

    Parameters
    ----------
    flux : float, np.ndarray, pd.Series
        unnormalised flux density value/s

    rad_dist : float, np.ndarray, pd.Series
        corresponding radial distance from Earth (in R_E)

    Returns
    -------
    flux_norm : float, np.ndarray, pd.Series
        normalised flux density values
    """

    assert target in ['Earth', 'Sun']

    
    earth_radius = 6378100.0 #m
    # earth_radius = astro_const.R_earth.value
    
    # au = astro_const.au.value
    au = 149597870700.0 # in m

    if target is 'Earth':
        dist_in_au = (rad_dist * earth_radius) / au
    
    elif target is 'Sun':
        dist_in_au = (au + (rad_dist * earth_radius)) / au
    
    
    scaling = np.power(dist_in_au, 2.)

    flux_norm = flux * scaling

    return flux_norm


def z_akr_flux(z_power, latitude, gain):
    """
    Given the power spectral density received by the Z antenna and the
    spacecraft GSE latitude, return the associated flux, assuming a point
    radio source located at Earth's center and with zero linear polarisation
    (following simplifications to MF80)

    Assumes z_power given in V^2Hz^-1
    """

    return np.divide(2. * z_power, gain * np.power(np.cos(latitude), 2))


def select_akr(df, select_label, out_label,
    threshold=0.1):

    df[out_label] = np.where(df['sigma_z'] >= threshold, df[select_label], np.nan)

    return df


def generate_empty_dataframe():
    n_spectra = 475

    freqs = default_freqs()

    n_freqs = freqs.shape[0]

    dt_add = pd.Timedelta(23, 'H') + \
        pd.Timedelta(59, 'min') + \
        pd.Timedelta(30, 's')

    dtimes = pd.date_range(date, date+dt_add, periods=n_spectra)

    spectrum_i = np.repeat(np.arange(n_spectra), n_freqs)

    dtimes = np.repeat(dtimes, n_freqs)

    freqs = np.array(np.tile(freqs, n_spectra), dtype=np.float64)

    df_dict = {
            'freq': freqs,
            'datetime_ut': dtimes,
            'spectrum_i': spectrum_i,
            'akr_flux_si_1au': np.repeat(np.nan, freqs.shape[0]),
            'snr_db': np.repeat(np.nan, freqs.shape[0])}

    df = pd.DataFrame(df_dict)

    return df


def main(date, version='01',
    save_calibrated_data=False):

    datadir = '../data/wind_waves_l2_hres/1999'
    
    # calls waves_rad1_l2_analysis (externally)
    spin_df = process_l2(date, datadir)

    if spin_df is None:
        # all data invalidated - generate empty dataframe

        out_df = generate_empty_dataframe()

    else:
        
        flux_df = cal.calibration(spin_df)
        print(flux_df.columns)
        
        flux_label = 'flux_si'

        if save_calibrated_data:
        
            # 3 minute resolution dataframe with calibrated flux
            out_flux_df = cal.create_calibrated_flux_dataframe(flux_df, flux_label)
            
            calibrated_data_file = '../data/calibrated_data/final/1999'

            calibrated_flux_fp = 'wi_wa_rad1_l3_no_selection_{}.csv'.format(
                date.strftime('%Y%m%d'),
                date.strftime('%j')
            )

            fp = os.path.join(calibrated_data_file, calibrated_flux_fp)

            out_flux_df.to_csv(fp, index=False)

        # if normalise_wind_flux:

        # normalise data (Wind to 1 AU for comparison)
        flux_df['flux_si_1au'] = normalise_flux(flux_df['flux_si'], flux_df['RADIUS'])

        flux_label = 'flux_si_1au'

        out_df = cal.create_calibrated_flux_dataframe(flux_df, flux_label, flux_label)

        # now selection on sigma_z        
        out_df = select_akr(out_df, flux_label, 'akr_flux_si_1au')

        print(out_df.columns)

    # Save output
    out_dir = '../data/masked_flux/1999'
    
    out_fn = 'wi_wa_rad1_l3_akr_{}_v{}.csv'.format(date.strftime('%Y%m%d'),
        version)

    keep_cols = ['datetime_ut', 'freq', 'snr_db', 'akr_flux_si_1au']
    out_df = out_df.loc[:, keep_cols]

    fp = os.path.join(out_dir, out_fn)
    out_df.to_csv(fp, na_rep='-1', index=False)

    return None

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('date')

    parser.add_argument('version')

    args = parser.parse_args()

    date = pd.Timestamp(args.date)

    version = args.version

    print(date.strftime('%d %b %Y - DOY %j'))

    main(date, version, save_calibrated_data=True)
    # for testing (eg background subtraction)
    # process_l2(date, '../data/wind_waves_l2_hres')
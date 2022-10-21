#!usr/bin/env python3
"""
18/02/2021

Program to contain all relevant Wind/Waves calibration
"""
import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.io.idl import readsav
import scipy.interpolate as interpolate


def add_position(source_df, source_dt_label,
    target_df, target_dt_label, source_data=['RADIUS']):

    def unix_time(df, index=False,
        dt_label='datetime_ut'):
        'For converting datetimes into Unix time for interpolation and eg adding position data'
        new_df = df
        if not index:
            new_df['unix'] = new_df[dt_label].astype(np.int64) / 1e9
        else:
            new_df['unix'] = new_df.index.astype(np.int64) / 1e9

        return new_df

    if 'unix' not in source_df.columns:

        source_df = unix_time(source_df, dt_label=source_dt_label)

    target_df = unix_time(target_df, dt_label=target_dt_label)

    print_minmax = lambda df: print(df['unix'].min(), df['unix'].max())
    print_minmax(source_df)
    print_minmax(target_df)

    for label in source_data:

        func = interpolate.interp1d(source_df['unix'], source_df[label])

        target_df[label] = func(target_df['unix'])

    return target_df



def gain_factor(gamma, length):
    'Calculate G, given Gamma ratio and antenna length in m'

    z0 = 376.730313668

    gain = (gamma**2) * (length**2) * z0

    return gain


def z_flux(z_power, gain):
    """
    Calibrate Z antenna power, omitting term due to Wind latitude
    (ie for Type IIIs - still assumption of point source implicit
    here)
    """

    return np.divide(2. * z_power, gain)


def z_akr_flux(z_power, latitude, gain):
    """
    Given the power spectral density received by the Z antenna and the
    spacecraft GSE latitude, return the associated flux, assuming a point
    radio source located at Earth's center and with zero linear polarisation
    (following simplifications to MF80)

    Assumes z_power given in V^2Hz^-1
    """

    return np.divide(2. * z_power, gain * np.power(np.cos(latitude), 2))


def calibrate_z(power, latitude):
    """
    Perform calibration on the Z antenna measurements using equation 3
    of Waters+2021

    Takes:
        Z power in V^2.Hz^-1
        GSE latitude in radians

    Should then give flux in W.m^2.Hz^-1 as per equation 1 of Zarka 2004
    """

    gamma = 0.4
    length_z = 4.65 # m - quoted in Wind Docs (paragraph before §1.4.1.2)

    # new value using Z antenna instrument values and Zarka+2004
    gain = gain_factor(gamma, length_z)

    flux = z_akr_flux(power, latitude, gain)

    return flux


def calibration(spin_df, dt_label='DATETIME_Z', year=None):

    if year is None:
        eph_data = '/data/jw27g13/ephemeris/wind_ephemeris_45s_2003.csv'
    else:
        
        eph_data = '../../../substorm_akr_statistical_study/data/ephemeris/wind_ephemeris_45s_{}_magnetic_viewing.csv'.format(year)
        
    eph_df = pd.read_csv(eph_data, parse_dates=['Epoch'])

    spin_df = add_position(eph_df, 'Epoch',
        spin_df, 'DATETIME_Z', ['GSE_LAT', 'RADIUS'])

    power = spin_df['AMPL_Z'] * 1.e-12
    lat = np.radians(spin_df['GSE_LAT'])

    spin_df['flux_si'] = calibrate_z(power, lat)

    return spin_df


def contiguous_sweep_datetimes(l3_df, sweep_label, dt_label,
    sweep_period=183, return_original=False):
    """
    Return appropriate datetime labels for each 3-minute spectrum of L2 data,
    as well as for non-existent sweeps to maintain contiguity (as with
    spectrograms)

    taken from spectrogram_plotting.spectrogram_array
    """

    # Get list of dates of recorded sweep cycles
    dates = np.array([pd.Timestamp(np.min(df[dt_label].values))
        for _, df in l3_df.groupby(sweep_label)])
    
    # get ideal number of sweeps for the period in the data, given sweep_period
    # in seconds
    n_ideal_sw = np.ceil((dates[-1] - dates[0]).total_seconds() / sweep_period)

    dtimes = pd.date_range(dates[0], dates[-1], periods=n_ideal_sw)

    if not return_original:
        return dtimes
    else:
        return dates, dtimes


def ft_array(csv_df, param_label, set_zero_nan=True, zero=True):
    """
    Create 2D array corresponding to frequency-time spectrogram of `param_label`

    `sweep_period` in seconds
    """
    sweep_label = 'SWEEP'
    freq_label = 'FREQ'
    dt_label = 'DATETIME_Z'
    # wind sweep period in seconds
    sweep_period = 183

    sweep_dates, ideal_dates = contiguous_sweep_datetimes(csv_df, sweep_label, dt_label,
        return_original=True)
    # print(sweep_dates.shape)

    # Store indices for sweeps that correspond to datetime position in array
    sweep_array_inds = np.array([pd.Index(ideal_dates).get_loc(d, method='nearest')
        for d in sweep_dates])
    # print(sweep_array_inds.shape)

    freqs = np.sort(csv_df[freq_label].dropna().unique())
    
    # array for sampled frequency range and sweep cycles
    if zero:
        out_array = np.zeros((freqs.shape[0], int(ideal_dates.shape[0])))
    else:
        out_array = np.ones((freqs.shape[0], int(ideal_dates.shape[0])))
    # Filling array with sampled frequency means
    for arr_i, (_, sweep_df) in zip(sweep_array_inds,
        csv_df.groupby(sweep_label)):

        # print(arr_i, sweep_df['SWEEP'].unique())


        # groupby to average parameters
        parameter_mean = sweep_df.groupby(freq_label).agg({param_label: np.nanmean})

        # apparently no need to manually flip - first output shows
        # upside-down spectra
        means = parameter_mean[param_label].values
        # if np.all(sweep_df['sweep_flag'] == 1):
        #     print(means.shape)

        # assuming array sizes less than n_freqs have missing vales/erroneous
        # sweeps, so set to NaN
        if len(means) < freqs.shape[0]:
            means = np.repeat(np.nan, freqs.shape[0])
        
        # out_array[:, int(sw_i)] = flipped_means
        out_array[:, arr_i] = means

        # if np.all(sweep_df['sweep_flag'] == 1):
        #     print(sweep_df.loc[:, ['SWEEP', 'sigma_z']])
        #     print(sweep_df.shape)

    
    if set_zero_nan:
        # prevent division by 0 in colorbar
        out_array = np.where(np.isclose(out_array, 0., atol=1e-31), np.nan, out_array)
    # print(np.min(out_array))
    # print(np.nanmin(out_array))
    # print(np.max(out_array))
    # print(np.nanmax(out_array))

    return out_array



def create_calibrated_flux_dataframe(data,
    flux_label='flux_si', out_flux_label='flux_si', sweep_label='SWEEP'):
    """
    Mostly identical to create_selected... but for calibrated data only

    Takes 3 minute arrays of SNR and sigma_z also
    """
    # NEW 14TH MARCH 
    # (setting freqs only as those from non-flagged data made arrays
    # different lengths - error on out_df creation)
    # data = data.loc[data['sweep_flag'] == 0, :]

    flux = ft_array(data, flux_label)
    snr = ft_array(data, 'SNR_dB', set_zero_nan=False)
    sigma_z = ft_array(data, 'sigma_z', set_zero_nan=False)

    sweep_flag = ft_array(data, 'sweep_flag', set_zero_nan=False, zero=False)

    real_dtimes, out_dtimes = contiguous_sweep_datetimes(data,
        'SWEEP', 'DATETIME_Z', return_original=True)
    freqs = np.sort(data['FREQ'].unique())
    # NEW 11TH MARCH - PREVENT WEIRD FREQUENCIES?
    # freqs = np.sort(data.loc[data['sweep_flag'] == 0, 'FREQ'].unique())

    series_length = flux.shape[0] * flux.shape[1]
    print(series_length) 

    # tile frequencies for each spectrum
    n_spectra = flux.shape[1]

    # get appropriate list of frequencies for output dataframe
    df_freqs = np.array(np.tile(freqs, n_spectra), dtype=np.float64)
    print(df_freqs.shape[0])
    n_freqs = flux.shape[0]

    # flatten all 3-minute ft arrays
    flux = pd.Series(flux.flatten(order='F'))

    snr = pd.Series(snr.flatten(order='F'))

    sigma_z = pd.Series(sigma_z.flatten(order='F'))

    sweep_flag = pd.Series(sweep_flag.flatten(order='F'))

    # get indexes for each spectrum, assigned to data from each frequency
    spectrum_i = np.repeat(np.arange(n_spectra), n_freqs)
    # 15th March - so not to overlap with bad sweep nos from L2 data
    # spectrum_i = np.repeat(data['SWEEP'].unique(), n_freqs)
    # 15th March - using original, as no problem with overlap if don't concat bad_dfs
    # (Now using the empty spectra of ft_array output to infer bad sweeps...)

    # dtimes = np.repeat(dtimes, n_freqs) # old - `dtimes` was the list of "fake" dts
    # dtimes = np.repeat(real_dtimes, n_freqs) # 15 March - actual min sweep times of "good" sweeps
    dtimes = np.repeat(out_dtimes, n_freqs) # 15 March - fake times of all sweeps - including empty spectra created in ft_array calls

    # if data['sweep_flag'].unique().shape[0] > 1:
    #     raise ValueError("Multiple values for sweep flag")
    # else:
        # sweep_flag = data['sweep_flag'].unique()[0]
        

    df_dict = {
        'freq': df_freqs,
        'datetime_ut': dtimes,
        'spectrum_i': spectrum_i,
        'snr_db': snr,
        'sigma_z': sigma_z,
        out_flux_label: flux,
        'sweep_flag': sweep_flag}

    out_df = pd.DataFrame(df_dict)

    return out_df


def apply_crosscal_spectrum(cal_df):
    """
    16/4/21

    Load and apply the scaling spectrum derived from comparing Type III bursts
    etc

    """
    fp = 'cross_calibration_scaling_spectrum.csv'
    
    # load scaling spectrum, setting frequency as index
    df = pd.read_csv(fp, index_col=0)

    cal_df['new_flux'] = [flux / df.loc[f, 'avg_offset']
        for f, flux in zip(cal_df['freq'], cal_df['flux_si'])]
    
    cal_df = cal_df.drop(columns=['flux_si'])
    cal_df = cal_df.rename(columns={'new_flux': 'flux_si'})

    return cal_df

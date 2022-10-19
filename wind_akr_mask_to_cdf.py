#!usr/bin/env python3
"""
12th april

Modified from https://gitlab.obspm.fr/maser/nda-cdf/-/blob/master/NDA_RT1toCDF.py

Convert .CSV files of AKR-masked, calibrated Wind data to .CDF
"""
from pathlib import Path
import os
import numpy as np 
import pandas as pd
from datetime import datetime

from spacepy import pycdf


def load_csv(date, csv_dir):
    """
    Load CSV file from csv_dir / wi_wa_rad1_l3_akr_YYYYMMDD_v01.csv

    `date` : pd.Timestamp
    `csv_dir` : pathlib.Path, str
    """

    filename = 'wi_wa_rad1_l3_akr_{}_v01.csv'.format(
        date.strftime('%Y%m%d'))
    filepath = os.path.join(csv_dir, filename)

    df = pd.read_csv(filepath, parse_dates=['datetime_ut'],
        dtype={'akr_flux_si_1au': np.float64}, na_values='-1')

    flux_arr = ft_array(df, 'akr_flux_si_1au')
    snr_arr = ft_array(df, 'snr_db')
    dt_arr = np.sort(df['datetime_ut'].unique())
    freq_arr = np.sort(df['freq'].unique())

    return dt_arr, freq_arr, flux_arr, snr_arr


def ft_array(csv_df, param_label):
    """
    Create 2D array corresponding to frequency-time spectrogram of `param_label`

    """
    sweep_label = 'spectrum_i'
    freq_label = 'freq'
    dt_label = 'datetime_ut'
    # wind sweep period in seconds
    sweep_period = 183

    freqs = np.sort(csv_df[freq_label].dropna().unique())

    if sweep_label not in csv_df.columns:

        n_freqs = freqs.shape[0]
        n_spectra = csv_df[dt_label].unique().shape[0]
        csv_df[sweep_label] = np.repeat(np.arange(n_spectra), n_freqs)

    # Get list of dates of recorded sweep cycles
    dates = np.array([pd.Timestamp(np.min(df[dt_label].values))
        for _, df in csv_df.groupby(sweep_label)])

    # get ideal number of sweeps for the period in the data
    n_ideal_sw = np.ceil((dates[-1] - dates[0]).total_seconds() / sweep_period)

    ideal_dates = pd.date_range(dates[0], dates[-1], periods=n_ideal_sw)

    # Store indices for sweeps that correspond to datetime position in array
    sweep_array_inds = np.array([pd.Index(ideal_dates).get_loc(d, method='nearest')
        for d in dates])
    
    # array for sampled frequency range and sweep cycles
    out_array = np.zeros((freqs.shape[0], int(n_ideal_sw)))
    
    # Filling array with sampled frequency means
    for arr_i, (_, sweep_df) in zip(sweep_array_inds,
        csv_df.groupby(sweep_label)):

        # groupby to average parameters
        parameter_mean = sweep_df.groupby(freq_label).agg({param_label: 'mean'})

        means = parameter_mean[param_label].values

        # assuming array sizes less than n_freqs have missing vales/erroneous
        # sweeps, so set to NaN
        if len(means) < freqs.shape[0]:
            means = np.repeat(np.nan, freqs.shape[0])
        
        out_array[:, arr_i] = means

    # prevent division by 0 in colorbar
    out_array = np.where(np.isclose(out_array, 0., atol=1e-31), np.nan, out_array)
    
    return out_array


def set_istp_attr(cdf, cdf_path, cdf_stem, date, version='01'):
    # cdf.attrs['Logical_file_id'] = cdf_path.stem

    # SETTING ISTP GLOBAL ATTRIBUTES
    cdf.attrs['TITLE'] = 'Wind/Waves flux density collection calibrated for Auroral Kilometric Radiation'
    cdf.attrs['Project'] = [
        'PADC>Paris Astronomical Data Centre',
        'OBSPM>Observatoire de Paris'
    ]
    cdf.attrs['Discipline'] = 'Space Physics>Magnetospheric Science'
    cdf.attrs['Data_type'] = 'L3>Level3'
    cdf.attrs['Descriptor'] = 'AKR'
    cdf.attrs['Data_version'] = version
    cdf.attrs['Instrument_type'] = 'Radio and Plasma Waves (space)'
    # logical_file_id should be filename without '.cdf' (or preceding filepath - see Path.stem)
    cdf.attrs['Logical_file_id'] = cdf_stem
    cdf.attrs['Logical_source'] = 'wi_wa_rad1_l3_akr'
    cdf.attrs['Logical_source_description'] = 'Wind/Waves RAD1 Level 3 Calibrated AKR Flux Density'
    cdf.attrs['File_naming_convention'] = 'source_type_descriptor_yyyymmdd_ver'
    cdf.attrs['Mission_group'] = 'Wind'
    cdf.attrs['PI_name'] = 'K. Issautier'
    cdf.attrs['PI_affiliation'] = [
        (
            'LESIA>LESIA, Observatoire de Paris, PSL Research University, '
            'CNRS, Sorbonne Universites, Universite de Paris, 92195 Meudon, France'
        )
    ]
    cdf.attrs['Source_name'] = 'wi_wa_rad1>Wind/Waves RAD1'
    cdf.attrs['TEXT'] = (
        "Wind/WAVES is a set of orthogonal radio antennae on board the "
        "spin-stabilised Wind spacecraft. The RAD1 receiver operates between "
        "20-1040 kHz and makes 8 radio measurements on the axial Z antenna and the "
        "synthetically-inclined S antenna at a given frequency during a ~3 " 
        "second spin. Z antenna measurements are calibrated based on a "
        "modified GP inversion that makes assumptions relevant to the AKR source " 
        "(see paper for details). The flux densities are resampled and then "
        "selected based on the spin-normalised variability of the axial antenna"
        "to select AKR-associated emissions. Finally, selected "
        "data are normalised to an observer's distance of 1 AU."
    )
    cdf.attrs['Generated_by'] = [
        "SEP>Space Environment Physics Group, University of Southampton",
        "LESIA>Laboratoire d'Etudes Spatiales et d'Instrumentation en Astrophysique"
    ]
    cdf.attrs['Generation_date'] = '{}'.format(
        np.datetime64('now').astype(str) + 'Z'
    )
    cdf.attrs['LINK_TEXT'] = 'The data is available at'
    cdf.attrs['LINK_TITLE'] = 'MASER/PADC'
    cdf.attrs['HTTP_LINK'] = 'https://doi.org/10.25935/wxv0-vr90'
    cdf.attrs['DOI'] = 'https://doi.org/10.25935/wxv0-vr90'
    cdf.attrs['MODS'] = [
        (
            'v01: First release'
        )
    ]
    cdf.attrs['Parents'] = 'wi_wa_rad_l2_{}_v01.dat'.format(date.strftime('%Y%m%d'))
    cdf.attrs['Rules_of_use'] = [
        (
            'The data is distributed under licence CC-BY-4.0 '
            '(Creative Commons Attribution 4.0)'
        ),
        (
            "We kindly request that the authors of any publications that use "
            "this data contact the authors of this dataset, and "
            "include a citation to the reference below as well as appropriate "
            "acknowledgements in their work."
        ),
        'References:',
        (   # MASER landing page reference
            "J.E. Waters, B. Cecconi, X. Bonnin, & L. Lamy (2021). Wind/Waves "
            "flux density collection calibrated for Auroral Kilometric "
            "Radiation (Version 1.0) [Data set]. PADC. "
            "https://doi.org/10.25935/wxv0-vr90 , "
            "Waters J. E. et al. (2021) J. Geophys. Res. (Submitted)."
        ),
        'Acknowledgements: see the acknowledgement field.'
    ]
    cdf.attrs['Skeleton_version'] = version
    cdf.attrs['Software_version'] = version
    # cdf.attrs['Software_language'] = 'python3'
    cdf.attrs['Software_name'] = 'WindWaves_AKR_calibration_selection'
    cdf.attrs['Time_resolution'] = '3 minutes'
    cdf.attrs['Acknowledgement'] = (
        "Support from Paris Astronomical Data Centre (PADC) is acknowledged "
        "for hosting and providing access to the data."
    )

    return cdf


def to_cdf(date, cdf_out_path, version='01', real_fillval=-1e31):

    dtimes, freq, flux, snr = load_csv(date)

    flux = np.where(np.isnan(flux), real_fillval, flux)
    snr = np.where(np.isnan(snr), real_fillval, snr)

    cdf_filename = 'wi_wa_rad1_l3_akr_{}_v{}.cdf'.format(date.strftime('%Y%m%d'), version)
    # with pathlib
    # cdf_path = Path('./1999')
    # cdf_path = cdf_path / cdf_filename

    cdf_path = os.path.join(cdf_out_path, cdf_filename)
    # print(cdf_path)
    cdf_stem = cdf_path[-34:-4]
    print(cdf_stem)

    if os.path.isfile(cdf_path):
        os.remove(cdf_path)

    # Opening CDF object
    pycdf.lib.set_backward(False) # this is setting the CDF version to be used

    cdf = pycdf.CDF(cdf_path, '')

    # required settings for ISTP and PDS compliance
    cdf.col_major(True)  # Column Major
    # leave this for now
    # cdf.compress(pycdf.const.NO_COMPRESSION)  # No file level compression   
    
    cdf = set_istp_attr(cdf, cdf_path, cdf_stem, date)

    # SPASE
    cdf.attrs['spase_DatasetResourceID'] = ' '

    # SETTING VESPA GLOBAL ATTRIBUTES
    cdf.attrs['VESPA_dataproduct_type'] = 'DS>Dynamic Spectrum'
    cdf.attrs['VESPA_target_class'] = 'planet'
    cdf.attrs['VESPA_target_region'] = 'magnetosphere'
    cdf.attrs['VESPA_target_name'] = 'Earth'
    cdf.attrs['VESPA_feature_name'] = 'AKR>Auroral Kilometric Radiation'
    cdf.attrs['VESPA_instrument_host_name'] = 'Wind'
    cdf.attrs['VESPA_instrument_name'] = 'Waves'
    cdf.attrs['VESPA_receiver_name'] = 'RAD1'
    cdf.attrs['VESPA_measurement_type'] = 'phys.flux.density;em.radio'
    cdf.attrs['VESPA_access_format'] = 'application/x-cdf'
    cdf.attrs['VESPA_bib_reference'] = '10.1029/2021JA029425'

    vespa_date_min = date.floor('D')
    cdf.attrs['VESPA_time_min'] = vespa_date_min.strftime('%Y-%m-%dT%H:%M:%S')
    cdf.attrs['VESPA_time_max'] = (vespa_date_min + pd.Timedelta(1, 'D')).strftime('%Y-%m-%dT%H:%M:%S')
    cdf.attrs['VESPA_time_origin'] = 'Wind'
    cdf.attrs['VESPA_time_scale'] = 'TT'

    cdf.attrs['VESPA_time_sampling_step'] = '183' # 183.0 s
    cdf.attrs['VESPA_time_sampling_step_unit'] = 's'

    cdf.attrs['VESPA_spectral_range_min'] = '20'
    cdf.attrs['VESPA_spectral_range_max'] = '1040'
    cdf.attrs['VESPA_spectral_range_unit'] = 'kHz'
    
    # throughout freq range, what is min/max sampling step (whether lin or log)
    cdf.attrs['VESPA_spectral_sampling_step_min'] = '4' 
    cdf.attrs['VESPA_spectral_sampling_step_max'] = '136'
    cdf.attrs['VESPA_spectral_sampling_step_unit'] = 'kHz'

    # SETTING VARIABLES
    # time in datetime.datetime or pd.Timestamp if works
    # have kept as UNIX time for now - not sure about eg 'FILLVAL' if date in datetime.datetime or pd.Timestamp etc
    ts = (dtimes - np.datetime64('1970-01-01T00:00:00')) / np.timedelta64(1, 's')
    time_correct_format = list(map(datetime.utcfromtimestamp, ts))
    cdf.new('Epoch',
            data=time_correct_format,
            type=pycdf.const.CDF_TIME_TT2000,
            compress=pycdf.const.NO_COMPRESSION)
    cdf['Epoch'].attrs.new('VALIDMIN', data=time_correct_format[0], type=pycdf.const.CDF_TIME_TT2000)
    cdf['Epoch'].attrs.new('VALIDMAX', data=time_correct_format[-1], type=pycdf.const.CDF_TIME_TT2000)
    cdf['Epoch'].attrs.new('SCALEMIN', data=time_correct_format[0], type=pycdf.const.CDF_TIME_TT2000)
    cdf['Epoch'].attrs.new('SCALEMAX', data=time_correct_format[-1], type=pycdf.const.CDF_TIME_TT2000)
    cdf['Epoch'].attrs['CATDESC'] = "UTC time of each spectrum."
    cdf['Epoch'].attrs['FIELDNAM'] = "Epoch"
    cdf['Epoch'].attrs.new('FILLVAL', data=-9223372036854775808, type=pycdf.const.CDF_TIME_TT2000)
    cdf['Epoch'].attrs['LABLAXIS'] = "Epoch"
    cdf['Epoch'].attrs['UNITS'] = "ns"
    cdf['Epoch'].attrs['FORM_PTR'] = "CDF_TIME_TT2000"
    cdf['Epoch'].attrs['VAR_TYPE'] = "support_data"
    cdf['Epoch'].attrs['SCALETYP'] = "linear"
    cdf['Epoch'].attrs['MONOTON'] = "INCREASE"
    cdf['Epoch'].attrs['REFERENCE_POSITION'] = "Earth"
    cdf['Epoch'].attrs['SI_CONVERSION'] = "1.0e-9>s"
    cdf['Epoch'].attrs['UCD'] = "time.epoch"
    cdf['Epoch'].attrs['TIME_BASE'] = 'J2000'
    cdf['Epoch'].attrs['TIME_SCALE'] = 'TT'

    cdf.new('Frequency',
        data=freq,
        type=pycdf.const.CDF_REAL4,
        compress=pycdf.const.NO_COMPRESSION,
        recVary=False)
    cdf['Frequency'].attrs['CATDESC'] = "Central frequency of each step of the spectral sweep." # Central frequency of each step of the spectral sweep.
    cdf['Frequency'].attrs['DICT_KEY'] = "frequency"
    cdf['Frequency'].attrs['FIELDNAM'] = 'Frequency'
    cdf['Frequency'].attrs.new('FILLVAL', data=real_fillval, type=pycdf.const.CDF_REAL4)
    cdf['Frequency'].attrs['FORMAT'] = "E12.4"
    cdf['Frequency'].attrs['LABLAXIS'] = 'Frequency'
    cdf['Frequency'].attrs['UNITS'] = "kHz"
    cdf['Frequency'].attrs.new('VALIDMIN', data=20., type=pycdf.const.CDF_REAL4)
    cdf['Frequency'].attrs.new('VALIDMAX', data=1040., type=pycdf.const.CDF_REAL4)
    cdf['Frequency'].attrs['VAR_TYPE'] = "support_data"
    cdf['Frequency'].attrs['SCALETYP'] = "linear"
    cdf['Frequency'].attrs.new('SCALEMIN', data=20., type=pycdf.const.CDF_REAL4)
    cdf['Frequency'].attrs.new('SCALEMAX', data=1040., type=pycdf.const.CDF_REAL4)
    cdf['Frequency'].attrs['UCD'] = "em.freq"
    cdf['Frequency'].attrs['SI_CONVERSION'] = "1.0e3>Hz"

    # modify to flux density
    cdf.new('FLUX_DENSITY',
            data=flux.T,
            type=pycdf.const.CDF_REAL4,
            compress=pycdf.const.NO_COMPRESSION)
    cdf['FLUX_DENSITY'].attrs['CATDESC'] = "Flux density spectrogram measured on the Wind/WAVES spin-axis Z antenna." 
    cdf['FLUX_DENSITY'].attrs['DEPEND_0'] = 'Epoch'
    cdf['FLUX_DENSITY'].attrs['DEPEND_1'] = 'Frequency'
    cdf['FLUX_DENSITY'].attrs['DICT_KEY'] = 'electric_field>power'
    cdf['FLUX_DENSITY'].attrs['DISPLAY_TYPE'] = 'spectrogram'
    cdf['FLUX_DENSITY'].attrs['FIELDNAM'] = 'FLUX_DENSITY'
    cdf['FLUX_DENSITY'].attrs.new('FILLVAL', data=real_fillval, type=pycdf.const.CDF_REAL4)
    cdf['FLUX_DENSITY'].attrs['FORMAT'] = "E12.4"
    cdf['FLUX_DENSITY'].attrs['LABLAXIS'] = 'Flux Density @ 1AU'
    cdf['FLUX_DENSITY'].attrs['UNITS'] = 'Wm2Hz-1'
    cdf['FLUX_DENSITY'].attrs.new('VALIDMIN', data=1e-30, type=pycdf.const.CDF_REAL4)
    cdf['FLUX_DENSITY'].attrs.new('VALIDMAX', data=1e-15, type=pycdf.const.CDF_REAL4)
    cdf['FLUX_DENSITY'].attrs['VAR_TYPE'] = "data"
    cdf['FLUX_DENSITY'].attrs['SCALETYP'] = "log"
    # get relevant min/max flux
    cdf['FLUX_DENSITY'].attrs.new('SCALEMIN', data=1.69e-26, type=pycdf.const.CDF_REAL4)
    cdf['FLUX_DENSITY'].attrs.new('SCALEMAX', data=2.4e-18, type=pycdf.const.CDF_REAL4)
    cdf['FLUX_DENSITY'].attrs['UCD'] = "phys.flux.density;em.radio"

    # modify to SNR
    cdf.new('SNR',
            data=snr.T,
            type=pycdf.const.CDF_REAL4,
            compress=pycdf.const.NO_COMPRESSION)
    cdf['SNR'].attrs['CATDESC'] = "SNR spectrogram measured on the Wind/WAVES spin-axis Z antenna using 24 hr 5% background."
    cdf['SNR'].attrs['DEPEND_0'] = 'Epoch'
    cdf['SNR'].attrs['DEPEND_1'] = 'Frequency'
    cdf['SNR'].attrs['DICT_KEY'] = 'electric_field>power'
    cdf['SNR'].attrs['DISPLAY_TYPE'] = 'spectrogram'
    cdf['SNR'].attrs['FIELDNAM'] = 'SNR'
    cdf['SNR'].attrs.new('FILLVAL', data=real_fillval, type=pycdf.const.CDF_REAL4)
    cdf['SNR'].attrs['FORMAT'] = "E12.4"
    cdf['SNR'].attrs['LABLAXIS'] = 'SNR'
    cdf['SNR'].attrs['UNITS'] = 'dB'
    # check sensible valid min/max limits
    cdf['SNR'].attrs.new('VALIDMIN', data=10, type=pycdf.const.CDF_REAL4)
    cdf['SNR'].attrs.new('VALIDMAX', data=100, type=pycdf.const.CDF_REAL4)
    cdf['SNR'].attrs['VAR_TYPE'] = "support_data"
    cdf['SNR'].attrs['SCALETYP'] = "linear"
    # get min/max limits
    cdf['SNR'].attrs.new('SCALEMIN', data=16, type=pycdf.const.CDF_REAL4)
    cdf['SNR'].attrs.new('SCALEMAX', data=84, type=pycdf.const.CDF_REAL4)
    cdf['SNR'].attrs['UCD'] = "stat.snr;phys.flux.density;em.radio"

    cdf.close()

    return None


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('date0')

    parser.add_argument('date1')

    parser.add_argument('csv_path')

    args = parser.parse_args()

    date0 = pd.Timestamp(args.date0)
    date1 = pd.Timestamp(args.date1)
    csv_path = str(args.csv_path)

    dates = pd.date_range(date0, date1, freq='D')

    for d in dates:

        to_cdf(d, csv_path)

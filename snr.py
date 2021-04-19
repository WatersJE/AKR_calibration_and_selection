#! usr/bin/env python3
"""
Contains routines to calculate SNR for a set of source parameters
"""
import numpy as np

import background

def snr_sigma(bg_df, freq, antenna,
              k=1, beta=3e3, tau_ssp=0.154, tau_z=0.308):
    """
    Compute the signal fluctuation at a given frequency for a fiven
    antenna
    """
    ant_lab = 'AMPL_' + antenna

    bg = bg_df.loc[int(freq), ant_lab]

    if (antenna == 'S') or (antenna == 'SP'):
        sigma = (k * bg) / np.sqrt(beta * tau_ssp)
    else:
        sigma = (k * bg) / np.sqrt(beta * tau_z)

    return sigma


def snr_single(bg_df, amplitude, freq, antenna):
    """
    Compute the SNR for a single measurement (or array of measurements)
    """
    sigma = snr_sigma(bg_df, freq, antenna)

    snr_sing = amplitude / sigma

    snr_sing = 10 * np.log10(snr_sing)

    return snr_sing


def snr_average(bg_df, amplitudes, freq, antenna):
    """
    Compute the average SNR for a single measurement step 
    (set of 8 measurements of one frequency for RAD1)
    """

    assert antenna in ['S', 'SP', 'Z'], \
        'Antenna must be one of [S, SP, Z]'

    n_ampl = len(amplitudes)

    sigma = snr_sigma(bg_df, freq, antenna)

    avg_ampl = np.mean(amplitudes)

    avg_sigma = sigma / np.sqrt(n_ampl)

    snr_avg = avg_ampl / avg_sigma

    snr_avg = 10 * np.log10(snr_avg)

    return snr_avg


        


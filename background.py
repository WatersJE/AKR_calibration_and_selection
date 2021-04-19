#! usr/bin/env python3
"""
Module to remove background (QTN and galactic sources) from L2 data prior to
GP inversion and least squares calculations

eg find closest frequency if not present in bg_df

NB in future should amend this with routines for calculating background 
values itself! given some other quiet period - i.e past just loading in
current background and subtracting
"""
import os
import numpy as np
import pandas as pd


def new_background_file(l2_df, quant=0.05,
    freq_label='FREQ', custom_group=None):
    """
    Create and save file of background data corresponding to the 24hour period
    of data found in `l2_df`. Background values computed for each frequency
    using the quantile given in `quant` (5% by default)
    """
    agg_func = lambda x: np.quantile(x, quant)

    if custom_group is not None and isinstance(custom_group, list):

        bg_df = l2_df.groupby(custom_group).agg({'AMPL_S': agg_func,
                                            'AMPL_SP': agg_func,
                                            'AMPL_Z': agg_func})
    else:
            
        bg_df = l2_df.groupby(freq_label).agg({'AMPL_S': agg_func,
                                            'AMPL_SP': agg_func,
                                            'AMPL_Z': agg_func})

    return bg_df


def subtract_background(bg_df, amplitude, freq, antenna):
    """
    Remove the background voltage spectral density value for a single frequency


    """
    assert isinstance(antenna, str), 'antenna must be given as string'
    
    freq = int(freq)

    assert antenna in ['S', 'SP', 'Z'], \
        'Antenna must be one of [S, SP, Z]'
           
    bg_label = 'AMPL_' + antenna.upper()

    bg = bg_df.loc[freq, bg_label]

    amplitude_sub = amplitude - bg
        
    if np.any(amplitude_sub <= 0):
        # print('Ineligible value (< 0) after background subtraction '
        #       '- using minimum background value')
        amplitude_sub = np.where(amplitude_sub > 0, amplitude_sub, bg)

    return amplitude_sub


def nearest_frequency(freq, bg_freqs,
                      accpt_diff=5):
    """
    Find nearest available background frequency if `freq` not present
    
    Parameters
    ----------
    ...

    accpt_diff : int
        maximum allowed difference between frequencies before
        removing background
    """
    f_diffs = np.abs(np.array(bg_freqs) - freq)

    if np.min(f_diffs) <= accpt_diff:

        new_freq = bg_freqs[np.argmin(f_diffs)]

    else:
        # print('TEMPORARY FIX: AUTO SET DIFFERENCE TO NECESSARY VALUE')
        new_accpt = f_diffs.min()
        # print('NEW DIFFERENCE {}kHz'.format(new_accpt))
        
        new_freq = nearest_frequency(freq, bg_freqs, new_accpt)

    return new_freq


if __name__ == "__main__":

    os.chdir('../data/windwaves_background')

    
    

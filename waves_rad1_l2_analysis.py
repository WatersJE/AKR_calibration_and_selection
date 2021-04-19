import os
import pandas as pd 
import numpy as np
import logging


def setup_log(log_filename):

    logging.basicConfig(filename=log_filename,
        level=logging.INFO)

    return None


def get_sweep_datetime(waves_l2, sweep_i):
    """
    Given a header/data block of wind/waves L2 data,
    extract the datetime object corresponding to the start of the sweep cycle
    given by `sweep_i`

    Parameters
    ----------
    waves_l2: Waves_data object
        see read_wind_waves_file.read_l2_hres() output

    sweep_i: int
        index of sweep cycle of interest
        0 <= sweep_i <= max_sweep-1

    Returns
    -------
    date_time: datetime object
        the datetime of the beginning of the sweep cycle (cycle # sweep_i+1)
    """
    
    year = waves_l2.header[sweep_i]['YEAR']
    month = waves_l2.header[sweep_i]['MONTH']
    day = waves_l2.header[sweep_i]['DAY']
    
    date = pd.Timestamp(year, month, day)

    msec = waves_l2.header[sweep_i]['MSEC_OF_DAY']

    datetime = date + pd.Timedelta(msec, 'milli')

    return datetime


def sampled_freqs(l2_obj, n_freqs=32, raise_err=False):
    """
    Checks for consistent frequency sampling across 24 H L2 object

    n_freqs is expected number of channels

    if raise_err is false, returns boolean array [0, 1] with False (0) if
    the succeeding set of sampled frequencies (i+1) differs 
    """

    # freqs = np.array(l2_obj.data[0]['FREQ'])
    freqs = [np.array(l2_obj.data[i]['FREQ']) for i in range(len(l2_obj.data))]

    freqs = [np.array([np.int(f) for f in np.sort(np.unique(freqs_i))])
        for freqs_i in freqs]
    print(len(freqs))

    same_freqs = np.zeros(len(freqs))

    for i in range(len(freqs) - 1):

        if freqs[i].shape == freqs[i+1].shape:

            same_freqs[i] = True if np.allclose(freqs[i],
                freqs[i+1]) else False

        else:

            same_freqs[i] = False

    if raise_err:

        if np.all(same_freqs):

            return np.array(freqs[0]), same_freqs

        else:

            raise ValueError('Day contains sweeps that sample '
                'different frequencies')

    else:

        return np.array(freqs[0]), same_freqs


def convert_sweep_cycle_block(waves_l2, sweep_i):
    """
    (from prev GP retrieval function)
    See 20190627_notes.txt for details of progression from D_0 and cos(gamma)
    
    Insert check for post 04 aug 2000 - R changes

    (new)
    Transfer data from dictionary into pandas dataframe for one 
    header/data sweep cycle block for all antenna
    
    List of antennae channels used:

        Z is spin-axis-aligned antenna
        S is synthetic antenna (Wind/Waves SUM mode)
        SP is synthetic antenna with phase shift
    
    Dictionary keys are:
        ('VZPAL/TZPAL' for Z, 'VSPAL/TSPAL' for S/Sp)

    Parameters
    ----------
    waves_l2: Waves_data class object

    sweep_i: int
        Indexer used to access appropriate header/data block of waves_l2
        sweep # - 1


    Returns
    -------
    block: pandas DataFrame
        columns and corresponding original wind/waves table names:
            'AMPL' : either 'VSPAL' or 'VZPAL'
            'ANT' : sampled antenna either 'Z' or 'S' (Sum mode)
                i.e 'IANTEN' = 2 for 'S'
            'FREQ' : 'PALKHZ'
            'TIME' : 'IAMJCY'+'IHMSCY'+ either 'TSPAL' or 'TZPAL'
    """

    # Extract frequency bins (PalkHz)
    freqs = np.array(waves_l2.data[sweep_i]['FREQ'])
    
    # Checking number of measurement steps correspond to recorded value
    #(Could move into test_ function - tests all these values between header/data)
    if np.shape(freqs)[0] != waves_l2.header[sweep_i]['NPALIF']:
        
        raise AssertionError('Line 150 assertion raised (make_wav_l3_rad1.py)')

    # Extract measurement times and amplitudes for each antenna
    # S
    times_s = np.array(waves_l2.data[sweep_i]['TSPAL'])[::2]
    ampl_s = np.array(waves_l2.data[sweep_i]['VSPAL'])[::2]
    N_measurements_s = waves_l2.header[sweep_i]['NSPALF'] / 2
    
    # S_prime
    times_sp = np.array(waves_l2.data[sweep_i]['TSPAL'])[1::2]
    ampl_sp = np.array(waves_l2.data[sweep_i]['VSPAL'])[1::2]
    N_measurements_sp = waves_l2.header[sweep_i]['NSPALF'] / 2
    
    # Z
    times_z = np.array(waves_l2.data[sweep_i]['TZPAL'])
    ampl_z = np.array(waves_l2.data[sweep_i]['VZPAL'])
    N_measurements_z = waves_l2.header[sweep_i]['NZPALF']
    
    assert (N_measurements_s == N_measurements_sp) and \
           (N_measurements_s == N_measurements_z) and \
           (N_measurements_sp == N_measurements_z), \
           "Unequal number of measurements made between each antenna."
           
    N_samples = N_measurements_s

    freqs = np.repeat(freqs, N_samples)
    freqs = np.array([int(f) for f in freqs])
    
    # Pull sun angle (see DSUNCY in docs) and spin rate (see DSSPIN in docs)
    sun_angle = waves_l2.header[sweep_i]['SUN_ANGLE']
    spin_rate = waves_l2.header[sweep_i]['SPIN_RATE']
    # 1-indexed sweep cycle number
    i_sweep = waves_l2.header[sweep_i]['ISWEEP']
    # L2 intensity units tag
    i_unit = waves_l2.header[sweep_i]['IUNIT']
    # solar top quality (?) tag
    k_spin = waves_l2.header[sweep_i]['KSPIN']
    # data acquisition mode
    mode = waves_l2.header[sweep_i]['MODE']
    # antenna configuration tag
    i_ant = waves_l2.header[sweep_i]['IANTEN']
    # polarisation flag
    i_pol = waves_l2.header[sweep_i]['IPOLA']
    # eq dipole used
    i_dipole_xy = waves_l2.header[sweep_i]['IDIPXY']
    # sweep cycles duration (seconds)
    sweep_dur = waves_l2.header[sweep_i]['SDURCY']
    # measurement step duration (duration to sample one frequency)
    step_dur = waves_l2.header[sweep_i]['SDURPA']
    
    # number of measurement steps (sets of 8 frequency samples)
    n_meas_steps = waves_l2.header[sweep_i]['NPALCY']
    # number of frequencies in measurement step (assumed 1)
    n_freq_per_step = waves_l2.header[sweep_i]['NFRPAL']
    # number of total measurement steps 
    n_meas_step_f = waves_l2.header[sweep_i]['NPALIF']
    # number of unique frequencies measured
    n_freq = waves_l2.header[sweep_i]['NFREQ']
    
    # get number of measurement steps for S/SP and Z antennae
    n_samp_per_step_S = waves_l2.header[sweep_i]['NSPALF']
    n_samp_per_step_Z = waves_l2.header[sweep_i]['NZPALF']    
    
    swcy_tag = 0
    try:
        assert n_freq_per_step == 1, \
            "More than one frequency sampled per measurement step"
    except AssertionError:
        swcy_tag = 1
    
    # Still have datetime of start of sweep cycle
    sweep_start_date = get_sweep_datetime(waves_l2, sweep_i)

    # unpacking datetimes from time arrays
    datetime_s, datetime_sp, datetime_z = [np.array([sweep_start_date + 
        pd.Timedelta(t, 's') for t in times]) for times in [times_s, times_sp, times_z]]
    
    block = pd.DataFrame({'FREQ': freqs,
                          'TIME_S': times_s,
                          'TIME_SP': times_sp,
                          'TIME_Z': times_z,
                          'AMPL_S': ampl_s,
                          'AMPL_SP': ampl_sp,
                          'AMPL_Z': ampl_z,
                          'DATETIME_S': datetime_s,
                          'DATETIME_SP': datetime_sp,
                          'DATETIME_Z': datetime_z})

    block['SWEEP'] = block.shape[0] * [sweep_i]
    
    sweep_out = {"sweep_tag": swcy_tag,
                 "i_sweep": i_sweep,
                 "i_units": i_unit,
                 "k_spin": k_spin,
                 "data_acq_mode": mode,
                 "ant_config": i_ant,
                 "pol_present": i_pol,
                 "eq_dipole": i_dipole_xy,
                 "sweep_dur": sweep_dur,
                 "step_dur": step_dur,
                 "start_date": sweep_start_date,
                 "sun_angle": sun_angle,
                 "spin_rate": spin_rate,
                 "n_s_samples_per_step": n_samp_per_step_S,
                 "n_z_samples_per_step": n_samp_per_step_Z,
                 "n_steps": n_meas_steps,
                 "n_freqs": n_freq,
                 "data": block}
    
    # Apply tag if sweep cycle differs from base expectation
    try:
        check_sweep_cycle(sweep_out)                          
    except ValueError:
        logging.info('Sweep {} header raised validation flag - see above.'.format(sweep_i))
        sweep_out['sweep_tag'] = 1
        
    return sweep_out
    

def check_sweep_cycle(sweep_cycle_dict):
    """
    Ensure that appropriate tags are present within the sweep cycle header
    file - e.g polarisation present, optimal antennae used, correct units etc
    
    Parameters
    ----------
    sweep_cycle_dict : dict
        contains relevant header parameters and amplitude/time data for each
        channel
        
    Returns
    -------
    None
    
    
    NB - may require changes, asserts not necessary for less stringent
    requirements eg eq_dipole
    
    """
    # Check receiver amplitudes in microV^2 Hz^-1
    if not sweep_cycle_dict['i_units'] == 3:

        logging.info('Receiver amplitude data are not in the correct units.')

        raise ValueError()

    # Check solar top quality
    if not sweep_cycle_dict['k_spin'] == 0:
        
        logging.info('Poor solar top quality - interpolation has been applied. (?)')
        raise ValueError()
    
    # Check data acquired using measure/list mode
    if not sweep_cycle_dict['data_acq_mode'] == 3:
        logging.info('Frequencies sampled with incorrect data acquisition mode')
        raise ValueError()
        
    # Check X and Z antenna used in SUM mode
    if not sweep_cycle_dict['ant_config'] == 2:
        logging.info('X and Z antennae not configured in SUM mode')
        raise ValueError()
        
    # Check polarisation present
    if not sweep_cycle_dict['pol_present'] == 1:
        logging.info('Polarisation not present for this sweep cycle.')
        raise ValueError()

    # Check appropriate, longer equatorial antenna used
    if not sweep_cycle_dict['eq_dipole'] == 1:
        logging.info('Shorter equatorial plane Y-antenna used')
        raise ValueError()
        
    # Check same number of samples per measurement step observed for S and SP
    # as for Z
    if not sweep_cycle_dict['n_s_samples_per_step'] == \
           2 * sweep_cycle_dict['n_z_samples_per_step']:
        
        logging.info('Unequal number of measurements from Z and S/Sprime antennae')
        raise ValueError()
    
    return None


def antenna_gain_ratio(date, R_initial=5.,
    ant_length_0=50., ant_length_1=32.):
    """
    Return the appropriate value for the antenna gain ratio R. Checks if date
    is before or after date of collision (3rd August 2000).

    NB thought to have secondary collision after this date, but
    unconfirmed/uncalibrated

    Parameters
    ----------
    date : str, pd.Timestamp

    Returns
    -------
    R : float
    """
    
    if isinstance(date, str):
        
        date = pd.Timestamp(date)

    elif not isinstance(date, pd.Timestamp):

        raise ValueError('Date must be either string or pandas Timestamp')

    R_pre = 5.
    R_post = R_pre * (ant_length_1 / ant_length_0)

    R = R_pre if date <= pd.Timestamp('20000803') else R_post

    return R


def concat_l2_sweeps(l2_object, n_sw, log_filename=None):
    """
    Convert sweep cycle data object into dataframe for the whole day
    
    NB will have to convert times to datetimes prior to concatenation, 
    as times are stored as seconds from sweep cycle start
    """
    if log_filename is not None:

        setup_log(log_filename)

    dfs = np.empty(n_sw, dtype=object)

    for i in range(n_sw):
        
        sweep_dict = convert_sweep_cycle_block(l2_object, i)
        
        sweep_df = sweep_dict['data']

        sweep_df['sweep_start_date'] = np.repeat(sweep_dict['start_date'],
            sweep_df.shape[0])

        sweep_df['sun_angle'] = np.repeat(sweep_dict['sun_angle'],
            sweep_df.shape[0])

        sweep_df['spin_rate'] = np.repeat(sweep_dict['spin_rate'],
            sweep_df.shape[0])
        # print(i, sweep_dict['i_sweep'])
        # sweep_df['SWEEP'] = np.repeat(sweep_dict['i_sweep'],
        #     sweep_df.shape[0])

        sweep_df['sweep_flag'] = np.repeat(sweep_dict['sweep_tag'],
        sweep_df.shape[0])

        if sweep_dict['sweep_tag'] == 0:

            dfs[i] = sweep_df    
    
        else:

            dfs[i] = None

    if np.all([d is None for d in dfs]):
        logging.info('\nAll sweeps for {} are invalid - storing empty dataframe'.format(
            sweep_dict['start_date'].strftime('%d %b %Y (DOY %j)')
        ))

        raise NotImplementedError('Raise error to ensure checking for this'
            'when calling concat_l2_sweeps, and creating appropriate data'
            'structure/message.')

    else:

        day_df = pd.concat(dfs, axis=0, ignore_index=True)
    
    return day_df

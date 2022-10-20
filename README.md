# Wind/Waves AKR Calibration and Selection Pipeline

This simple pipeline takes L2 raw Wind/Waves RAD1 data, performs a full calibration assuming an AKR source in the direction of Earth, and selects and retains data that meets a numerical threshold based on the variability of the spin signal. Full details of the calibration and selection can be found in _(DOI)_. This documentation describes the scripts and relevant inputs and data needed to process L2 binary data in this way.

The repository is self-contained, but requires Wind ephemeris data (namely GSE latitude and radial distance from Earth) as well as L2 data for the relevant observations. **The filepaths for this data should be changed manually in the following files _(cmd line args in future)_:** 

* `wind_waves_akr_flux.py`
  * Input filepath for L2 data (at line 274 (filepath) and line 64 (filename) ... include this detail?)
  * Output filepath for L3 data
* `calibration.py`
  * Input filepath for Wind ephemeris data
  
With filepaths specified, the pipeline can be run from the command line with a call to `wind_waves_akr_flux.py`, as described below.

## Programs
* `wind_waves_akr_flux.py` - the main program that creates a .csv file of masked AKR data for a single L2 binary file - corresponding to a 24 hour observation. After changing relevant input (L2 binary data, ephemeris) and output (L3 .csv file) filepaths, this can be run from the terminal with the date (YYYYMMDD) and version (?) as follows:
	`python3 wind_waves_akr_flux.py 20000101 01`

* `waves_rad1_l2_analysis.py` - data manipulation and validation routines. Mostly helper functions. 

* `background.py` - routines for creating background spectra for substraction prior to calibration

* `calibration.py` - routines for calibration and further data manipulation. Requires manual change to filepath for ephemeris data _(change to control by cmd input)_. Needs specific csv file and data labels at present.

* `snr.py` - routines to calculate signal-to-noise ratio (SNR) for each spin and resulting average spectrum.

* `read_wind_waves_l2.py` - reads binary files of RAD1 L2 data

* `wind_akr_mask_to_cdf.py` 

## Dependencies
* `pathlib`
* `numpy`
* `pandas`
* `spacepy`


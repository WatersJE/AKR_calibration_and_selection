#Wind/Waves AKR Calibration and Selection Pipeline

This pipeline takes L2 raw Wind/Waves RAD1 data, performs a full calibration assuming an AKR source in the direction of Earth, and selects and retains data that meets a numerical threshold based on the variability of the spin signal. Full details of the calibration and selection can be found in _(DOI)_. This documentation describes the scripts and relevant inputs and data needed to process L2 binary data in this way.

The repository is self-contained, but requires Wind ephemeris data (namely GSE latitude and radial distance from Earth) as well as L2 data for the relevant observations. **The filepaths for this data should be changed manually in the following files:** 

* `wind_waves_akr_flux.py`
* `calibration.py`
* _(check for others)_

##Programs
* `wind_waves_akr_flux.py` is the main program that creates a CSV file of masked AKR data for a single L2 binary file - corresponding to a ~24 hour observation. After changing relevant input (L2 binary data, ephemeris) and output (L3 .csv file) filepaths, this can be run from the terminal with the date (YYYYMMDD) and version (?) as follows:
	`python3 wind_waves_akr_flux.py 20000101 01`

* `waves_rad1_l2_analysis.py`

* `background.py`

* `calibration.py`

* `snr.py`

* `read_wind_waves_l2.py`



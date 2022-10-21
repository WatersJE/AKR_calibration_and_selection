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

## Installation

Notes and to do:

* git clone etc - describe installation
* include example scripts for running program
* simplest regression test needs to involve single L2 file -> L3 output
* how is ephemeris that has been used retrieved? (SPDF etc) - needs description - what to retrieve etc, how stored/expected to be read

  * The download of ephemeris could be automatic in the script - future

* remove 'version' from cmd line argument - this indicates the versionof the prorgam itself, so should be set witin program

* If wanted to store large cdf for all Wind ephemeris in repository - can use git LFS for large files - but description needed for ephemeris whichever solution used

- baptiste will fork repo and include files for eg date to test
- provide ephemeris file and expected output for eg date (20000101)

- CDF file creation - VESPA attributes for scalemin/max relates to plotting and colorbar scales - should be representative of all data

## Retrieving Wind Ephemeris

* Visit the [NASA SPDF spacecraft locator form](https://sscweb.gsfc.nasa.gov/cgi-bin/Locator.cgi) and complete with the following:

### Required Settings
#### Spacecraft/Time Range Selection
* In panel "Satellites" choose Wind. Set the relevant time range (note both inclusive)
* Note that L2 files often contain spectra that are measured in the day preceding the one labelled by the L2 file - so ephemeris data from the previous day will be required
* for e.g. 1st Jan 2000:
  *  Start Time: 1999/12/31
  *  Stop Time:	2000/01/01


#### Output Options
 
* For the purposes of the pipeline, only the GSE latitude and the radial distance from Earth are required - from the output options screen select at least the following:
   * "GSE LAT" in main options
   * "Radial Distance" under values in additional options

### Optional Settings
#### Output Format and Units Options   

* Select date format *(need to insure suggestion here works with date parser in pipeline)*
  * Ensure time selected to second resolution e.g. hh:mm:ss
 
* Choose "Earth Radii" (default) for the Distance units, with $\geq$ 2 decimal places

* Degrees Format, Direction/Range can both be kept as default values (-dd.ddd... with 2 decimal places; Lat (-90,+90), Long(0,360) respectively)

* Choose "Text" for Output Format Options

### Input Summary
The data selection can be checked via the Input Summary page

### Execution Options
Use these buttons to access the data.

*For now - will need to save the ephemeris file in a place and specify the filepath specifically in calibration.py.*

*Check also that the csv file saved with default text headers etc reads without trouble (as determined by the read_csv call in calibration.py)*
	

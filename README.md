# Wind/Waves AKR Calibration and Selection Pipeline

This simple pipeline takes L2 raw Wind/Waves RAD1 data, performs a full calibration assuming an AKR source in the direction of Earth, and selects and retains data that meets a numerical threshold based on the variability of the spin signal. Full details of the calibration and selection can be found in [Waters et al. (2021)](https://onlinelibrary.wiley.com/doi/10.1029/2021JA029425). This documentation describes the scripts and relevant inputs and data needed to process L2 binary data in this way.

The repository is self-contained, but requires Wind ephemeris data (namely GSE latitude and radial distance from Earth) as well as L2 data for the relevant observations. The directories, for L2 (raw, input) and L3 (AKR-selected, output) data, and the full filepath for Wind ephemeris data are passed as command line arguments to `wind_waves_akr_flux.py`, as well as the date to be processed in format `YYYYMMDD`. Instructions for installation and usage, as well as descriptions of each program, are given below. 

## Installation

Copy the files of this GitHub repository to your local machine using the green "<> Code" dropdown button on the [repository webpage](https://github.com/WatersJE/WindWaves_AKR_calibration_selection). You can:
* Clone the repository via HTTPS, SSH or with the GitHub CLI
* Open the repository with the GitHub Desktop app
* Download a ZIP file of the repository

## Usage Instructions

* To do - include an example script showing basic usage in command line
Once the Wind/WAVES L2-level file(s) relevant for processing has been accessed and stored locally, as well as the ephemeris files (see **Retrieving Wind Ephemeris** below), the AKR selection program can be run by entering the following in the command line (**need to include specific files for example**):
``` python3 wind_waves_akr_flux.py YYYYMMDD /path/to/L2/data/directory /path/to/ephemeris/data /path/to/AKR/output/directory ```

### Retrieving Wind Ephemeris

* Visit the [NASA SPDF spacecraft locator form](https://sscweb.gsfc.nasa.gov/cgi-bin/Locator.cgi) and complete with the following:

#### Required Settings
##### Spacecraft/Time Range Selection
* In panel "Satellites" choose Wind. Set the relevant time range (note both inclusive)
* Note that L2 files often contain spectra that are measured in the day preceding the one labelled by the L2 file - so ephemeris data from the previous day will be required
* for e.g. 1st Jan 2000:
  *  Start Time: 1999/12/31
  *  Stop Time:	2000/01/01


##### Output Options
 
* For the purposes of the pipeline, only the GSE latitude and the radial distance from Earth are required - from the output options screen select at least the following:
   * "GSE LAT" in main options
   * "Radial Distance" under values in additional options

#### Optional Settings
##### Output Format and Units Options   

* Select date format as follows
  * "Date": Choose yy/mm/dd
  * "Time": Choose hh:mm:ss
 
* Choose "Earth Radii" (default) for the Distance units, with $\geq$ 2 decimal places

* Degrees Format, Direction/Range can both be kept as default values (-dd.ddd... with 2 decimal places; Lat (-90,+90), Long(0,360) respectively)

* Choose "CDF" for Output Format Options

#### Input Summary
The data selection can be checked via the Input Summary page

#### Execution Options
Use these buttons to access the data. Save the CDF file appropriately and provide the filepath to the `wind_waves_akr_flux` program as a command line argument (see Usage instructions)

### Dependencies
The following software and Python libraries are required to be installed to use these programs:
* `CDF`
* `Python` $\geq$ 3.8
* `numpy`
* `pandas`
* `spacepy`

## Programs
* `wind_waves_akr_flux.py` - the main program that creates a .csv file of masked AKR data for a single L2 binary file - corresponding to a 24 hour observation. After changing relevant input (L2 binary data, ephemeris) and output (L3 .csv file) filepaths, this can be run from the terminal with the date (YYYYMMDD) and version (?) as follows:
	`python3 wind_waves_akr_flux.py 20000101 01`

* `waves_rad1_l2_analysis.py` - data manipulation and validation routines. Mostly helper functions.

* `background.py` - routines for creating background spectra for substraction prior to calibration

* `calibration.py` - routines for calibration and other relevant data processing. 

* `snr.py` - routines to calculate signal-to-noise ratio (SNR) for each spin and resulting average spectrum.

* `read_wind_waves_l2.py` - reads binary files of RAD1 L2 data

* `wind_akr_mask_to_cdf.py` 
	
---
Notes and to do:

* include example scripts for running program
* simplest regression test needs to involve single L2 file -> L3 output
* The download of ephemeris could be automatic in the script - future
  * If wanted to store large cdf for all Wind ephemeris in repository - can use git LFS for large files - but description needed for ephemeris whichever solution used

- baptiste will fork repo and include files for eg date to test
- provide ephemeris file and expected output for eg date (20000101)

- CDF file creation - VESPA attributes for scalemin/max relates to plotting and colorbar scales - should be representative of all data
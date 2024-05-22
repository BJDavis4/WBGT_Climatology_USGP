The scripts in this directory were used to generate the data and analysis for Davis et al 2024: Climatology of Wet Bulb Globe Temperature in the United States Great Plains, JAMC.

1) Needed Data <br>
Hourly Oklahoma Mesonet data for 1998-2020 for all hours in .mdf format. <br>
	Variables Needed: <br>
	2 meter temperature (TAIR) <br>
	2 meter relative humidity (RELH) <br>
	2 meter wind speed (WS2M) <br>
	10 meter wind speed (WSPD) <br>
	Solar radiation (SRAD) <br>
	Surface pressure (PRES) <br>
        geoinfo.csv <br>

Hourly ERA5 Reanalysis data for 1960-2020 for all hours in NetCDF format (Can use grib if you want to modify the code that reads in files)
	Variables needed:
	2 meter temperature (t2m)
	2 meter dewpoint temperature (d2m)
	10 meter u component of the wind (u10)
        10 meter v component of the wind (v10)
	Mean surface downward short-wave radiation flux (msdwswrf)
	Surface pressure (sp)

Needed packages
Python3 (I use the anaconda distibution)
	numpy
	xarray
	matplotlib
	mpl_toolkits.basemap; Basemap
	glob
	scipy
	netCDF4
        datetime
        multiprocessing
        dask
        numba
        mkl
Liljegren Cython code from [link cited in paper]

You must process the Oklahoma Mesonet data before you can calculate WBGT from the ERA5 reanalysis. The 2 meter wind speed is calculated based on the climatology of the roughness length from the Oklahoma Mesonet Observations.

2) Process Oklahoma Mesonet Data
Run the following scripts in order. Modify file paths and date ranges to match the data you have.
python mod_mdf.py --> this script calculated WBGT at each station and connects the lat/lon of each mesonet 
python grid_mdf.py --> this script regrids the Oklahoma Mesonet observations to match the ERA5 grid
python Z0Climo.py --> this script calculates the roughness length (z0) from the gridded Oklahoma Mesonet data and calculates the climatology for each hour of each month across the data availability period (1998-2020)

3) Process ERA5 reanalysis
Run the following scripts in order. Modify file paths and date ranges to match the data you have. The date range does not have to match that of the Oklahoma Mesonet data
python Calc2mWind.py --> This script calculates 2 meter wind speed from the 10 meter wind components and the roughness length
python intperp_SRAD.py --> This script interpolates the solar radiation (msdwswrf) to be centered on HH:00. Without this msdwswrf is averaged over an hour centered on HH:30 and results in a temporal mismatch with the other data
python WBGT_ERA5.py --> Calculates WBGT from the ERA5 reanalysis using Dimiceli
python WBGT_Liljegren.py --> Calculates WBGT from the ERA5 reanalysis using Liljegren
python WBGT2NA.py --> Trims WBGT dataset to the North American domain to reduce memory when reading in files during analysis
python WBGT2Cat.py --> converts WBGT to WBGT uniform categories. Need to run separately for Dimiceli and Liljegren
python WBGT2ModCat.py --> converts WBGT to WBGT regional categories. Need to run separately for Dimiceli and Liljegren

4) Analyses
python CalcClimo.py --> Run this script for each ERA5 variable to calculate the climatology for each variable. This will create two output files for each month of the year, one for the mean, and one for the standard deviation
python Plot_climo.py --> This script calculates the bias in the ERA5 vs Oklahoma Mesonet WBGT, and assumes the Oklahoma Mesonet is the "truth" and plots a map of each dataset and the difference between the two.
python RMSD.py --> This script calculates the RMSD between the ERA5 WBGT and the Oklahoma Mesonet WBGT
python Plot_Time_series.py --> This script plots the diurnal cycle of the mean and standard deviation of two WBGT datasets for each month as shown in Figure 2
python WBGT_pct.py --> Calculates various percentiles of the WBGT climatology for each day of the year
Python CatClimo.py --> This script does the remaining analysis in the paper for Figures 3-13 as well as several analyses that did not make it into the final draft of the manuscript.

5) Other necessary scripts that are called from those listed above
read_mdf.py
WBGT_func.py
netCDF_mods.py
Python_mods.py
MonteCarlo.py
wet_bulb.py
Black_Globe.py
WBGT.cpython-38-x86_64-linux-gnu.so (This is a compiled cython file. Your file name may vary)
coszenith.cpython-38-x86_64-linux-gnu.so (This is a compiled cython file. Your file name may vary) 


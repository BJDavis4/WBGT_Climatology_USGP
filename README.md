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

Hourly ERA5 Reanalysis data for 1960-2020 for all hours in NetCDF format (Can use grib if you want to modify the code that reads in files) <br>
	Variables needed: <br>
	2 meter temperature (t2m) <br>
	2 meter dewpoint temperature (d2m) <br>
	10 meter u component of the wind (u10) <br>
        10 meter v component of the wind (v10) <br>
	Mean surface downward short-wave radiation flux (msdwswrf) <br>
	Surface pressure (sp) <br>

Needed packages <br>
Python3 (I use the anaconda distibution) <br>
	numpy <br>
	xarray <br>
	matplotlib <br>
	mpl_toolkits.basemap; Basemap <br>
	glob <br>
	scipy <br>
	netCDF4 <br>
        datetime <br>
        multiprocessing <br>
        dask <br>
        numba <br>
        mkl <br>
Liljegren Cython code from https://github.com/QINQINKONG/PyWBGT

You must process the Oklahoma Mesonet data before you can calculate WBGT from the ERA5 reanalysis. The 2 meter wind speed is calculated based on the climatology of the roughness length from the Oklahoma Mesonet Observations. <br>

2) Process Oklahoma Mesonet Data <br>
Run the following scripts in order. Modify file paths and date ranges to match the data you have. <br>
python mod_mdf.py --> this script calculated WBGT at each station and connects the lat/lon of each mesonet <br>
python grid_mdf.py --> this script regrids the Oklahoma Mesonet observations to match the ERA5 grid <br>
python Z0Climo.py --> this script calculates the roughness length (z0) from the gridded Oklahoma Mesonet data and calculates the climatology for each hour of each month across the data availability period (1998-2020) <br>

3) Process ERA5 reanalysis <br>
Run the following scripts in order. Modify file paths and date ranges to match the data you have. The date range does not have to match that of the Oklahoma Mesonet data <br>
python Calc2mWind.py --> This script calculates 2 meter wind speed from the 10 meter wind components and the roughness length <br>
python intperp_SRAD.py --> This script interpolates the solar radiation (msdwswrf) to be centered on HH:00. Without this msdwswrf is averaged over an hour centered on HH:30 and results in a temporal mismatch with the other data <br>
python WBGT_ERA5.py --> Calculates WBGT from the ERA5 reanalysis using Dimiceli<br>
python WBGT_Liljegren.py --> Calculates WBGT from the ERA5 reanalysis using Liljegren <br>
python WBGT2NA.py --> Trims WBGT dataset to the North American domain to reduce memory when reading in files during analysis <br>
python WBGT2Cat.py --> converts WBGT to WBGT uniform categories. Need to run separately for Dimiceli and Liljegren<br>
python WBGT2ModCat.py --> converts WBGT to WBGT regional categories. Need to run separately for Dimiceli and Liljegren<br>

4) Analyses<br>
python CalcClimo.py --> Run this script for each ERA5 variable to calculate the climatology for each variable. This will create two output files for each month of the year, one for the mean, and one for the standard deviation<br>
python Plot_climo.py --> This script calculates the bias in the ERA5 vs Oklahoma Mesonet WBGT, and assumes the Oklahoma Mesonet is the "truth" and plots a map of each dataset and the difference between the two.<br>
python RMSD.py --> This script calculates the RMSD between the ERA5 WBGT and the Oklahoma Mesonet WBGT<br>
python Plot_Time_series.py --> This script plots the diurnal cycle of the mean and standard deviation of two WBGT datasets for each month as shown in Figure 2<br>
python WBGT_pct.py --> Calculates various percentiles of the WBGT climatology for each day of the year<br>
Python CatClimo.py --> This script does the remaining analysis in the paper for Figures 3-13 as well as several analyses that did not make it into the final draft of the manuscript.<br>

5) Other necessary scripts that are called from those listed above <br>
read_mdf.py<br>
WBGT_func.py<br>
netCDF_mods.py<br>
Python_mods.py<br>
MonteCarlo.py<br>
wet_bulb.py<br>
Black_Globe.py<br>
WBGT.cpython-38-x86_64-linux-gnu.so (This is a compiled cython file. Your file name may vary)<br>
coszenith.cpython-38-x86_64-linux-gnu.so (This is a compiled cython file. Your file name may vary)<br>


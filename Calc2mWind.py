import numpy as np
from netCDF4 import Dataset
from Python_mods import Calc_ws, calc_uk, interp_ws_log
from netCDF_mods import *

#Variables for reading in data
Start_year=1959
End_year=2020	
Start_month=1
End_month=12
input_path='/data/deluge/scratch/ERA5-Ben/2D/hourly'
z0_path='/data/deluge/observations/OKMesonet'
#output_path='/data/deluge/scratch/ERA5-Ben/2D/hourly'
u_var='u10'
v_var='v10'
z0_var='z0'

#Variables for interpolation
#z0=0.1
z2=2.0
z10=10.0

#Variables for saving to file
Var_name='ws2m'
Var_unit="m s**-1"
Var_long_name="2 metre wind speed"
output_path='/data/deluge/scratch/ERA5-Ben/2D/hourly'


for y in range(Start_year,End_year+1):
    for month in range(Start_month,End_month+1):
        m='{:02d}'.format(month)

#Read in 10m wind components
        f='{}/{}/{}.{}{}.nc'.format(input_path,u_var,u_var,y,m)
        Varn, Varnm, lat, lon = ReadNetCDF(f,['time',u_var])
        time=Varn[0]
        u=Varn[1]

        f='{}/{}/{}.{}{}.nc'.format(input_path,v_var,v_var,y,m)
        Varn, Varnm, lat, lon = ReadNetCDF(f,['time',v_var])
        v=Varn[1]

#Calculate and test wind speed
        #u=np.array([3])
        #v=np.array([4])
        #w=np.array([1])
        C=Calc_ws(u,v)
        
        print(u[0,52,72])
        print(v[0,52,72])
        print(C[0,52,72])

#Read in z0 climo for month
        f='{}/{}/{}.{}.mean.nc'.format(z0_path,'z0','z0',m)
        Varn, Varnm, latz0, lonz0 = ReadNetCDF(f,['time',z0_var])
        z0=Varn[1]
#        print(z0)
#        print(z0.shape)
#        exit()

        if month==2:
            hh=C.shape[0]
            z0=z0[:hh,:,:]
#interpolate to 2 m
        ws2m=interp_ws_log(z2,z0,z10,C)

        print(np.nanmean(C[0]))
        print(np.nanmean(ws2m[0]))
        print(z0[0][52][72])
        #exit()
#write to file

        print('Saving ws2m to file')
        lat=lat[:]
        #lat.shape=((1,lat.shape[0]))
        lon=lon[:]
        #lon.shape=((1,lon.shape[0]))
        #Month_time.shape=((1,Month_time.shape[0]))
        Dims=[lon,lat,time]
        DimsName=['longitude','latitude','time']
        DimsUnits=['degrees_east','degrees_north','hours since 1900-01-01 00:00:00.0']
        DimsLN=['longitude','latitude','time']
        Vars=[ws2m]
        VarsNames=[Var_name]
        VarsUnits=[Var_unit]
        VarsLN=[Var_long_name]
        VarsDims=[('time','latitude','longitude')]
        output_dir='{}/{}'.format(output_path,Var_name)
        filename='{}.{}{}.nc'.format(Var_name,y,m)

        WriteNetCDF(Dims,DimsName,DimsUnits,DimsLN,Vars,VarsNames,VarsUnits,VarsLN,VarsDims,output_dir,filename)



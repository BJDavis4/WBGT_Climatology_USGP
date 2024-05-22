#######################
# Need to make sure changes for chunks are good
import os
os.environ["OMP_NUM_THREADS"] = "1" # export OMP_NUM_THREADS=1
os.environ["OPENBLAS_NUM_THREADS"] = "1" # export OPENBLAS_NUM_THREADS=1
os.environ["MKL_NUM_THREADS"] = "1" # export MKL_NUM_THREADS=1
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" # export VECLIB_MAXIMUM_THREADS=1
os.environ["NUMEXPR_NUM_THREADS"] = "1" # export NUMEXPR_NUM_THREADS=1

import multiprocessing; multiprocessing.Pool(processes=1)
import numpy as np
from Python_mods import *
from WBGT_func import *
#from netCDF_mods import *
import glob
from datetime import datetime

import dask
import dask.array as da
from dask.distributed import Client
import xarray as xr
from numba import njit,vectorize
#import pyximport; pyximport.install(pyimport=True)
from WBGT import WBGT_Liljegren,fdir
from coszenith import coszda, cosza
import mkl; mkl.set_num_threads(1)

Start_year=1960
End_year=2020
Start_month=1
End_month=12
Num_hours=24
Num_lat=361
Num_lon=720
Chunk_size=4
num_chunks=8
max_workers=16

input_path='/data/deluge/scratch/ERA5-Ben/2D/hourly'
WS2m_var='ws2m'
t2m_var='t2m'
d2m_var='d2m'
Srad_var='msdwswrf'
pres_var='sp'

CorrectSRAD=False
SRADavgb=False #Set to false if using linear fit
CorrectWS2M=True
WS2MTF=True #True if wind speed is at 2 meters, False if at 10 meters

biaspath="/data/deluge/scratch/ERA5-Ben/2D/hourly/bias"
BiasVars=['SRADbias']

Var_name='WBGT'
Var_unit="K"
Var_long_name="Wet Bulb Globe Temperature"
output_path='/data/deluge/scratch/ERA5-Ben/2D/hourly'

def WriteNetCDF(Array,var_name,fill_value,attrs,coord_attrs,filepath):
    print('Computing dask array')
    time_start=datetime.now()
    #from dask.distributed import Client
    #if __name__ == '__main__':
     #client=Client(memory_limit='32GB',n_workers=1, threads_per_worker=2)
     #list_workers=list(client.scheduler_info()['workers'])

    data=Array.compute(num_workers=max_workers)#workers=list_workers[:4])
    #data=client.compute(Array)
    print(data[var_name].values[6::24,110,525])
    time_end=datetime.now()
    time_total=time_end-time_start
    print('Done computing dask array. Time:',time_total)
    #exit()
    data[var_name]=xr.where(data[var_name].isnull(),fill_value,data[var_name])
    data[var_name]=data[var_name].assign_attrs(attrs)

    for dim in data[var_name].dims:
        data.coords[dim]=data.coords[dim].assign_attrs(coord_attrs[dim])

    dataMax=np.nanmax(data[var_name].values)
    dataMin=np.nanmin(data[var_name].values)
    add_offset=dataMin
    scale_factor=(dataMax-dataMin)/((2**15)-1)

    print('------------------')
    print('add_offset:',add_offset)
    print('scale_factor',scale_factor)
    print('missing_value',fill_value)
    print('------------------')

    encoding={}
    encoding[var_name]={'add_offset':add_offset,'scale_factor':scale_factor,'missing_value':fill_value,'dtype':'int16'}
#    print(encoding)
#    print(data[var_name].attrs)
#    print(data[var_name].coords['latitude'].attrs)

    print('writing file')
    time_start=datetime.now()
    data.to_netcdf(filepath,encoding=encoding)
    time_end=datetime.now()
    time_total=time_end-time_start
    print('Done writing. Time:',time_total)

start_datetime=datetime.now()
for y in range(Start_year,End_year+1):
    if y==1959:
        Start_month=2
    else:
        Start_month=1
    for month in range(Start_month,End_month+1):
        m='{:02d}'.format(month)
        print(m,y)
        if month in [1,3,5,7,8,10,12]:
            M_days=31
        elif month in [4,6,9,11]:
            M_days=30
        elif month==2:
            if y%4==0:
                M_days=29
            else:
                M_days=28

        #print('starting empty monthly arrays')
        #Month_time=np.ones((M_days*Num_hours,))*np.nan
        #Month_WBGT=np.ones((M_days*Num_hours,Num_lat,Num_lon),dtype=np.float32)*np.nan
        #Read in input variables
        f='{}/{}/{}.{}{}.nc'.format(input_path,WS2m_var,WS2m_var,y,m)
        ws2m=xr.load_dataset(f)

        f='{}/{}/{}.{}{}.nc'.format(input_path,d2m_var,d2m_var,y,m)
        Td=xr.load_dataset(f)
        
        f='{}/{}/{}.{}{}.nc'.format(input_path,t2m_var,t2m_var,y,m)
        T=xr.load_dataset(f)

        f='{}/{}/{}.{}{}.nc'.format(input_path,Srad_var,Srad_var,y,m)
        inSrad=xr.load_dataset(f)

        f='{}/{}/{}.{}{}.nc'.format(input_path,pres_var,pres_var,y,m)
        pres=xr.load_dataset(f)
        pres[pres_var]=(('time','latitude','longitude'),pres[pres_var].values)

        #Get lat and lon
        Grid_lat=T.coords['latitude'].values
        Grid_lon=T.coords['longitude'].values

        #Bias corrections - This is probably not needed
        print('Correcting SRAD')
        CorrectSRAD==False
        CorrectWS2M==False
        if CorrectSRAD: #This section needs work to integrate it
            biasfile="{}/biasinfo.{}.ltm.nc".format(biaspath,m)
#            BiasVars=['SRADbias']
            SRADbias_grid=np.ones(inSrad.shape,dtype=np.float32)
            Varn, Varnm, lat, lon= ReadNetCDF(biasfile,BiasVars)
            SRADbias=Varn[0].astype(np.float32)            
 
            #print(SRADbias)
            if SRADavgb:
                SRADmb=np.array([np.nanmean(SRAD[i::24]) for i in range(24)])
                #print(SRADmb)
                SRADmb_month=np.ones(SRAD.shape)*np.nan
                for i in range(24):
                    SRADmb_month[i::24]=SRADmb[i]
                    #print(SRADmb_month)
            else: SRADmb_month=SRADbias
            
            for hour in range(SRADbias_grid.shape[0]):
                SRADbias_grid[hour,::,::]=SRADmb_month[hour]            

            inSrad=inSrad-SRADbias_grid
        inSrad[Srad_var]=xr.where(inSrad[Srad_var]<0,0.0,inSrad[Srad_var])

        if CorrectWS2M:
            ws2m=ws2m*1.25

        #Creating 2d lat/lon grid for sun angle calculations
        lon2d,lat2d=np.meshgrid(Grid_lon,Grid_lat)
        lat2d=lat2d*np.pi/180
        lon2d=lon2d*np.pi/180

        #Creating date array for sun angle calculations
        date=xr.DataArray(T.coords['time'].values,dims=('time'),coords={'time':T.coords['time']}).chunk({'time':num_chunks})
        print(type(date['time'].data[0]))

        #Sun angle calculations
        cza=da.map_blocks(cosza,date.data,lat2d,lon2d,1,dtype=np.float32,chunks=(num_chunks,lat2d.shape[0],lat2d.shape[1]),new_axis=[1,2])
        cza=xr.DataArray(cza,dims=('time','latitude','longitude'),coords={'time':T.coords['time'].values,'latitude':T.coords['latitude'].values,'longitude':T.coords['longitude'].values})
        print(cza[:,:,:])
        #print(type(cza))
        #cza=np.array(cza,dtype=np.float32)


        print(date.data,lat2d,lon2d)
        exit()

        czda=da.map_blocks(coszda,date.data,lat2d,lon2d,1,dtype=np.float32,chunks=(num_chunks,lat2d.shape[0],lat2d.shape[1]),new_axis=[1,2])
        print(type(cza))
        czda=xr.DataArray(czda,dims=('time','latitude','longitude'),coords={'time':T.coords['time'].values,'latitude':T.coords['latitude'].values,'longitude':T.coords['longitude'].values})
        czda=xr.where(czda<=0,-0.5,czda)

        #Convert Td to RH
        #RH=xr.apply_ufunc(Convert_moisture,Td[d2m_var].values,'Td','K','RH','Fraction',T[t2m_var].values,'K',dask="parallelized",output_dtypes=[np.float32])
        RH=xr.apply_ufunc(Convert_moisture,Td[d2m_var].values,'Td','K','RH','Fraction',T[t2m_var].values,'K',dask="parallelized",output_dtypes=[np.float32])
        #print(RH)
        #RH=RH.compute(num_workers=max_workers)
        RH=xr.DataArray(RH*100,dims=('time','latitude','longitude'),coords={'time':T.coords['time'].values,'latitude':T.coords['latitude'].values,'longitude':T.coords['longitude'].values})
        RH=xr.where(RH>100,100,RH)

        print(np.nanmean(RH[6::24,110,525]))
        print(np.nanmean(RH[18::24,110,525]))
        print(np.nanmean(inSrad[Srad_var].values[6::24,110,525]))
        print(np.nanmean(T[t2m_var].values[6::24,110,525]))
        print(np.nanmean(pres[pres_var].values[6::24,110,525]))
        print(np.nanmean(ws2m[WS2m_var].values[6::24,110,525]))
        print(np.nanmean(inSrad[Srad_var].values[18::24,110,525]))
        print(np.nanmean(T[t2m_var].values[18::24,110,525]))
        print(np.nanmean(pres[pres_var].values[18::24,110,525]))
        print(np.nanmean(ws2m[WS2m_var].values[18::24,110,525]))
        print(T.coords['latitude'].values[110])
        print(T.coords['longitude'].values[525]-180)
        #exit()

        #Get fraction of solar radiation information
        f=xr.apply_ufunc(fdir,cza,czda,inSrad,date,dask="parallelized",output_dtypes=[np.float32])
        #f=xr.apply_ufunc(fdir,cza,czda,inSrad,date,dask="allowed",output_dtypes=[np.float32])

        f=f.compute(num_workers=max_workers)
        cza=cza.compute(num_workers=max_workers)
        f=xr.where(cza<=np.cos(89.5/180*np.pi),0,f)
        f=xr.where(f>0.9,0.9,f)
        f=xr.where(f<0,0,f)
        f=xr.where(inSrad<=0,0,f)

        cza=xr.where(cza<=0,-np.cos(89.5/180*np.pi),cza)

        print(np.nanmax(cza),np.nanmin(cza))
        exit()



        #Calculate WBGT
        wbgt_liljegren=xr.apply_ufunc(WBGT_Liljegren,T[t2m_var].values,RH,pres[pres_var].values,ws2m[WS2m_var].values,inSrad[Srad_var].values,f,czda,WS2MTF,dask="parallelized",output_dtypes=[np.float32])
        print('--------------------------------------')

        #Set up varnames and coordinates for the NetCDF
        wbgt_liljegren=wbgt_liljegren.rename_vars({Srad_var:Var_name})
        wbgt_liljegren.coords['latitude']=T.coords['latitude']
        wbgt_liljegren.coords['longitude']=T.coords['longitude']
        #filename='{}/{}-Liljegren/{}.{}{}.nc'.format(output_path,Var_name,Var_name,y,m)
        filename='{}/{}-Liljegren/{}.{}{}.nc'.format(output_path,Var_name,Var_name,y,m)
        #Set up attributes
        coord_attrs={}
        for dim in T[t2m_var].dims:
            coord_attrs[dim]=T[t2m_var].coords[dim].attrs
        attrs=T[t2m_var].attrs
        attrs['long_name']='Wet Bulb Globe Temperature'
        attrs['units']=Var_unit

        #####################
        ### Write to file ###
        #####################
        WriteNetCDF(wbgt_liljegren,'WBGT',-999,attrs,coord_attrs,filename)

end_datetime=datetime.now()
print(end_datetime-start_datetime)                         

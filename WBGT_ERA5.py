#######################
# Need to make sure changes for chunks are good


import numpy as np
from Python_mods import *
from WBGT_func import *
from netCDF_mods import *
import glob
#import dask
#import dask.array as da
from datetime import datetime

Start_year=1960
End_year=2020
Start_month=1
End_month=12
Num_hours=24
Num_lat=361
Num_lon=720
Chunk_size=8

input_path='/data/deluge/scratch/ERA5-Ben/2D/hourly'
WS2m_var='ws2m'
t2m_var='t2m'
d2m_var='d2m'
Srad_var='msdwswrf'
pres_var='sp'

CorrectSRAD=False
SRADavgb=False #Set to false if using linear fit
CorrectWS2M=False
use_dask_test=False

biaspath="/data/deluge/scratch/ERA5-Ben/2D/hourly/bias"
BiasVars=['SRADbias']

Var_name='WBGT'
Var_unit="F"
Var_long_name="Wet Bulb Globe Temperature"
output_path='/data/deluge/scratch/ERA5-Ben/2D/hourly'

start_datetime=datetime.now()
for y in range(Start_year,End_year+1):
    if y==1959:
        Start_month=2
    elif y==Start_year:
        Start_month=Start_month
    else:
        Start_month=1
    for month in range(Start_month,End_month+1):
        m='{:02d}'.format(month)

        if month in [1,3,5,7,8,10,12]:
            M_days=31
        elif month in [4,6,9,11]:
            M_days=30
        elif month==2:
            if y%4==0:
                M_days=29
            else:
                M_days=28

        print('starting empty monthly arrays')
        Month_time=np.ones((M_days*Num_hours,))*np.nan
        Month_WBGT=np.ones((M_days*Num_hours,Num_lat,Num_lon),dtype=np.float32)*np.nan

        #for c in range(5):
        #    print('starting month', m, 'chunk', str(c+1)+'/5')
        #    if c==0:
        #        chunk_start=0
        #        chunk_end=6*Num_hours
        #    elif c==1:
        #        chunk_start=6*Num_hours
        #        chunk_end=12*Num_hours
        #    elif c==2:
        #        chunk_start=12*Num_hours
        #        chunk_end=18*Num_hours
        #    elif c==3:
        #        chunk_start=18*Num_hours
        #        chunk_end=24*Num_hours
        #    elif c==4:
        #        chunk_start=24*Num_hours
        #        chunk_end=M_days*Num_hours
        #    C_hours=chunk_end-chunk_start
        #    C_Data=np.ones((C_hours,Num_lat,Num_lon))*np.nan

        #Read in input variables components
        #Adjust to only read in for the chunk
        f='{}/{}/{}.{}{}.nc'.format(input_path,WS2m_var,WS2m_var,y,m)
        Varn, Varnm, lat, lon = ReadNetCDF(f,['time',WS2m_var])#,chunk_start,chunk_end)
        time=Varn[0]
        ws2m=Varn[1].astype(np.float32)
#        print(time.shape[0],time.shape[0]/4)
#        exit()

        f='{}/{}/{}.{}{}.nc'.format(input_path,d2m_var,d2m_var,y,m)
        Varn, Varnm, lat, lon = ReadNetCDF(f,['time',d2m_var])#,chunk_start,chunk_end)
        Td=Varn[1].astype(np.float32)
        Td=Convert_T_Unit(Td,'K','C')

        f='{}/{}/{}.{}{}.nc'.format(input_path,t2m_var,t2m_var,y,m)
        Varn, Varnm, lat, lon = ReadNetCDF(f,['time',t2m_var])#,chunk_start,chunk_end)
        T=Varn[1].astype(np.float32)
        T=Convert_T_Unit(T,'K','C')

        f='{}/{}/{}.{}{}.nc'.format(input_path,Srad_var,Srad_var,y,m)
        Varn, Varnm, lat, lon = ReadNetCDF(f,['time',Srad_var])#,chunk_start,chunk_end)
        inSrad=Varn[1].astype(np.float32)

        f='{}/{}/{}.{}{}.nc'.format(input_path,pres_var,pres_var,y,m)
        Varn, Varnm, lat, lon = ReadNetCDF(f,['time',pres_var])#,chunk_start,chunk_end)
        pres=Varn[1].astype(np.float32)/100

        Grid_lat=lat
        Grid_lon=lon
        #Bias corrections
        #This is no longer necessary
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
        inSrad[np.where(inSrad<0)]=0.0

        #This is no longer necessary
        if CorrectWS2M:
            ws2m=ws2m*1.25

        #Do not use. It does not work
        if use_dask_test:
            print('*********************')
            print('WARNING: EXPERIMENTAL')
            print('*********************')
            Chunk_size=24
            Data_size=24

            for h in range(int(time.shape[0]/Data_size)):
                hstarttime=datetime.now()
                start=h*Data_size#*Chunk_size*12
                end=(h+1)*Data_size#*Chunk_size*12
                if end > time.shape[0]:
                    end=time.shape[0]
                print('Hour {}/{} for {}/{}'.format(h*Data_size,time.shape[0],m,y))
                Twbgt=da.map_blocks(Calc_WBGT_from_obs,T[start:end],Td[start:end],inSrad[start:end],ws2m[start:end],pres[start:end],'Td',chunks=(Chunk_size,T.shape[1],T.shape[2]),new_axis=[0,1,2])
                Month_WBGT[start:end,::,::]=np.array(Twbgt,dtype=np.float32)
                hendtime=datetime.now()
                #print(hendtime-hstarttime)
                print(((hendtime-hstarttime)/Data_size)*24*30)
                exit()
        
        else:
            for h in range(int(time.shape[0]/Chunk_size)):
                hstarttime=datetime.now()
                start=h*Chunk_size
                end=(h+1)*Chunk_size
                if end > time.shape[0]:
                    end=time.shape[0]
                print('Hour {}/{} for {}/{}'.format(h*Chunk_size,time.shape[0],m,y))
                Twbgt=Calc_WBGT_from_obs(T[start:end],Td[start:end],inSrad[start:end],ws2m[start:end],pres[start:end],'Td')
                Month_WBGT[start:end,::,::]=Twbgt
                hendtime=datetime.now()
                print(((hendtime-hstarttime)/Chunk_size)*24*30)


        #####################
        ### Write to file ###
        #####################
        lat=Grid_lat[:]
        lon=Grid_lon[:]
        time=time[:]
        Dims=[lon,lat,time]
        DimsName=['longitude','latitude','time']
        DimsUnits=['degrees_east','degrees_north','hours since 1900-01-01 00:00:00.0']
        DimsLN=['longitude','latitude','time']

        Vars=[Month_WBGT]
        VarsNames=[Var_name]
        VarsUnits=[Var_unit]
        VarsLN=[Var_long_name]
        VarsDims=[('time','latitude','longitude')]
        #output_dir=output_path#'{}/OKMesonet.{}{}.nc'.format(output_path,y,m)
        filename='{}-uncorrected-SRADnew/{}.{}{}.nc'.format(Var_name,Var_name,y,m)

        WriteNetCDF(Dims,DimsName,DimsUnits,DimsLN,Vars,VarsNames,VarsUnits,VarsLN,VarsDims,output_path,filename)
end_datetime=datetime.now()
print(end_datetime-start_datetime)                         

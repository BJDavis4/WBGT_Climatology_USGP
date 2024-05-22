import numpy as np
from datetime import datetime
from scipy import interpolate
from read_mdf import *
from WBGT_func import Calc_WBGT_from_obs as Calc_WBGT
from netCDF_mods import WriteNetCDF
import glob


#mdf_files=["./mesonet/202109202100.mdf"]
year_start=2018
year_end=2020
month_start=1
month_end=12

Num_hours=24
Start_lat=32.5
End_lat=38.0
lat_inc=0.5
Num_lat=int((End_lat-Start_lat)/lat_inc)

Start_lon=-104.5
End_lon=-93.5
lon_inc=0.5
Num_lon=int((End_lon-Start_lon)/lon_inc)


input_path='/data/deluge/observations/OKMesonet/latlon'
output_path='/data/deluge/observations/OKMesonet/grid'


#geofile="./mesonet/geoinfo.csv"
#geoinfo=np.loadtxt(geofile,dtype='str',delimiter=',')
#StID_geo=geoinfo[::,1]
#Lat_geo=geoinfo[::,7]
#Lon_geo=geoinfo[::,8]

#need to rerun with these parameters
yi=np.arange(Start_lat,End_lat,lat_inc)#[::-1]
xi=np.arange(Start_lon,End_lon,lon_inc)
Xi,Yi=np.meshgrid(xi,yi)



print('starting loop')
for y in range(year_start,year_end+1):
    #if y == 2017:
    #    month_start=5
    #else: month_start=1
    for month in range(month_start,month_end+1):
        m='{:02d}'.format(month)
        input_dir='{}/{}/{}'.format(input_path,y,m)
        print(input_dir)
        mdf_files=sorted(glob.glob(input_dir+'/*'))
        print(mdf_files)

        if month in [1,3,5,7,8,10,12]:
            M_days=31
        elif month in [4,6,9,11]:
            M_days=30
        elif month==2:
            M_days=29
        M_hours=M_days*Num_hours
        Month_time=np.ones((M_hours,))*np.nan
        Month_T=np.ones((M_hours,Num_lat,Num_lon))*np.nan
        Month_RH=np.ones((M_hours,Num_lat,Num_lon))*np.nan
        Month_inSrad=np.ones((M_hours,Num_lat,Num_lon))*np.nan
        Month_WS_2m=np.ones((M_hours,Num_lat,Num_lon))*np.nan
        Month_p=np.ones((M_hours,Num_lat,Num_lon))*np.nan
        Month_WBGT=np.ones((M_hours,Num_lat,Num_lon))*np.nan
        Month_WS10m=np.ones((M_hours,Num_lat,Num_lon))*np.nan

        if (len(mdf_files)!=M_hours)&(len(mdf_files)!=28*Num_hours):
            exit('Possible missing data')

        for i,mdf_file in enumerate(mdf_files):
            file_name=mdf_file.split('/')[-1]
            day=file_name[6:8]
            hour=file_name[8:10]
            start_date=datetime(1900,1,1,0)
            end_date=datetime(int(y),int(m),int(day),int(hour))
            delta=end_date-start_date
            print(delta)
            Month_time[i]=int(delta.total_seconds()/3600)
                
            print('Reading:',file_name)
            header,Var_names,data=read_mdf(mdf_file,"All")
            #missing_val=np.where(data<=-900)
            #data[missing_val]=np.nan
            #print(header)
            StID_mdf=data[::,0]

            ###################################
            ########### Regrid data ###########
            ###################################
            print('Regriding data')

            #read in lat and lon
            latind=np.where(Var_names=='NLAT')[0]
            lats_mdf=data[::,latind].astype('float')
            lats_mdf.shape=(lats_mdf.shape[0],)

            lonind=np.where(Var_names=='ELON')[0]
            lons_mdf=data[::,lonind].astype('float')
            lons_mdf.shape=(lons_mdf.shape[0],)

            #read in and interpolate individual variables
            Tind=np.where(Var_names=='TAIR')[0]
            T=data[::,Tind].astype('float')
            T.shape=(T.shape[0],)
            nonnan=~np.isnan(T)
            isanan=np.isnan(T)
            try:
                T_grid=interpolate.griddata((lons_mdf[nonnan],lats_mdf[nonnan]),T[nonnan],(Xi,Yi))
                Month_T[i]=T_grid
            except:
                Month_T[i][:,:]=np.nan
                print('exception in T')

            RHind=np.where(Var_names=='RELH')[0]
            RH=data[::,RHind].astype('float')
            RH.shape=(RH.shape[0],)
            nonnan=~np.isnan(RH)
            isanan=np.isnan(RH)
            try:
                RH_grid=interpolate.griddata((lons_mdf[nonnan],lats_mdf[nonnan]),RH[nonnan],(Xi,Yi))
                Month_RH[i]=RH_grid
            except:
                Month_RH[i][:,:]=np.nan
                print('exception in RH')

            inSradind=np.where(Var_names=='SRAD')[0]
            inSrad=data[::,inSradind].astype('float')
            inSrad.shape=(inSrad.shape[0],)
            nonnan=~np.isnan(inSrad)
            isanan=np.isnan(inSrad)
            try:
                inSrad_grid=interpolate.griddata((lons_mdf[nonnan],lats_mdf[nonnan]),inSrad[nonnan],(Xi,Yi))
                Month_inSrad[i]=inSrad_grid
            except:
                Month_inSrad[i][:,:]=np.nan
                print('exception in Srad')

            WS_2mind=np.where(Var_names=='WS2M')[0]
            WS_2m=data[::,WS_2mind].astype('float')
            WS_2m.shape=(WS_2m.shape[0],)
            nonnan=~np.isnan(WS_2m)
            isanan=np.isnan(WS_2m)
            try:
                WS_2m_grid=interpolate.griddata((lons_mdf[nonnan],lats_mdf[nonnan]),WS_2m[nonnan],(Xi,Yi))
                Month_WS_2m[i]=WS_2m_grid
            except:
                Month_WS_2m[i][:,:]=np.nan
                print('exception in WS_2m')

            pind=np.where(Var_names=='PRES')[0]
            p=data[::,pind].astype('float')
            p.shape=(p.shape[0],)
            nonnan=~np.isnan(p)
            isanan=np.isnan(p)
            try:
                p_grid=interpolate.griddata((lons_mdf[nonnan],lats_mdf[nonnan]),p[nonnan],(Xi,Yi))
                Month_p[i]=p_grid
            except:
                Month_p[i][:,:]=np.nan
                print('exception in p')

            wbgtind=np.where(Var_names=='TWBG')[0]
            WBGT=data[::,wbgtind].astype('float')
            WBGT.shape=(WBGT.shape[0],)
            nonnan=~np.isnan(WBGT)
            isanan=np.isnan(WBGT)
            try:
                WBGT_grid=interpolate.griddata((lons_mdf[nonnan],lats_mdf[nonnan]),WBGT[nonnan],(Xi,Yi))
                Month_WBGT[i]=WBGT_grid
            except:
                Month_WBGT[i][:,:]=np.nan
                print('exception in WBGT')

            WS10mind=np.where(Var_names=='WSPD')[0]
            WS10m=data[::,WS10mind].astype('float')
            WS10m.shape=(WS10m.shape[0],)
            nonnan=~np.isnan(WS10m)
            isanan=np.isnan(WS10m)
            try:
                WS10m_grid=interpolate.griddata((lons_mdf[nonnan],lats_mdf[nonnan]),WS10m[nonnan],(Xi,Yi))
                Month_WS10m[i]=WS10m_grid
            except:
                Month_WS10m[i][:,:]=np.nan
                print('exception in WS10m')



            ##################################
            ########## Prepare data ########## need to modify to prepare for writing to netcdf
            ##################################
            #print('finding lat/lon')
        
            #mdf_lats=np.ones(StID_mdf.shape)*np.nan
            #mdf_lons=np.ones(StID_mdf.shape)*np.nan
        
            #for i in range(StID_mdf.shape[0]):
            #    StID=StID_mdf[i]
            #    StID_geo_ind=np.where(StID_geo==StID)[0]
            #    mdf_lats[i]=Lat_geo[StID_geo_ind]
            #    mdf_lons[i]=Lon_geo[StID_geo_ind]

        lat=yi[:]
        lon=xi[:]

        Dims=[lon,lat,Month_time]
        DimsName=['longitude','latitude','time']
        DimsUnits=['degrees_east','degrees_north','hours since 1900-01-01 00:00:00.0']
        DimsLN=['longitude','latitude','time']


        Vars=[Month_T,Month_RH,Month_inSrad,Month_WS_2m,Month_p,Month_WBGT,Month_WS10m]
        VarsNames=['TAIR','RELH','SRAD','WS2M','PRES','TWBG','WSPD']
        VarsUnits=['C','RH%','W m**-2','m s**-2','hPa','F','m s**-2']
        VarsLN=['temperature','relative humidity','incoming solar radiation','2 meter wind speed','surface pressure','wet bulb globe temperature','10 meter wind speed']
        VarsDims=[('time','latitude','longitude'),('time','latitude','longitude'),('time','latitude','longitude'),('time','latitude','longitude'),('time','latitude','longitude'),('time','latitude','longitude'),('time','latitude','longitude')]
        #output_dir=output_path#'{}/OKMesonet.{}{}.nc'.format(output_path,y,m)
        filename='OKMesonet.{}{}.nc'.format(y,m)

        WriteNetCDF(Dims,DimsName,DimsUnits,DimsLN,Vars,VarsNames,VarsUnits,VarsLN,VarsDims,output_path,filename)


            ##################################
            ######### Write new file ######### swap out mdf for netcdf
            ##################################
            #print('writing file')
        
            #mdf_lats.shape=(mdf_lats.shape[0],1)
            #mdf_lons.shape=(mdf_lons.shape[0],1)
            #WBGT.shape=(WBGT.shape[0],1)
            #New_data=[WBGT,mdf_lats,mdf_lons]

            #Var_names.shape=(1,Var_names.shape[0])
            #New_vars=np.array([['TWBG','NLAT','ELON']])
        
            #outfile='{}/{}/{}/{}'.format(output_path,y,m,file_name)
        
            #Var_names,data=append_mdf(Var_names,data,New_vars,New_data,True,outfile,header)
        
            #outfile='./mesomod/202109202100.mdf'
            #write_mdf(outfile,header,Var_names,data)
        
            #for i in range(data.shape[0]):
            #    print(data[i,::])

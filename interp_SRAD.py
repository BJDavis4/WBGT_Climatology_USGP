import numpy as np
from netCDF_mods import ReadNetCDF, WriteNetCDF

Cat0=80
Cat1=85
Cat2=88
Cat3=90
Cat4=212

ERA5_path='/data/deluge/scratch/ERA5-Ben/2D/hourly'
Var_name='msdwswrf'

Start_year=1960
End_year=2020
Start_month=1
End_month=12

num_days=366
num_hours=24

min_lat=25
max_lat=50
min_lon=230
max_lon=295

expected_ei=[31*24,60*24,91*24,121*24,152*24,182*24,213*24,244*24,274*24,305*24,335*24,366*24]
num_lats=361#((max_lat-min_lat)*2)+1
num_lons=720#((max_lon-min_lon)*2)+1
for y in range(Start_year,End_year+1):
    WBGT_year=np.ones((num_days*num_hours,num_lats,num_lons),dtype=np.float32)*np.nan
    time_year=np.ones((num_days*num_hours,))*np.nan
    start_ind=0
    f='{}/{}_native/{}.{}{}.nc'.format(ERA5_path,Var_name,Var_name,y,'{:02d}'.format(1))
    Varn, Varnm, lat, lon = ReadNetCDF(f,['time',Var_name])
    time=Varn[0]
    WBGT_globe=Varn[1].astype(np.float32)
    lat=lat[:]
    lon=lon[:]
    del Varn
    del Varnm
    for month in range(Start_month,End_month+1):
        m='{:02d}'.format(month)
        mp1='{:02d}'.format(month+1)
        if month in [1,3,5,7,8,10,12]:
            M_days=31
        elif month in [4,6,9,11]:
            M_days=30
        else: M_days=29

        M_hours=M_days*num_hours

        if month==End_month:
            f='{}/{}_native/{}.{}{}.nc'.format(ERA5_path,Var_name,Var_name,str(y+1),'{:02d}'.format(1))
        else:
            f='{}/{}_native/{}.{}{}.nc'.format(ERA5_path,Var_name,Var_name,y,mp1)
        Varn, Varnm, lat1, lon1 = ReadNetCDF(f,['time',Var_name])#,chunk_start,chunk_end)
        time1=Varn[0]
        WBGT_globe1=Varn[1].astype(np.float32)
        del Varn
        del Varnm
        WBGT_globe=np.append(WBGT_globe,[WBGT_globe1[0,::,::]],axis=0)
        WBGT_globe=(WBGT_globe[:-1:,::,::]+WBGT_globe[1::,::,::])*0.5
       
        #lons,lats=np.meshgrid(lon,lat)

        #print(lat.shape)
        #print(lon.shape)
        #print(WBGT_globe.shape)
        #trim to desired domain

        ##US_lats=np.where((lat<=max_lat)&(lat>=min_lat))[0]
        ##US_WBGT=WBGT_globe[::,US_lats]
        ##lats=lats[US_lats]
        ##lons=lons[US_lats]
        #print(US_WBGT.shape)
        ##US_lons=np.where((lon<=max_lon)&(lon>=min_lon))[0]
        ##WBGT=US_WBGT[::,::,US_lons]
        ##lats=lats[::,US_lons]
        ##lons=lons[::,US_lons]
        #print(WBGT.shape)
        #print(lats.shape,lons.shape)

        #print(lons)

        #Categorize data
#        WBGT[np.where(WBGT<Cat0)]=0
#        WBGT[np.where((WBGT>=Cat0)&(WBGT<Cat1))]=1
#        WBGT[np.where((WBGT>=Cat1)&(WBGT<Cat2))]=2
#        WBGT[np.where((WBGT>=Cat2)&(WBGT<Cat3))]=3
#        WBGT[np.where((WBGT>=Cat3)&(WBGT<Cat4))]=4

        #print(time.shape[0])

        #either 1) write to file or 2) put in annual array
        #Somehow there is a day missing when adding data to the array by counting
        ##end_ind=start_ind+time.shape[0]
        #print('expected end',expected_ei[month-1])
        #print(time.shape[0])
        ##print(start_ind,end_ind)
        ##WBGT_year[start_ind:end_ind,::,::]=WBGT_globe#WBGT
        ##time_year[start_ind:end_ind]=time
        ##start_ind=start_ind+M_hours
        
        ##print(WBGT_year.shape)
        ##print(time_year)

        ##time=time1
        ##WBGT_globe=WBGT_globe1

        #write to file
        #####################
        ### Write to file ###
        #####################
        #print(lats[:,0],lons[0])
        lat=lat[:]
        lon=lon[:]
        print(lat)
        print(lon)
        time=time[:]
        Dims=[lon,lat,time]
        DimsName=['longitude','latitude','time']
        DimsUnits=['degrees_east','degrees_north','hours since 1900-01-01 00:00:00.0']
        DimsLN=['longitude','latitude','time']

        Vars=[WBGT_globe]
        VarsNames=[Var_name]
        VarsUnits=['W m**-2']
        VarsLN=['Mean surface downward short-wave radiation flux']
        VarsDims=[('time','latitude','longitude')]
            #output_dir=output_path#'{}/OKMesonet.{}{}.nc'.format(output_path,y,m)
        #filename='{}/{}_NA.{}.nc'.format(Var_name,Var_name,y)
        filename='{}/{}.{}{}.nc'.format(Var_name,Var_name,y,m)

        WriteNetCDF(Dims,DimsName,DimsUnits,DimsLN,Vars,VarsNames,VarsUnits,VarsLN,VarsDims,ERA5_path,filename)

        time=time1.copy()
        WBGT_globe=WBGT_globe1.copy()











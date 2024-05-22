import numpy as np
from netCDF_mods import ReadNetCDF, WriteNetCDF
from Python_mods import Convert_T_Unit

Liljegren=True

R1Cat0=76
R1Cat1=81
R1Cat2=84
R1Cat3=86
#R1Cat4=212

R2Cat0=80
R2Cat1=85
R2Cat2=88
R2Cat3=90
#R2Cat4=212

R3Cat0=82
R3Cat1=87
R3Cat2=90
R3Cat3=92
#R3Cat4=212

if Liljegren:
    R1=84#82
    R2=88#86
else:
    R1=82
    R2=86
#R3=212

ERA5_path='/data/deluge/scratch/ERA5-Ben/2D/hourly'
Var_name='WBGT'
Var_climo='MJJAS90_Max'

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

f='{}/{}-Liljegren/{}Climo.nc'.format(ERA5_path,Var_name,Var_name)
Varn, Varnm, lat90, lon90 = ReadNetCDF(f,['time',Var_climo])#,chunk_start,chunk_end)
time90=Varn[0]
WBGTMax90=Varn[1][0].astype(np.float32)
if Liljegren:
    WBGTMax90=Convert_T_Unit(WBGTMax90,'K','F')
del Varn
del Varnm
lat90=lat90[:]
lon90=lon90[:]

num_lats=((max_lat-min_lat)*2)+1
num_lons=((max_lon-min_lon)*2)+1
for y in range(Start_year,End_year+1):
    WBGT_year=np.ones((num_days*num_hours,num_lats,num_lons),dtype=np.float32)*np.nan
    time_year=np.ones((num_days*num_hours,))*np.nan
    start_ind=0
    for month in range(Start_month,End_month+1):
        m='{:02d}'.format(month)
        if month in [1,3,5,7,8,10,12]:
            M_days=31
        elif month in [4,6,9,11]:
            M_days=30
        else: M_days=29

        M_hours=M_days*num_hours

        f='{}/{}-Liljegren/{}.{}{}.nc'.format(ERA5_path,Var_name,Var_name,y,m)
        Varn, Varnm, lat, lon = ReadNetCDF(f,['time',Var_name])#,chunk_start,chunk_end)
        time=Varn[0]
        WBGT_globe=Varn[1].astype(np.float32)
        del Varn
        del Varnm
        lat=lat[:]
        lon=lon[:]

        lons,lats=np.meshgrid(lon,lat)
        #print(lat.shape)
        #print(lon.shape)
        #print(WBGT_globe.shape)
        #trim to desired domain

        US_lats=np.where((lat<=max_lat)&(lat>=min_lat))[0]
        US_WBGT=WBGT_globe[::,US_lats]
        lats=lats[US_lats]
        lons=lons[US_lats]
        #print(US_WBGT.shape)
        US_lons=np.where((lon<=max_lon)&(lon>=min_lon))[0]
        WBGT=US_WBGT[::,::,US_lons]
        lats=lats[::,US_lons]
        lons=lons[::,US_lons]
        #print(WBGT.shape)
        #print(lats.shape,lons.shape)

        #print(lons)
#        end_ind=start_ind+time.shape[0]
        WBGT90=np.ones(WBGT.shape,dtype=np.float32)*WBGTMax90

        if Liljegren:
            WBGT=Convert_T_Unit(WBGT,'K','F')

        #Categorize data
#        print(np.where((WBGT90<R1)&(WBGT<R1Cat0)))
        WBGT[np.where((WBGT90<R1)&(WBGT<R1Cat0))]=0
        WBGT[np.where((WBGT90<R1)&((WBGT>=R1Cat0)&(WBGT<R1Cat1)))]=1
        WBGT[np.where((WBGT90<R1)&((WBGT>=R1Cat1)&(WBGT<R1Cat2)))]=2
        WBGT[np.where((WBGT90<R1)&((WBGT>=R1Cat2)&(WBGT<R1Cat3)))]=3
        WBGT[np.where((WBGT90<R1)&(WBGT>=R1Cat3))]=4

        WBGT[np.where(((WBGT90>=R1)&(WBGT90<R2))&(WBGT<R2Cat0))]=0
        WBGT[np.where(((WBGT90>=R1)&(WBGT90<R2))&((WBGT>=R2Cat0)&(WBGT<R2Cat1)))]=1
        WBGT[np.where(((WBGT90>=R1)&(WBGT90<R2))&((WBGT>=R2Cat1)&(WBGT<R2Cat2)))]=2
        WBGT[np.where(((WBGT90>=R1)&(WBGT90<R2))&((WBGT>=R2Cat2)&(WBGT<R2Cat3)))]=3
        WBGT[np.where(((WBGT90>=R1)&(WBGT90<R2))&(WBGT>=R2Cat3))]=4

        WBGT[np.where((WBGT90>=R2)&(WBGT<R3Cat0))]=0
        WBGT[np.where((WBGT90>=R2)&((WBGT>=R3Cat0)&(WBGT<R3Cat1)))]=1
        WBGT[np.where((WBGT90>=R2)&((WBGT>=R3Cat1)&(WBGT<R3Cat2)))]=2
        WBGT[np.where((WBGT90>=R2)&((WBGT>=R3Cat2)&(WBGT<R3Cat3)))]=3
        WBGT[np.where((WBGT90>=R2)&(WBGT>=R3Cat3))]=4

        #print(time.shape[0])

        #either 1) write to file or 2) put in annual array
        #Somehow there is a day missing when adding data to the array by counting
        end_ind=start_ind+time.shape[0]
        #print('expected end',expected_ei[month-1])
        #print(time.shape[0])
        print(start_ind,end_ind)
        WBGT_year[start_ind:end_ind,::,::]=WBGT
        time_year[start_ind:end_ind]=time
        start_ind=start_ind+M_hours
        
        print(WBGT_year.shape)
        print(time_year)

    #write to file
    #####################
    ### Write to file ###
    #####################
    print(lats[:,0],lons[0])
    lat=lats[:,0]
    lon=lons[0]
    print(lat)
    print(lon)
    time=time_year[:]
    Dims=[lon,lat,time]
    DimsName=['longitude','latitude','time']
    DimsUnits=['degrees_east','degrees_north','hours since 1900-01-01 00:00:00.0']
    DimsLN=['longitude','latitude','time']

    Vars=[WBGT_year]
    VarsNames=['WBGT_cat']
    VarsUnits=['Category']
    VarsLN=['WBGT heat risk category']
    VarsDims=[('time','latitude','longitude')]
        #output_dir=output_path#'{}/OKMesonet.{}{}.nc'.format(output_path,y,m)
    filename='{}-Liljegren/{}Modrisk.{}.nc'.format(Var_name,Var_name,y)

    WriteNetCDF(Dims,DimsName,DimsUnits,DimsLN,Vars,VarsNames,VarsUnits,VarsLN,VarsDims,ERA5_path,filename)

import numpy as np
from netCDF4 import Dataset
import glob
from netCDF_mods import ReadNetCDF, WriteNetCDF
import matplotlib.pyplot as plt
#from calc_z0 import Calc_z0

Input_path="/data/deluge/observations/OKMesonet/grid"
Output_path="/data/deluge/observations/OKMesonet"
Num_years=25
Num_hours=24
max_lat=37.0
min_lat=33.5
max_lon=-97.0
min_lon=-103.0
Num_lat=361
Num_lon=720
Num_lat_meso=11
Num_lon_meso=22
#Num_lat=Num_lat_meso
#Num_lon=Num_lon_meso
save=True

def Calc_z0(ws10m,ws2m):
    ws_ratio=ws2m/ws10m
    ln2=np.log(2)
    ln10=np.log(10)
    lnz0=(ln2-(ws_ratio*ln10))/(1-ws_ratio)
    z0=np.exp(lnz0)
    return z0

for month in range(1,13):
    m='{:02d}'.format(month)
    print('Beginning month: {}'.format(m))
    files=sorted(glob.glob('{}/OKMesonet.*{}.nc'.format(Input_path,m)))
    if month in [1,3,5,7,8,10,12]:
        M_days=31
    elif month in [4,6,9,11]:
        M_days=30
    elif month==2:
        M_days=29

    WS10m=np.ones((Num_years,Num_hours*M_days,Num_lat_meso,Num_lon_meso))*np.nan
    WS2m=np.ones(WS10m.shape)*np.nan
    z0_mean_grid=np.ones((Num_hours*M_days,Num_lat,Num_lon))*np.nan
    z0_std_grid=np.ones(z0_mean_grid.shape)*np.nan
    
    for y, f in enumerate(files):
        Varn, Varnm, lat, lon = ReadNetCDF(f,['time', 'WSPD', 'WS2M'])
        lat_list=lat[:]
        lon_list=lon[:]
        lon,lat=np.meshgrid(lon_list,lat_list)
        time=Varn[0]
        WS10m_meso=Varn[1]
        WS2m_meso=Varn[2]
        print(WS10m[y].shape,WS10m_meso.shape)
        notGPInds=np.where((lat<min_lat)&(lat>max_lat)&(lon<min_lon)&(lon>max_lon))
        WS10m_meso[notGPInds]=np.nan
        WS2m_meso[notGPInds]=np.nan
        WS10m[y]=WS10m_meso
        WS2m[y]=WS2m_meso

    #This might need to be reorganized
    Y,H,I,J=WS10m.shape
    WS10m.shape=(Y*H,I,J)
    WS2m.shape=(Y*H,I,J)

    WS10m[np.where(WS10m<=0)]=np.nan
    WS2m[np.where(WS2m<=0)]=np.nan
    WS2m[np.where(WS2m>=WS10m)]=np.nan

    z0=Calc_z0(WS10m,WS2m)
#    z0[np.where(z0>2)]=np.nan
    print(z0.shape)
    z0_mean=np.ones((24,))
    z0_std=np.ones((24,))

    for h in range(0,24):
        times=z0[h::24]
        print(times)
        mean=np.nanmean(times)
        std=np.nanstd(times)
        print(mean,std)
        z0_mean[h]=np.nanmean(times)
        z0_std[h]=np.nanstd(times)

    #print(np.arange(0,24).shape,z0_mean.shape)

    plt.figure()
    plt.plot(z0_mean,color='b',label='Mean')
    plt.plot(z0_std,color='r',label='SD')
    plt.ylabel('Roughness length (m)')
    plt.xlabel('Time (UTC)')
    plt.yticks([0,0.1,0.2,0.3,0.4,0.5,0.6])
    plt.ylim([0.0,0.62])
    plt.legend()
    plt.title('z_0 climatology for Month: {}'.format(m))
    plt.savefig('Z0Climo_{}'.format(m))
    plt.show()

#while False:
    for hd in range(z0_mean_grid.shape[0]):
        hour=hd%Num_hours
        z0_mean_grid[hd][:,:]=z0_mean[hour]
        z0_std_grid[hd][:,:]=z0_std[hour]

    print(np.arange(0,360,0.5).shape,np.arange(-90,90.1,0.5).shape,time.shape)
    print(z0_mean_grid.shape)
    Dims=[np.arange(0,360,0.5),np.arange(-90,90.1,0.5),time]
    DimsName=['longitude','latitude','time']
    DimsUnits=['degrees_east','degrees_north','hours since 1900-01-01 00:00:00.0']
    DimsLN=['longitude','latitude','time']
    Vars=[z0_mean_grid]
    VarsNames=['z0']
    VarsUnits=['m']
    VarsLN=['Roughness length']
    VarsDims=[('time','latitude','longitude')]
    output_dir='{}/{}'.format(Output_path,'z0')
    filename='{}.{}.mean.nc'.format('z0',m)

    WriteNetCDF(Dims,DimsName,DimsUnits,DimsLN,Vars,VarsNames,VarsUnits,VarsLN,VarsDims,output_dir,filename)

    print('Saving standard deviation to file')
    Vars=[z0_std_grid]
    filename='{}.{}.std.nc'.format('z0',m)
    WriteNetCDF(Dims,DimsName,DimsUnits,DimsLN,Vars,VarsNames,VarsUnits,VarsLN,VarsDims,output_dir,filename)




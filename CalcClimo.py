import numpy as np
from netCDF4 import Dataset
import glob
from netCDF_mods import ReadNetCDF, WriteNetCDF
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

Input_path="/data/deluge/scratch/ERA5-Ben/2D/hourly"
Var_name="WBGT"
Var_unit="K"
Var_long_name="Wet Bulb Globe Temperature"
Num_years=25
Num_days=366
Num_hours=24
Num_lat=361
Num_lon=720

Num_time=Num_days*Num_hours

#Year_Mean=np.ones((Num_days*Num_hours,Num_lat,Num_lon))
#Year_Std=np.ones((Num_days*Num_hours,Num_lat,Num_lon))
#Year_start=0
#Year_end=0

for month in range(1,13):
    m='{:02d}'.format(month)
    print('Beginning month: {}'.format(m))
    print('{}/{}-Liljegren/{}.*{}.nc'.format(Input_path,Var_name,Var_name,m))
    files=sorted(glob.glob('{}/{}-Liljegren/{}.*{}.nc'.format(Input_path,Var_name,Var_name,m)))
    files=files[-Num_years:]
    print(files)
    #exit()
    if month in [1,3,5,7,8,10,12]:
        M_days=31
    elif month in [4,6,9,11]:
        M_days=30
    elif month==2:
        M_days=29

    print('starting empty monthly arrays')
    Month_time=np.ones((M_days*Num_hours,))*np.nan
    Month_Mean=np.ones((M_days*Num_hours,Num_lat,Num_lon),dtype=np.float32)*np.nan
    Month_Std=np.ones((M_days*Num_hours,Num_lat,Num_lon),dtype=np.float32)*np.nan
    print('finished empty monthly arrays')

    for c in range(5):
        print('starting month', m, 'chunk', str(c+1)+'/5')
        if c==0:
            chunk_start=0
            chunk_end=6*Num_hours
        elif c==1:
            chunk_start=6*Num_hours
            chunk_end=12*Num_hours
        elif c==2:
            chunk_start=12*Num_hours
            chunk_end=18*Num_hours
        elif c==3:
            chunk_start=18*Num_hours
            chunk_end=24*Num_hours
        elif c==4:
            chunk_start=24*Num_hours
            chunk_end=M_days*Num_hours
        C_hours=chunk_end-chunk_start
        C_Data=np.ones((25,C_hours,Num_lat,Num_lon),dtype=np.float32)*np.nan
        #Year_start=Year_end
        #Year_end=Year_end+C_hours

    ####### INSERT CODE TO LOOP THROUGH FILES AND READ AND STORE TO DATA ########
        for y, f in enumerate(files):
            print(f)

            Varn, Varnm, lat, lon = ReadNetCDF(f,['time',Var_name],chunk_start,chunk_end)
            time=Varn[0]
            Var=Varn[1].astype(np.float32)
            del Varn
            print(np.nanmax(Var))
            print(np.nanmin(Var))
            print(Var.shape)

            if month == 2 and c == 4:
                VarHours=Var.shape[0]
                C_Data[y,:VarHours]=Var
            else:
                C_Data[y]=Var
            #T,I,J=Var.shape
            #Var.shape=(1,T,I,J)
            #if y==0:
            #    M_data=Var
            #else:
            #    M_data=np.append(M_data,Var,axis=0)
        #print(y,M_data.shape)
#        print(time)
#        print(Var)
        

        print('calculating mean')
        Mean=np.nanmean(C_Data,axis=0)
        print("C_Data:", C_Data[:,0,52,72])
        print("mean",Mean[0,52,72])
        print(np.nanmax(Mean))
        print(np.nanmin(Mean))
        print(Mean.shape)
        print('calculating std')
        Std=np.nanstd(C_Data,axis=0)
        print(Std.shape)
        print('storing chunk mean')
        Month_Mean[chunk_start:chunk_end,:,:]=Mean
        print('storing chunk std')
        Month_Std[chunk_start:chunk_end,:,:]=Std
        print('storing chunk time')
        Month_time[chunk_start:chunk_end]=time

    print('Saving mean to file')
    lat=lat[:]
    #lat.shape=((1,lat.shape[0]))
    lon=lon[:]
    #lon.shape=((1,lon.shape[0]))
    #Month_time.shape=((1,Month_time.shape[0]))
    Dims=[lon,lat,Month_time]
    DimsName=['longitude','latitude','time']
    DimsUnits=['degrees_east','degrees_north','hours since 1900-01-01 00:00:00.0']
    DimsLN=['longitude','latitude','time']
    Vars=[Month_Mean]
    VarsNames=[Var_name]
    VarsUnits=[Var_unit]
    VarsLN=[Var_long_name]
    VarsDims=[('time','latitude','longitude')]
    output_dir='{}/{}-Liljegren'.format(Input_path,Var_name)
    filename='{}.hourly.{}.ltm'.format(Var_name,m)

    WriteNetCDF(Dims,DimsName,DimsUnits,DimsLN,Vars,VarsNames,VarsUnits,VarsLN,VarsDims,output_dir,filename)

    print('Saving standard deviation to file')
    Vars=[Month_Std]
    filename='{}.hourly.{}.ltstd'.format(Var_name,m)
    WriteNetCDF(Dims,DimsName,DimsUnits,DimsLN,Vars,VarsNames,VarsUnits,VarsLN,VarsDims,output_dir,filename)

exit()

















print('plotting data')

minlat=-85
maxlat=85
minlon=0
maxlon=360

lon,lat=np.meshgrid(lon,lat)
latlines=15
lonlines=30

plt.figure()

plt.subplot(231)
m=Basemap(resolution='i', projection='merc', llcrnrlat=minlat, urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon, lat_ts=10)
m.drawcoastlines(linewidth=0.50)
m.drawcountries(linewidth=0.25)
m.drawstates(linewidth=0.15)
m.drawmeridians(np.arange(0,360,lonlines),labels=[False,False,False,True])
m.drawparallels(np.arange(-90,90,latlines),labels=[True,False,False,False])
#lon,lat=np.meshgrid(lon,lat)
x,y=m(lon,lat)

cmap=plt.get_cmap('RdYlGn')
m.contourf(x,y,Month_Mean[0,::,::],cmap=cmap,extend='both')
plt.title('Mean 0z')
plt.colorbar()

plt.subplot(232)
m=Basemap(resolution='i', projection='merc', llcrnrlat=minlat, urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon, lat_ts=10)
m.drawcoastlines(linewidth=0.50)
m.drawcountries(linewidth=0.25)
m.drawstates(linewidth=0.15)
m.drawmeridians(np.arange(0,360,lonlines),labels=[False,False,False,True])
m.drawparallels(np.arange(-90,90,latlines),labels=[True,False,False,False])
#lon,lat=np.meshgrid(lon,lat)
x,y=m(lon,lat)

cmap=plt.get_cmap('RdYlGn')
m.contourf(x,y,Month_Mean[12,::,::],cmap=cmap,extend='both')
plt.title('Mean 12z')
plt.colorbar()

plt.subplot(233)
m=Basemap(resolution='i', projection='merc', llcrnrlat=minlat, urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon, lat_ts=10)
m.drawcoastlines(linewidth=0.50)
m.drawcountries(linewidth=0.25)
m.drawstates(linewidth=0.15)
m.drawmeridians(np.arange(0,360,lonlines),labels=[False,False,False,True])
m.drawparallels(np.arange(-90,90,latlines),labels=[True,False,False,False])
#lon,lat=np.meshgrid(lon,lat)
x,y=m(lon,lat)

cmap=plt.get_cmap('RdYlGn')
m.contourf(x,y,Month_Mean[12,::,::]-Month_Mean[0,::,::],cmap=cmap,extend='both')
plt.title('12z-0z mean')
plt.colorbar()

plt.subplot(234)
m=Basemap(resolution='i', projection='merc', llcrnrlat=minlat, urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon, lat_ts=10)
m.drawcoastlines(linewidth=0.50)
m.drawcountries(linewidth=0.25)
m.drawstates(linewidth=0.15)
m.drawmeridians(np.arange(0,360,lonlines),labels=[False,False,False,True])
m.drawparallels(np.arange(-90,90,latlines),labels=[True,False,False,False])
#lon,lat=np.meshgrid(lon,lat)
x,y=m(lon,lat)

cmap=plt.get_cmap('RdYlGn')
m.contourf(x,y,Month_Std[0,::,::],cmap=cmap,extend='both')
plt.title('Std 0z')
plt.colorbar()

plt.subplot(235)
m=Basemap(resolution='i', projection='merc', llcrnrlat=minlat, urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon, lat_ts=10)
m.drawcoastlines(linewidth=0.50)
m.drawcountries(linewidth=0.25)
m.drawstates(linewidth=0.15)
m.drawmeridians(np.arange(0,360,lonlines),labels=[False,False,False,True])
m.drawparallels(np.arange(-90,90,latlines),labels=[True,False,False,False])
#lon,lat=np.meshgrid(lon,lat)
x,y=m(lon,lat)

cmap=plt.get_cmap('RdYlGn')
m.contourf(x,y,Month_Std[12,::,::],cmap=cmap,extend='both')
plt.title('Std 12z')
plt.colorbar()

plt.subplot(236)
m=Basemap(resolution='i', projection='merc', llcrnrlat=minlat, urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon, lat_ts=10)
m.drawcoastlines(linewidth=0.50)
m.drawcountries(linewidth=0.25)
m.drawstates(linewidth=0.15)
m.drawmeridians(np.arange(0,360,lonlines),labels=[False,False,False,True])
m.drawparallels(np.arange(-90,90,latlines),labels=[True,False,False,False])
#lon,lat=np.meshgrid(lon,lat)
x,y=m(lon,lat)

cmap=plt.get_cmap('RdYlGn')
m.contourf(x,y,Month_Std[12,::,::]-Month_Std[0,::,::],cmap=cmap,extend='both')
plt.title('12z-0z Std')
plt.colorbar()

plt.savefig('interesting_figure.png')
plt.show()




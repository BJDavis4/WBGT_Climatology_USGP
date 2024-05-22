import numpy as np
from netCDF_mods import ReadNetCDF, WriteNetCDF
from datetime import datetime

ERA5_path='/data/deluge/scratch/ERA5-Ben/2D/hourly'
Var_name='WBGT'

Start_year=1960
End_year=2020
Start_month=1
End_month=12

num_days=366
num_hours=24
num_years=End_year-Start_year+1

min_lat=25
max_lat=50
min_lon=230
max_lon=295

Calc_Daily=True#False
Calc_RM=True
Calc_Monthly=True#False
Calc_Season=True#False

expected_si=[0*24,31*24,60*24,91*24,121*24,152*24,182*24,213*24,244*24,274*24,305*24,335*24]
expected_ei=[31*24,60*24,91*24,121*24,152*24,182*24,213*24,244*24,274*24,305*24,335*24,366*24]
num_lats=((max_lat-min_lat)*2)+1
num_lons=((max_lon-min_lon)*2)+1
WBGT_all=np.ones((num_years,num_days*num_hours,num_lats,num_lons),dtype=np.float32)*np.nan
WBGT_Daily_Max=np.ones((num_years,num_days,num_lats,num_lons),dtype=np.float32)*np.nan
WBGT_Daily_Min=np.ones((num_years,num_days,num_lats,num_lons),dtype=np.float32)*np.nan
WBGT_Daily_Mean=np.ones((num_years,num_days,num_lats,num_lons),dtype=np.float32)*np.nan
for y in range(Start_year,End_year+1):
    #WBGT_year=np.ones((num_days*num_hours,num_lats,num_lons))*np.nan
    time_year=np.ones((num_days*num_hours,))*np.nan
    start_ind=0

    f='{}/{}-Liljegren/{}_NA.{}.nc'.format(ERA5_path,Var_name,Var_name,y)
    Varn, Varnm, lat, lon = ReadNetCDF(f,['time',Var_name])#,chunk_start,chunk_end)
    time=Varn[0]
    WBGT=Varn[1].astype(np.float32)
    del Varn
    del Varnm
    lat=lat[:]
    lon=lon[:]

    lons,lats=np.meshgrid(lon,lat)

    WBGT_all[y-Start_year,::,::,::]=WBGT
    time_year=time

#    for month in range(Start_month,End_month+1):
#        m='{:02d}'.format(month)
#        if month in [1,3,5,7,8,10,12]:
#            M_days=31
#        elif month in [4,6,9,11]:
#            M_days=30
#        else: M_days=29
#
#        M_hours=M_days*num_hours

#        #either 1) write to file or 2) put in annual array
#        #Somehow there is a day missing when adding data to the array by counting
#        end_ind=start_ind+time.shape[0]
#        #print('expected end',expected_ei[month-1])
#        #print(time.shape[0])
#        print(start_ind,end_ind)
#        WBGT_all[y-Start_year,start_ind:end_ind,::,::]=WBGT
#        time_year[start_ind:end_ind]=time
#        start_ind=start_ind+M_hours
#        
#        print(WBGT_all.shape)
#        print(time_year)
    
    for d in range(num_days):
        WBGT_Daily_Max[y-Start_year,d,::,::]=np.nanmax(WBGT_all[y-Start_year,d*num_hours:(d+1)*num_hours,::,::],axis=0)
        WBGT_Daily_Min[y-Start_year,d,::,::]=np.nanmin(WBGT_all[y-Start_year,d*num_hours:(d+1)*num_hours,::,::],axis=0)
        WBGT_Daily_Mean[y-Start_year,d,::,::]=np.nanmean(WBGT_all[y-Start_year,d*num_hours:(d+1)*num_hours,::,::],axis=0)

    print(WBGT_Daily_Max[y-Start_year,::,30,65])
    print(WBGT_Daily_Min[y-Start_year,::,30,65])
    print(WBGT_Daily_Mean[y-Start_year,::,30,65])

print('Finding Daily Percentiles')
if Calc_Daily:
    WBGT90_Max=np.nanpercentile(WBGT_Daily_Max,90,axis=0,interpolation='linear')
    WBGT90_Min=np.nanpercentile(WBGT_Daily_Min,90,axis=0,interpolation='linear')
    WBGT90_Mean=np.nanpercentile(WBGT_Daily_Mean,90,axis=0,interpolation='linear')
else:
    f='{}/{}/{}Climo.nc'.format(ERA5_path,Var_name,Var_name)
    Varn, Varnm, lat, lon = ReadNetCDF(f,['time','WBGT90D','WBGT10D','WBGT50D'])
    WBGT90_Max=Varn[1].astype(np.float32)
    WBGT90_Min=Varn[2].astype(np.float32)
    WBGT90_Mean=Varn[3].astype(np.float32)
    del Varn
    del Varnm

#    f='{}/{}/{}Climo.nc'.format(ERA5_path,Var_name,Var_name)
#    Varn, Varnm, lat, lon = ReadNetCDF(f,['time','WBGT90D','WBGT10D','WBGT50D'])#,chunk_start,chunk_end)
#    time=Varn[0]
#    WBGT90_Max=Varn[1]
#    WBGT90_Min=Varn[2]
#    WBGT90_Mean=Varn[3]

print(WBGT_Daily_Max.shape)
print(WBGT_Daily_Min.shape)
print(WBGT_Daily_Mean.shape)
print(WBGT90_Max[::,30,65])
print(WBGT90_Min[::,30,65])
print(WBGT90_Mean[::,30,65])
print(WBGT90_Max.shape)
print(WBGT90_Min.shape)
print(WBGT90_Mean.shape)


print('Starting 31 day percentiles')
time_start=datetime.now()

if Calc_RM:
    WBGT_Daily_Max_RM=np.ones((num_years,num_days+30,num_lats,num_lons),dtype=np.float32)*np.nan
    WBGT_Daily_Min_RM=np.ones((num_years,num_days+30,num_lats,num_lons),dtype=np.float32)*np.nan
    WBGT_Daily_Mean_RM=np.ones((num_years,num_days+30,num_lats,num_lons),dtype=np.float32)*np.nan

    WBGT90_Max_RM=np.ones(WBGT90_Max.shape,dtype=np.float32)*np.nan
    WBGT_Daily_Max_RM[::,:15:,::,::]=WBGT_Daily_Max[::,-15::,::,::]
    WBGT_Daily_Max_RM[::,15:-15:,::,::]=WBGT_Daily_Max
    WBGT_Daily_Max_RM[::,-15::,::,::]=WBGT_Daily_Max[::,:15:,::,::]

    WBGT90_Min_RM=np.ones(WBGT90_Min.shape,dtype=np.float32)*np.nan
    WBGT_Daily_Min_RM[::,:15:,::,::]=WBGT_Daily_Min[::,-15::,::,::]
    WBGT_Daily_Min_RM[::,15:-15:,::,::]=WBGT_Daily_Min
    WBGT_Daily_Min_RM[::,-15::,::,::]=WBGT_Daily_Min[::,:15:,::,::]

    WBGT90_Mean_RM=np.ones(WBGT90_Mean.shape,dtype=np.float32)*np.nan
    WBGT_Daily_Mean_RM[::,:15:,::,::]=WBGT_Daily_Mean[::,-15::,::,::]
    WBGT_Daily_Mean_RM[::,15:-15:,::,::]=WBGT_Daily_Mean
    WBGT_Daily_Mean_RM[::,-15::,::,::]=WBGT_Daily_Mean[::,:15:,::,::]

    for d in range(num_days):
        print(d)
        WBGT90_Max_RM[d,::,::]=np.nanpercentile(WBGT_Daily_Max_RM[::,d:d+31:,::,::],90,axis=(0,1),interpolation='linear')
        WBGT90_Min_RM[d,::,::]=np.nanpercentile(WBGT_Daily_Min_RM[::,d:d+31:,::,::],90,axis=(0,1),interpolation='linear')
        WBGT90_Mean_RM[d,::,::]=np.nanpercentile(WBGT_Daily_Mean_RM[::,d:d+31:,::,::],90,axis=(0,1),interpolation='linear')

#WBGT_Daily_Max
#WBGT_Daily_Min
#WBGT_Daily_Mean

else:
    f='{}/{}/{}Climo.nc'.format(ERA5_path,Var_name,Var_name)
    Varn, Varnm, lat, lon = ReadNetCDF(f,['time','WBGT90RM','WBGT10RM','WBGT50RM'])
    WBGT90_Max_RM=Varn[1].astype(np.float32)
    WBGT90_Min_RM=Varn[2].astype(np.float32)
    WBGT90_mean_RM=Varn[3].astype(np.float32)
    del Varn
    del Varnm

time_end=datetime.now()
time_total=time_end-time_start
print('Finished with 31 day percentiles in {}'.format(str(time_total)))


WBGT90_Max_monthly=np.ones(WBGT90_Max.shape,dtype=np.float32)*np.nan
WBGT90_Min_monthly=np.ones(WBGT90_Min.shape,dtype=np.float32)*np.nan
WBGT90_Mean_monthly=np.ones(WBGT90_Mean.shape,dtype=np.float32)*np.nan


print('Finding Monthly percentiles')
if Calc_Monthly:
    for m in range(12):
        #WBGTMaxMonth=WBGT_Daily_Max[::,expected_si[m]:expected_ei[m],::,::]
        #WBGTMaxMonth=WBGT_Daily_Min[::,expectes_si[m]:expected_ei[m],::,::]
        #Y,D,I,J=WBGTMaxMonth.shape
        #WBGTMaxMonth.shape=(Y*D,I,J)
        WBGT90_Max_monthly[int(expected_si[m]/24):int(expected_ei[m]/24),::,::]=np.nanpercentile(WBGT_Daily_Max[::,int(expected_si[m]/24):int(expected_ei[m]/24),::,::],90,axis=(0,1))
        WBGT90_Min_monthly[int(expected_si[m]/24):int(expected_ei[m]/24),::,::]=np.nanpercentile(WBGT_Daily_Min[::,int(expected_si[m]/24):int(expected_ei[m]/24),::,::],90,axis=(0,1))
        WBGT90_Mean_monthly[int(expected_si[m]/24):int(expected_ei[m]/24),::,::]=np.nanpercentile(WBGT_Daily_Mean[::,int(expected_si[m]/24):int(expected_ei[m]/24),::,::],90,axis=(0,1))

else:
    f='{}/{}/{}Climo.nc'.format(ERA5_path,Var_name,Var_name)
    Varn, Varnm, lat, lon = ReadNetCDF(f,['time','WBGT90M','WBGT10M','WBGT50M'])
    WBGT90_Max_monthly=Varn[1].astype(np.float32)
    WBGT90_Min_monthly=Varn[2].astype(np.float32)
    WBGT90_Mean_monthly=Varn[3].astype(np.float32)
    del Varn
    del Varnm


print(WBGT90_Max_monthly[::,30,65])
print(WBGT90_Min_monthly[::,30,65])
print(WBGT90_Mean_monthly[::,30,65])

if Calc_Season:
    MJJAS90_Max=np.ones(WBGT90_Max.shape,dtype=np.float32)*np.nan
    MJJAS90_all=np.ones(WBGT90_Min.shape,dtype=np.float32)*np.nan
    JJA90_Max=np.ones(WBGT90_Max.shape,dtype=np.float32)*np.nan
    JJA90_all=np.ones(WBGT90_Min.shape,dtype=np.float32)*np.nan


    MJJAS90_Max[::,::,::]=np.nanpercentile(WBGT_Daily_Max[::,int(expected_si[4]/24):int(expected_ei[8]/24),::,::],90,axis=(0,1))
    JJA90_Max[::,::,::]=np.nanpercentile(WBGT_Daily_Max[::,int(expected_si[5]/24):int(expected_ei[7]/24),::,::],90,axis=(0,1))

    MJJAS90_all[::,::,::]=np.nanpercentile(WBGT_all[::,expected_si[4]:expected_ei[8],::,::],90,axis=(0,1))
    JJA90_all[::,::,::]=np.nanpercentile(WBGT_all[::,expected_si[5]:expected_ei[7],::,::],90,axis=(0,1))
else:
    f='{}/{}/{}Climo.nc'.format(ERA5_path,Var_name,Var_name)
    Varn, Varnm, lat, lon = ReadNetCDF(f,['time','MJJAS90_Max','MJJAS90','JJA90_Max','JJA90'])
    MJJAS90_Max=Varn[1].astype(np.float32)
    MJJAS90_all=Varn[2].astype(np.float32)
    JJA90_Max=Varn[3].astype(np.float32)
    JJA90_all=Varn[4].astype(np.float32)
    del Varn
    del Varnm



for i in [WBGT90_Max,WBGT90_Min,WBGT90_Max_monthly,WBGT90_Min_monthly,MJJAS90_Max,MJJAS90_all,JJA90_Max,JJA90_all]:
    print(i.shape)

print(lat.shape)
print(lon.shape)
print(time_year.shape)




#write to file
#####################
### Write to file ###
#####################
print(lats[:,0],lons[0])
lat=lats[:,0]
lon=lons[0]
print(lat)
print(lon)
time=time_year[::24]
Dims=[lon,lat,time]    
DimsName=['longitude','latitude','time']
DimsUnits=['degrees_east','degrees_north','hours since 1900-01-01 00:00:00.0']
DimsLN=['longitude','latitude','time']

Vars=[WBGT90_Max,WBGT90_Min,WBGT90_Mean,WBGT90_Max_monthly,WBGT90_Min_monthly,WBGT90_Mean_monthly,WBGT90_Max_RM,WBGT90_Min_RM,WBGT90_Mean_RM,MJJAS90_Max,MJJAS90_all,JJA90_Max,JJA90_all]
VarsNames=['WBGT90D','WBGT10D','WBGT50D','WBGT90M','WBGT10M','WBGT50M','WBGT90RM','WBGT10RM','WBGT50RM','MJJAS90_Max','MJJAS90','JJA90_Max','JJA90']
#VarsUnits=['F','F','F','F','F','F','F','F','F','F','F','F','F']
VarsUnits=['K']*len(VarsNames)
print(VarsUnits)
VarsLN=['Daily 90th pct Max WBGT','Daily 90th pct Min WBGT','Daily 90th pct Mean WBGT','Monthly 90th pct Max WBGT','Monthly 90th pct Min WBGT','Monthly 90th pct Mean WBGT','31 day 90th pct Max WBGT','31 day 90th pct Min WBGT','31 day 90th pct Mean WBGT','MJJAS 90th pct Max WBGT','MJJAS 90th pct WBGT','JJA 90th pct Max WBGT','JJA 90th pct WBGT']
VarsDims=[('time','latitude','longitude'),('time','latitude','longitude'),('time','latitude','longitude'),('time','latitude','longitude'),('time','latitude','longitude'),('time','latitude','longitude'),('time','latitude','longitude'),('time','latitude','longitude'),('time','latitude','longitude'),('time','latitude','longitude'),('time','latitude','longitude'),('time','latitude','longitude'),('time','latitude','longitude')]
    #output_dir=output_path#'{}/OKMesonet.{}{}.nc'.format(output_path,y,m)
filename='{}-Liljegren/{}Climo_temp.nc'.format(Var_name,Var_name)

WriteNetCDF(Dims,DimsName,DimsUnits,DimsLN,Vars,VarsNames,VarsUnits,VarsLN,VarsDims,ERA5_path,filename)


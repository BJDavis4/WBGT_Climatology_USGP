import numpy as np
from netCDF_mods import ReadNetCDF, WriteNetCDF
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from numba import jit
from MonteCarlo import SigTrend
from Python_mods import Convert_T_Unit
import xarray as xr
#from LinRegress import LinRegress
#define necessary variables

ERA5_path='/data/deluge/scratch/ERA5-Ben/2D/hourly'
WBGT_type='WBGT'#['Cat','WBGT','ModCat']
Var_name_path='WBGT'
Formula='Mesonet' #[Liljegren,Mesonet]

if Formula=='Liljegren':
    form_name='Liljegren'
elif Formula=='Mesonet':
    form_name='nobc'
if WBGT_type=='Cat':
    Var_name='WBGT_cat'
    Save_dir='figures-{}-WBGTcat-RM'.format(form_name)
#    Save_dir='figures-nobc-WBGTcat-RM'
    Plot_var='Cat'
elif WBGT_type=='WBGT':
    Var_name='WBGT'
    Save_dir='figures-{}-WBGT-RM'.format(form_name)
#    Save_dir='figures-nobc-WBGT-RM'
    Plot_var='WBGT'
elif WBGT_type=='ModCat':
    Var_name='WBGT_cat'
    Save_dir='figures-{}-WBGTmodcat-RM'.format(form_name)
#    Save_dir='figures-nobc-WBGTmodcat-RM'
    Plot_var='Cat'
elif WBGT_type=='input':
    Var_name=Var_name_path
    Save_dir='figures/figures-SRAD'
    Plot_var='SRAD'
Start_year=1960
End_year=2020
Start_month=1
End_month=12

num_days=366
num_months=End_month-Start_month+1
num_hours=24
num_leap_years=16#7
num_lat=51
num_lon=131
num_years=End_year-Start_year+1
num_cats=5
HWMaxBreak=1

region='South'
if region=='South':
    ref_lat=30
    ref_lon=65
elif region=='Central':
    ref_lat=20#8#30
    ref_lon=60#65
elif region=='North':
    ref_lat=8#30
    ref_lon=60#65

trend_scale=10
plot_time=20
Pre_calc=True
plt_bounds='Manual' #['Manual','Auto']
pltsize=(4,5)#(6.4,4.8)
figdpi=1000#1000 #Set DPI to easily switch between pub quality and draft

#Set start (in days) and end index (in hours), and label months (sorry its different units I needed them for different things at first)
expected_ei=[31*24,60*24,91*24,121*24,152*24,182*24,213*24,244*24,274*24,305*24,335*24,366*24]
month_start_days=[0,31,60,91,121,152,182,213,244,274,305,335]
month_labels=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
minlats=[26,14,0]#[37,43,50]
maxlats=[36,26,14]#[32,37,43]
maxlons=[72,72,72]#[S,C,N]
minlons=[48,48,48]#[S,C,N]
Region_Labels=['South','Central','North']

#Set which sets of figures to plot
HourlyClimo=False
MonthlyClimo=True
DailyClimo=False
ShiftedDailyClimo=False
MonthlyFullClimo=False
RegionalMean=False
#PlotHeatDays=True


#@jit
def LinRegress(inputs,y):
    E=np.ones((len(inputs[0]),len(inputs)+1))
    for i in range(len(inputs)):
        E[:,i]=inputs[i]
    #E[:,-1]=np.random.normal(size=E[:,-1].shape[0])
    xhat = np.dot(np.dot(np.linalg.inv(np.dot(np.transpose(E),E)),np.transpose(E)),y)
    return E,xhat

#######Code from general/heat_class.py to identify heat days#######
#######Modify for this data format
#######HWtemp is threshold variables
#######Full_data is the data array
#######Might not need all the reshaping idk
#######Make sure that "default" category (24 hours, 0 int, 0 cat) isnt considered heat day
#@jit
def CalcHeatDays(Full_data,HWtemp,Side_of_thresh='Above',MaxBreak=1,HeatType='Wave'):
    #Full_data --> data array of variable used in determining heat threshold
    #HWtemp --> Threshold values that must be achieved to be a heat day
    #Side_of_thresh --> controls whether a value must be above or below the 
    #                   threshold to be considered a heat day. 
    #MaxBreak --> maximum number of consecutive days between heat days for them to be
    #             considered the same event
    if HeatType not in ['Days','Wave']:
        exit('Invalid HeatType. Valid types:', ['Day','Wave'])

    HW_data=np.ones(Full_data.shape)*np.nan
    for i in range(Full_data.shape[0]):
        HW_data[i,:,:,:]=Full_data[i,:,:,:]-HWtemp

    if Side_of_thresh=='Below':
        HW_data=HW_data*-1

    #Y,T,I,J=Full_data.shape
    print(HW_data.shape)
    HW_data.shape=(Y*num_days,I,J)
    HW_data.shape=(Y*num_days,I*J)

    Year=np.ones((num_years,num_days))
    for i in range(num_days):
        Year[:,i]=np.arange(Start_year,End_year+1)
    Year[np.where(np.arange(Start_year,End_year+1)%4!=0),59]=np.nan
    print(Year[:,59])
    Year.shape=(Y*num_days)
    print(Year)
    #Month.shape=(Y*T)
    #Day.shape=(Y*T)

    HW_data=np.rint(HW_data) #why? I guess to account for rounding errors in categories
    HWinds=np.where(HW_data[0]>=0)
    nonHWinds=np.where(HW_data[0]<0)
    HW_data[0,HWinds]=1
    HW_data[0,nonHWinds]=0

    if MaxBreak==0:
        for i in range(1,Y*num_days):
            if np.isnan(Year[i-1]):
                prev_ind=i-2
                #HWinds=np.where(HW_data[i]>=0)
                #nonHWinds=np.where(HW_data[i]<0)
                #HW_data[i,HWinds]=HW_data[i-2,HWinds]+1
                #HW_data[i,nonHWinds]=0
            else: prev_ind=i-1

            HWinds=np.where(HW_data[i]>=0)
            #nonHWinds=np.where(HW_data[i]<0)

            HW_data[i,HWinds]=HW_data[prev_ind,HWinds]+1
            #HW_data[i,nonHWinds]=0 #Replaced with if/elif below to clear single day waves

            if HeatType=='Days':
                nonHWinds=np.where(HW_data[i]<0)
            elif HeatType=='Wave':
                #print('HeatType=Wave check')
                nonHWinds=np.where((HW_data[i]<0)&(HW_data[prev_ind]!=1))
                nonHWInds1=np.where((HW_data[i]<0)&(HW_data[prev_ind]==1))
                #print(len(nonHWInds1),len(nonHWInds1[0]),print(nonHWInds1))
                HW_data[prev_ind,nonHWInds1]=0
                HW_data[i,nonHWInds1]=0
            #else: print('Did not enter wave')
            HW_data[i,nonHWinds]=0

    elif MaxBreak==1:
        for i in range(1,Y*num_days):
            if np.isnan(Year[i-1]):
                prev_ind=i-2
            else: prev_ind=i-1
            #Note, to require 2+ consecutive days to start a heat wave
            #just adjust the threshold for prev_ind in the nonHWinds to 1
            #However, this will eliminate individual heat days.
            #To keep individual heat days create nonHWinds_one
            HWInds_pos=np.where((HW_data[i]>=0)&(HW_data[prev_ind]>=0))
            HWInds_neg=np.where((HW_data[i]>=0)&(HW_data[prev_ind]<0))
            #nonHWinds_pos=np.where((HW_data[i]<0)&(HW_data[prev_ind]>0))
            nonHWinds_neg=np.where((HW_data[i]<0)&(HW_data[prev_ind]<=0))

            HW_data[i,HWInds_pos]=HW_data[prev_ind,HWInds_pos]+1

            HW_data[prev_ind,HWInds_neg]=HW_data[prev_ind,HWInds_neg]*-1
            HW_data[i,HWInds_neg]=HW_data[prev_ind,HWInds_neg]+1

            #HW_data[i,nonHWinds_pos]=(HW_data[pre_ind,nonHWinds_pos]+1)*-1

            HW_data[i,nonHWinds_neg]=0
            HW_data[prev_ind,nonHWinds_neg]=0

            #Set nonHWinds_pos as case may need and HeatType=='Wave' specific uses
            #Note: For HeatType=='Wave' this logic requires 2+ consucutive heat 
            #      days to start and then 2+ consecutive non heat days to end. 
            #      A single heat day followed by a single non-heat day followed 
            #      by 2+ heat days will start at the 2+ heat days, not the single 
            #      heat day
            if HeatType=='Days':
                nonHWinds_pos=np.where((HW_data[i]<0)&(HW_data[prev_ind]>0))

            elif HeatType=='Wave':
                nonHWinds_pos=np.where((HW_data[i]<0)&(HW_data[prev_ind]>1))
                nonHWinds_one=np.where((HW_data[i]<0)&(HW_data[prev_ind]==1))

                #HW_data[i,nonHWinds_pos]=(HW_data[pre_ind,nonHWinds_pos]+1)*-1
                #This does the same thing for nonHWinds_one as it does nonHWinds_neg
                HW_data[i,nonHWinds_one]=0
                HW_data[prev_ind,nonHWinds_one]=0

            HW_data[i,nonHWinds_pos]=(HW_data[prev_ind,nonHWinds_pos]+1)*-1

    #del HWinds
    #del HWInds_pos
    #del HWInds_neg
    #del noHWinds_one
    #del nonHWinds

    HW_data.shape=(Y*num_days,I,J)
    HW_data.shape=(Y,num_days,I,J)
    return HW_data

#Plot data on a map
def PlotMap(Data,lon_2D,lat_2D,levs,exttype,plttitle,pltsave,colormap='RdYlGn_r',units=None,UseSig=False,Sig=None):
    if plt_bounds=='Auto':
        plt.figure()
    elif plt_bounds=='Manual':
        plt.figure(figsize=(4.5,5))
    m=Basemap(resolution='i', projection='merc', llcrnrlat=minlat, urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon, lat_ts=10)
    try:
        m.drawcoastlines(linewidth=0.50)
    except:
        pass
    m.drawcountries(linewidth=0.25)
    m.drawstates(linewidth=0.15)
    m.drawmeridians(np.arange(0,360,5),labels=[False,False,False,True])
    m.drawparallels(np.arange(-90,90,5),labels=[True,False,False,False])
    x,y=m(lon_2D,lat_2D)
    #xi,yi=m(lon[ref_lat,ref_lon],lat[ref_lat,ref_lon])
    cmap=plt.get_cmap(colormap)
    m.contourf(x,y,Data,cmap=cmap,extend=exttype,levels=levs)
    cbar=plt.colorbar()
    cbar.set_label(units,fontsize=20)
    cbar.ax.tick_params(labelsize=16)
    if UseSig:
        try:
            if Sig==None: SigPts=None
        except: SigPts=np.where(~Sig)
        plt.scatter(x[SigPts],y[SigPts],s=4,c='k')
    #plt.title(plttitle,fontsize=24)
    plt.tight_layout()
    plt.savefig(pltsave)
    plt.show()
    plt.close()



#read in wbgt categories
	#If possible, fit into 25xdata_size array
Cats=np.ones((num_years,num_days*num_hours,num_lat,num_lon),dtype=np.float32)*np.nan
for i,y in enumerate(range(Start_year,End_year+1)):
    print(i,y)
    if WBGT_type=='Cat':
        f='{}/{}/{}risk.{}.nc'.format(ERA5_path,Var_name_path,Var_name_path,y)
    elif WBGT_type=='WBGT':
        f='{}/{}/{}_NA.{}.nc'.format(ERA5_path,Var_name_path,Var_name_path,y)
    elif WBGT_type=='ModCat':
        f='{}/{}/{}Modrisk.{}.nc'.format(ERA5_path,Var_name_path,Var_name_path,y)
    elif WBGT_type=='input':
        f='{}/{}/{}_NA.{}.nc'.format(ERA5_path,Var_name_path,Var_name_path,y)

    #Varn, Varnm, lat, lon = ReadNetCDF(f,['time',Var_name])#,chunk_start,chunk_end)
    Varn=xr.load_dataset(f)
    Varn_units=Varn[Var_name].units
    mv=Varn[Var_name].encoding['missing_value']
    #time=Varn[0]
    time_data=Varn.coords['time'].values
    lat=Varn.coords['latitude'].values
    lon=Varn.coords['longitude'].values

    Varn[Var_name]=xr.where(Varn[Var_name]==mv,np.nan,Varn[Var_name])
    #print(Varn[Var_name].units)
    #print(time_data)
    #print(lat)
    #print(lon)
    #exit()
    if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
        Cats[i]=np.rint(Varn[Var_name].values).astype(np.float32)
    elif (WBGT_type=='WBGT') or (WBGT_type=='input'):
        Cats[i]=Varn[Var_name].values.astype(np.float32)
        Cats[i]=Convert_T_Unit(Cats[i],Varn_units,'F')

    del Varn
    #del Varnm
    lat_1D=lat[:]
    lon_1D=lon[:]
    #exit()

print(np.nanmin(Cats),np.nanmax(Cats))

#lat=np.arange(50,24.9,-0.5)
#lon=np.arange(230,295.1,0.5)
print(lat_1D,lat_1D.shape)
print(lon_1D,lon_1D.shape)

#set domain for plotting in degrees. Western hemisphere is >180 because the ERA5 files are in that format. Generally sould match format of your data
if plt_bounds=='Auto':
    minlat=np.nanmin(lat_1D)
    maxlat=np.nanmax(lat_1D)
    minlon=np.nanmin(lon_1D)
    maxlon=np.nanmax(lon_1D)
elif plt_bounds=='Manual':
    minlat=32#np.nanmin(lat_1D)
    maxlat=50#p.nanmax(lat_1D)
    minlon=254#np.nanmin(lon_1D)
    maxlon=266#np.nanmax(lon_1D)
print(minlat,maxlat)
print(minlon,maxlon)

#exit()

print(lat_1D.shape,lon_1D.shape)
print(Cats.shape)
Y,T,I,J=Cats.shape
print(Y,T,I,J)
#for i in range(20,T,24):
#    print(i)
#    print(Cats[::,i,ref_lat,ref_lon])
#exit()

#Calculate data that are generally hourly (or dependent on hourly data)
#
#Data calculated in this section:
#	Trend in WBGT for each hour of the year
#	Smoothed trend in WBGT for each hour of the year
#	Trend of smoothed WBGT for each hour of the year
#
#Plots made in this section:
#	Smoothed trend in WBGT for specified hour of the day for each day of the year for a single point
#	Trend in Smoothed WBGT for specified hour of the day for each day of the year for a single point (and 3 points)
#	Map of smoothed trend of WBGT for a specified date and time
#	Map of trend of smoothed WBGT for a specified date and time

if HourlyClimo:
    print('###############')
    print('# HourlyClimo #')
    print('###############')

    #Set variables only needed for this psuedo-function
    time=plot_time #Set hour of the day we want to plot
    period=31 #Window of running mean in days
    plotHourlyClimo=False
    #calc trend by hour
    Cats.shape=(Y,T,I*J)
    Cats.shape=(Y,T*I*J)
    print(Cats.shape)

    #calculate trend of unsmoothed WBGT
    E,xhat=LinRegress([np.arange(Start_year,End_year+1)],Cats)
    print(xhat.shape,xhat.dtype)
    Cats.shape=(Y,T,I*J)
    xhat.shape=(2,T,I*J)
    Cats.shape=(Y,T,I,J)
    xhat.shape=(2,T,I,J)

    print(np.nanmax(xhat[0]))
    print(np.where(xhat[0]==np.nanmax(xhat[0])))
#    ind=np.where(xhat[0,:,ref_lat,ref_lon]==np.nanmax(xhat[0,:,ref_lat,ref_lon]))[0][0]
    #ind=ind
    #print(ind)
    data=xhat[0,time::24,::,::]*trend_scale#,30,65]*25
    data_sm=np.ones(data.shape,dtype=np.float32)*np.nan #smoothed after trends
    Cats_sm=np.ones((Y,num_days,I,J),dtype=np.float32)*np.nan #smoothed before trend
    Cats_time=Cats[::,time::24,::,::]

    print(Cats.shape)
    print(minlats[0],maxlats[0],minlons[0],maxlons[0])

    Cats_S=Cats[::,time::24,minlats[0]:maxlats[0]:,minlons[0]:maxlons[0]]
    Cats_C=Cats[::,time::24,minlats[1]:maxlats[1]:,minlons[1]:maxlons[1]]
    Cats_N=Cats[::,time::24,minlats[2]:maxlats[2]:,minlons[2]:maxlons[2]]
    print(Cats_S)

    Cats_S=np.nanmean(Cats_S,axis=(2,3))
    Cats_C=np.nanmean(Cats_C,axis=(2,3))
    Cats_N=np.nanmean(Cats_N,axis=(2,3))
    print(Cats_S)
    Cats_S_sm=np.ones((Y,num_days))*np.nan
    Cats_C_sm=np.ones((Y,num_days))*np.nan
    Cats_N_sm=np.ones((Y,num_days))*np.nan

    #smoothe trends and WBGT
    if period==31:
        print(period)
        for i in range(15,352):
            print(i)
            data_sm[i]=np.nanmean(data[i-15:i+16],axis=0)
            Cats_sm[::,i,::,::]=np.nanmean(Cats_time[::,i-15:i+16:,::,::],axis=1) 
            #print(Cats_sm.shape,np.nanmean(Cats_time[::,i-15:i+16:,::,::],axis=1).shape)
            #print(Cats_sm[::,i,ref_lat,ref_lon])
            #print(Cats_time[::,i,ref_lat,ref_lon])
            Cats_S_sm[::,i]=np.nanmean(Cats_S[::,i-15:i+16:],axis=1)
            Cats_C_sm[::,i]=np.nanmean(Cats_C[::,i-15:i+16:],axis=1)
            Cats_N_sm[::,i]=np.nanmean(Cats_N[::,i-15:i+16:],axis=1)
            print(Cats_S_sm)
            
    elif period==15:
        print(period)
        for i in range(7,359):
            data_sm[i,::,::]=np.nanmean(data[i-7:i+8,::,::],axis=0)
            Cats_sm[::,i,::,::]=np.nanmean(Cats_time[::,i-7:i+8:,::,::],axis=1)
    elif period==1:
        print(period)
        data_sm=data
        Cats_sm=Cats_time

    #Calculate trend of smoothed WBGT
    Cats_sm.shape=(Y,num_days,I*J)
    Cats_sm=Cats_sm.reshape(Y,num_days*I*J)
    E2,xhat2=LinRegress([np.arange(Start_year,End_year+1)],Cats_sm)
    ES,xhatS=LinRegress([np.arange(Start_year,End_year+1)],Cats_S_sm)
    EC,xhatC=LinRegress([np.arange(Start_year,End_year+1)],Cats_C_sm)
    EN,xhatN=LinRegress([np.arange(Start_year,End_year+1)],Cats_N_sm)
    Cats_sm.shape=(Y,num_days,I*J)
    xhat2.shape=(2,num_days,I*J)
    Cats_sm.shape=(Y,num_days,I,J)
    xhat2.shape=(2,num_days,I,J)
    Cat_trend=xhat2[0,::,::,::]*trend_scale
    Cat_trend_S=xhatS[0,::]*trend_scale
    Cat_trend_C=xhatC[0,::]*trend_scale
    Cat_trend_N=xhatN[0,::]*trend_scale
    print(xhatS)
    print(Cat_trend_S)

    del Cats_sm
    del E2
    del xhat2

    #Find where the max trend is
    ind=np.where(data_sm[120:240,ref_lat,ref_lon]==np.nanmax(data_sm[120:240,ref_lat,ref_lon]))[0][0]
    ind=200#224
    print("Max index for {} region: {}".format(region,str(ind)))
    #print(data_sm)

    #Plot trend for specified daily hour. Smoothed after trend calculated
    plt.figure()
    plt.plot(data_sm[::,ref_lat,ref_lon])
    plt.ylabel(r'$\Delta${}/{} years'.format(Plot_var,str(trend_scale)))
    plt.xticks(month_start_days,month_labels)
#    plt.title('{} day average trend in WBGT category'.format(str(period)))
    if period==1:
        if WBGT_type=='WBGT':
            plt.title('Daily trend in WBGT for day {} at {} UTC'.format(str(ind),str(time)))
        else:
            plt.title('Daily trend in WBGT category for day {} at {} UTC'.format(str(ind),str(time)))
    else:
        if WBGT_type=='WBGT':
            plt.title('{} day mean trend in WBGT for day {} at {} UTC'.format(str(period),str(ind),str(time)))            
        else:   
            plt.title('{} day mean trend in WBGT category for day {} at {} UTC'.format(str(period),str(ind),str(time)))
    plt.savefig('{}/cat_test_{}_{}UTC_{}-{}.png'.format(Save_dir,str(period),str(time),str(ref_lat),str(ref_lon)))
    plt.show()
    plt.close()

    #plot trend of smoothed WBGT
    plt.figure()
    plt.plot(Cat_trend[::,ref_lat,ref_lon])
    plt.ylabel(r'$\Delta${}/{} years'.format(Plot_var,str(trend_scale)))
    plt.xticks(month_start_days,month_labels)
#    plt.title('{} day average trend in WBGT category'.format(str(period)))
    if period==1:
        if WBGT_type=='WBGT':
            plt.title('Daily trend in WBGT for day {} at {} UTC'.format(str(ind),str(time)))
        else:
            plt.title('Daily trend in WBGT category for day {} at {} UTC'.format(str(ind),str(time)))
    else:
        if WBGT_type=='WBGT':
            plt.title('{} day mean trend in mean WBGT for day {} at {} UTC'.format(str(period),str(ind),str(time)))
        else:
            plt.title('{} day mean trend in mean WBGT category for day {} at {} UTC'.format(str(period),str(ind),str(time)))
    plt.savefig('{}/cat_test_presmoothed_{}_{}UTC_{}-{}.png'.format(Save_dir,str(period),str(time),str(ref_lat),str(ref_lon)))
    plt.show()
    plt.close()

    #plot trend of smoothed WBGT for all regions
    plt.figure()
    plt.axhline(0,color='k')
    plt.plot(Cat_trend[::,30,65],color='blue',label='South')
    plt.plot(Cat_trend[::,20,60],color='red',label='Central')
    plt.plot(Cat_trend[::,8,60],color='green',label='North')
    if WBGT_type=='WBGT':
        plt.ylabel(r'$\Delta${}/{} years (F)'.format(Plot_var,str(trend_scale)),fontsize=18)
    else:
        plt.ylim([-0.04,1.0])
        plt.yticks(np.arange(0,1.01,0.2))
        plt.ylabel(r'$\Delta${}/{} years'.format(Plot_var,str(trend_scale)),fontsize=18)

    plt.xticks(month_start_days,month_labels)
    plt.tick_params(labelsize=14)
#    plt.title('{} day average trend in WBGT category'.format(str(period)))
    #if period==1:
    #    if WBGT_type=='WBGT':
    #        plt.title('Daily trend in WBGT for day {} at {} UTC'.format(str(ind),str(time)))
    #    else:
    #        plt.title('Daily trend in WBGT category for day {} at {} UTC'.format(str(ind),str(time)))
    #else:
    #    if WBGT_type=='WBGT':
    #        plt.title('{} day mean trend in mean WBGT for day {} at {} UTC'.format(str(period),str(ind),str(time)))
    #    else:
    #        plt.title('{} day mean trend in mean WBGT category for day {} at {} UTC'.format(str(period),str(ind),str(time)))
    plt.tight_layout()
    plt.savefig('{}/cat_test_presmoothed_{}_{}UTC_all_nolegend.png'.format(Save_dir,str(period),str(time)),dpi=figdpi)
    plt.legend(loc=0,fontsize=20)
    plt.savefig('{}/cat_test_presmoothed_{}_{}UTC_all.png'.format(Save_dir,str(period),str(time)),dpi=figdpi)
    plt.show()
    plt.close()

    #Plot trend of smoothed regionally averaged WBGT for all regions
    plt.figure()
    plt.axhline(0,color='k')
    plt.plot(Cat_trend_S,color='blue',label='South')
    plt.plot(Cat_trend_C,color='red',label='Central')
    plt.plot(Cat_trend_N,color='green',label='North')
    print(Cat_trend_S)
    if WBGT_type=='WBGT':
        plt.ylim([-0.665,1.125])
        plt.yticks(np.arange(-0.6,1.01,0.2))
        plt.ylabel(r'$\Delta${}/{} years (F)'.format(Plot_var,str(trend_scale)),fontsize=18)
    else:
        plt.ylim([-0.04,1.0])
        plt.yticks(np.arange(0,1.01,0.2))
        plt.ylabel(r'$\Delta${}/{} years'.format(Plot_var,str(trend_scale)),fontsize=18)

    plt.xticks(month_start_days,month_labels)
    plt.tick_params(labelsize=14)
#    plt.title('{} day average trend in WBGT category'.format(str(period)))
    #if period==1:
    #    if WBGT_type=='WBGT':
    #        plt.title('Daily trend in WBGT for day {} at {} UTC'.format(str(ind),str(time)))
    #    else:
    #        plt.title('Daily trend in WBGT category for day {} at {} UTC'.format(str(ind),str(time)))
    #else:
    #    if WBGT_type=='WBGT':
    #        plt.title('{} day mean trend in mean WBGT for day {} at {} UTC'.format(str(period),str(ind),str(time)))
    #    else:
    #        plt.title('{} day mean trend in mean WBGT category for day {} at {} UTC'.format(str(period),str(ind),str(time)))
    plt.tight_layout()
    plt.savefig('{}/cat_test_presmoothed_regavg_{}_{}UTC_all_nolegend.png'.format(Save_dir,str(period),str(time)),dpi=figdpi)
    plt.legend(loc=0,fontsize=20)
    plt.savefig('{}/cat_test_presmoothed_regavg_{}_{}UTC_all.png'.format(Save_dir,str(period),str(time)),dpi=figdpi)
    plt.show()
    plt.close()

    #exit()


    print(np.nanmax(np.abs(data_sm[::,ref_lat,ref_lon]-Cat_trend[::,ref_lat,ref_lon])))

    if plotHourlyClimo:
        #Plot map of trend on the selected day (such as max trend at a point) for the post trend smoothed variable
        lon_2D,lat_2D=np.meshgrid(lon_1D,lat_1D)

        plt.figure()
        m=Basemap(resolution='i', projection='merc', llcrnrlat=minlat, urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon, lat_ts=10)
        m.drawcoastlines(linewidth=0.50)
        m.drawcountries(linewidth=0.25)
        m.drawstates(linewidth=0.15)
        m.drawmeridians(np.arange(0,360,5),labels=[False,False,False,True])
        m.drawparallels(np.arange(-90,90,5),labels=[True,False,False,False])
        #Xi,Yi=m(Xi,Yi)
        x,y=m(lon_2D,lat_2D)
        xi,yi=m(lon_2D[ref_lat,ref_lon],lat_2D[ref_lat,ref_lon])
        cmap=plt.get_cmap('bwr')
#        m.contourf(x,y,xhat[0][ind,:,:]*25,cmap=cmap,extend='both',levels=np.arange(-3.5,3.51,0.25))
#        m.contourf(x,y,data_sm[ind,:,:],cmap=cmap,extend='both',levels=np.arange(-3.5,3.51,0.25))
        if WBGT_type=='WBGT':
            cflevels=np.arange(-10,10.1,0.5)
            clevels=np.arange(-10,10.1,1)
        else:
            cflevels=np.arange(-6,6.1,1)
            clevels=np.arange(-8,8.1,1)
        m.contourf(x,y,data_sm[ind,:,:],cmap=cmap,extend='both',levels=cflevels)
        plt.colorbar()
        plt.scatter(xi,yi)
#        CS=m.contour(x,y,data_sm[ind,:,:],levels=np.arange(-4,4.1,1))
        CS=m.contour(x,y,data_sm[ind,:,:],colors='k',levels=clevels)
        plt.clabel(CS)
        plt.tick_params(labelsize=16)
        #Set Titles
        #if period==1:
        #    if WBGT_type=='WBGT':
        #        plt.title('Trend in WBGT for day {} at {} UTC'.format(str(ind),str(time)))
        #    else:
        #        plt.title('Trend in WBGT category for day {} at {} UTC'.format(str(ind),str(time)))
        #else:
        #    if WBGT_type=='WBGT':
        #        plt.title('{} day mean trend in WBGT for day {} at {} UTC'.format(str(period),str(ind),str(time)))
        #    else:
        #        plt.title('{} day mean trend in WBGT category for day {} at {} UTC'.format(str(period),str(ind),str(time)))

        plt.savefig('{}/cat_test_hourly_map_{}_day_{}_{}_{}'.format(Save_dir,str(period),str(ind),str(ref_lat),str(ref_lon)),dpi=figdpi)
        plt.show()
        plt.close()
        print(ind)

        #Plot map of trend of smoothed WBGT
        plt.figure()
        m=Basemap(resolution='i', projection='merc', llcrnrlat=minlat, urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon, lat_ts=10)
        m.drawcoastlines(linewidth=0.50)
        m.drawcountries(linewidth=0.25)
        m.drawstates(linewidth=0.15)
        m.drawmeridians(np.arange(0,360,5),labels=[False,False,False,True])
        m.drawparallels(np.arange(-90,90,5),labels=[True,False,False,False])
        #Xi,Yi=m(Xi,Yi)
        x,y=m(lon_2D,lat_2D)
        xi,yi=m(lon_2D[ref_lat,ref_lon],lat_2D[ref_lat,ref_lon])
        cmap=plt.get_cmap('bwr')
#        m.contourf(x,y,xhat[0][ind,:,:]*25,cmap=cmap,extend='both',levels=np.arange(-3.5,3.51,0.25))
#        m.contourf(x,y,data_sm[ind,:,:],cmap=cmap,extend='both',levels=np.arange(-3.5,3.51,0.25))
        m.contourf(x,y,Cat_trend[ind,:,:],cmap=cmap,extend='both',levels=cflevels)
        plt.colorbar()
        plt.scatter(xi,yi)
#        CS=m.contour(x,y,data_sm[ind,:,:],levels=np.arange(-4,4.1,1))
        CS=m.contour(x,y,Cat_trend[ind,:,:],colors='k',levels=clevels)
        plt.clabel(CS)
        plt.tick_params(labelsize=16)
        #Set Titles
        #if period==1:
        #    if WBGT_type=='WBGT':
        #        plt.title('Trend in WBGT for day {} at {} UTC'.format(str(ind),str(time),str(ref_lat),str(ref_lon)))
        #    else:
        #        plt.title('Trend in WBGT category for day {} at {} UTC'.format(str(ind),str(time),str(ref_lat),str(ref_lon)))
        #else:
        #    if WBGT_type=='WBGT':
        #        plt.title('{} day trend in mean WBGT for day {} at {} UTC'.format(str(period),str(ind),str(time)))
        #    else:
        #        plt.title('{} day trend in mean WBGT category for day {} at {} UTC'.format(str(period),str(ind),str(time)))

        plt.savefig('{}/cat_test_hourly_map_{}_mean_day_{}_{}_{}'.format(Save_dir,str(period),str(ind),str(ref_lat),str(ref_lon)),dpi=figdpi)
        plt.show()
        plt.close()
        print(ind)
    #exit()
    del data_sm
#    del Cats_sm
#    del E2
#    del xhat2
    del Cat_trend

#sum data by hour

if (WBGT_type!='WBGT') & (WBGT_type!='input'):
    #set array for each category which will be set to 1 if that hour was in that category and 0 if not. This will make statistics easier later
    HourlyCount=np.zeros((5,T,I,J),dtype=np.float32)
    for y in range(num_years):
        print(y)
        for c in range(num_cats):
            print(c)
            indc=np.where(Cats[y]==c)
            HourlyCount[c][indc]=HourlyCount[c][indc]+1
            del indc
#    print(1)
#    ind1=np.where(Cats[y]==1)
#    HourlyCount[1][ind0]=HourlyCount[1][ind0]+1
#    print(2)
#    ind2=np.where(Cats[y]==2)
#    HourlyCount[2][ind0]=HourlyCount[2][ind0]+1
#    print(3)
#    ind3=np.where(Cats[y]==3)
#    HourlyCount[3][ind0]=HourlyCount[3][ind0]+1
#    print(4)
#    ind4=np.where(Cats[4]==0)
#    HourlyCount[4][ind0]=HourlyCount[4][ind0]+1   

#sum data by hour within each month in each year

#calc trend by hour within month

#sum data by hour within month

#Calculate data that are generally summariezed by month (or based on monthly data)
#
#Data calculated in this section:
#       a
#       b
#       c
#
#Plots made in this section:
#       a
#       b
#       c

if MonthlyClimo:
    print('################')
    print('# MonthlyClimo #')
    print('################')

    time=21
    if WBGT_type=='WBGT':
        MonthlyMeanCat=np.ones((num_months*num_hours,I,J))*np.nan
        AnnMonthlyMeanCat=np.ones((num_years,num_months*num_hours,I,J))*np.nan
        MonthStart=0
        for m in range(num_months):
            MonthEnd=expected_ei[m]
            for h in range(num_hours):
                indc=np.nanmean((Cats[::,MonthStart+h:MonthEnd:24,::,::]),axis=(0,1))
                #print("indc shape",indc.shape)
                MonthlyMeanCat[m*num_hours+h,::,::]=indc
                #MonthlyCount[c,m*num_hours+h,::,::]=np.nansum(HourlyCount[c,MonthStart+h:MonthEnd:24,::,::],axis=0)
                del indc
                for y in range(num_years):
                    indm=np.nanmean((Cats[y,MonthStart+h:MonthEnd:24,::,::]),axis=0)
                    AnnMonthlyMeanCat[y,m*num_hours+h,::,::]=indm
                    del indm
            MonthStart=MonthEnd*1
        AnnMonthlyMeanCat.shape=(num_years,num_months*num_hours,I*J)
        AnnMonthlyMeanCat.shape=(num_years,num_months*num_hours*I*J)
        E,xhat=LinRegress([np.arange(Start_year,End_year+1)],AnnMonthlyMeanCat)
        MonthlyTrend=xhat[0]        

#        AnnMonthlyMeanCat.shape=(num_years,num_months*num_hours,I*J)
        MonthlyTrend.shape=(num_months*num_hours,I*J)
#        AnnMonthlyMeanCat.shape=(num_years,num_months*num_hours,I,J)
        MonthlyTrend.shape=(num_months*num_hours,I,J)

        del AnnMonthlyMeanCat
        del E
        del xhat

        plt.figure()
        plt.plot(MonthlyTrend[time::24,ref_lat,ref_lon]*trend_scale)
        plt.xticks(np.arange(12),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
        plt.title('{} UTC trend of each category by month'.format(str(time)))
        plt.savefig('{}/cat_test_monthly_trend_{}z_{}_{}.png'.format(Save_dir,str(time),str(ref_lat),str(ref_lon)))
        plt.show()
        plt.close()

    else:
        AnnMonthlyCount=np.zeros((num_cats,num_years,num_months*num_hours,I,J))
        for y in range(num_years):
            MonthStart=0
            for m in range(num_months):
                MonthEnd=expected_ei[m]
                for h in range(num_hours):
                    for c in range(num_cats):
                        indc=np.nansum((Cats[y,MonthStart+h:MonthEnd:24,::,::]==c),axis=0)
                        #print("indc shape",indc.shape)
                        AnnMonthlyCount[c,y,m*num_hours+h,::,::]=indc
                        #MonthlyCount[c,m*num_hours+h,::,::]=np.nansum(HourlyCount[c,MonthStart+h:MonthEnd:24,::,::],axis=0)
                        del indc
                MonthStart=MonthEnd*1

        #print(AnnMonthlyCount)
        MonthlyCount=np.nansum(AnnMonthlyCount,axis=1)
        print(MonthlyCount.shape)
        print(MonthlyCount[0,::,30,65])

        AnnMonthlyCount.shape=(num_cats,num_years,num_months*num_hours,I*J)
        AnnMonthlyCount.shape=(num_cats,num_years,num_months*num_hours*I*J)
        MonthlyTrend=np.ones((num_cats,num_months*num_hours*I*J))*np.nan
        for c in range(num_cats):
            E,xhat=LinRegress([np.arange(Start_year,End_year+1)],AnnMonthlyCount[c,::,::])
            MonthlyTrend[c,::]=xhat[0]
            del E
            del xhat
#        AnnMonthlyCount.shape=(num_cats,num_years,num_months*num_hours,I*J)
        MonthlyTrend.shape=(num_cats,num_months*num_hours,I*J) 
#        AnnMonthlyCount.shape=(num_cats,num_years,num_months*num_hours,I,J)
        MonthlyTrend.shape=(num_cats,num_months*num_hours,I,J) 

        del AnnMonthlyCount

        MonthlyCount=MonthlyCount/(31*num_years)
        MonthlyCount[::,24:48,::,::]=MonthlyCount[::,24:48,::,::]*(31*num_years)/((28*num_years)+num_leap_years)
        MonthlyCount[::,72:96,::,::]=MonthlyCount[::,72:96,::,::]*(31*num_years)/(30*num_years)
        MonthlyCount[::,120:144,::,::]=MonthlyCount[::,120:144,::,::]*(31*num_years)/(30*num_years)
        MonthlyCount[::,192:216,::,::]=MonthlyCount[::,192:216,::,::]*(31*num_years)/(30*num_years)
        MonthlyCount[::,240:264,::,::]=MonthlyCount[::,240:264,::,::]*(31*num_years)/(30*num_years)

        plt.figure()
        plt.plot(MonthlyCount[0,time::24,ref_lat,ref_lon],label=0)
        plt.plot(MonthlyCount[1,time::24,ref_lat,ref_lon],label=1)
        plt.plot(MonthlyCount[2,time::24,ref_lat,ref_lon],label=2)
        plt.plot(MonthlyCount[3,time::24,ref_lat,ref_lon],label=3)
        plt.plot(MonthlyCount[4,time::24,ref_lat,ref_lon],label=4)
        plt.legend()
        plt.xticks(np.arange(12),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
        plt.title('Mean {} UTC frequency of each category by month'.format(str(time)))
        plt.savefig('{}/cat_test_monthly_{}z_{}_{}.png'.format(Save_dir,str(time),str(ref_lat),str(ref_lon)))
        plt.show()
        plt.close()

        plt.figure()
        plt.plot(MonthlyCount[0,::,ref_lat,ref_lon],label=0)
        plt.plot(MonthlyCount[1,::,ref_lat,ref_lon],label=1)
        plt.plot(MonthlyCount[2,::,ref_lat,ref_lon],label=2)
        plt.plot(MonthlyCount[3,::,ref_lat,ref_lon],label=3)
        plt.plot(MonthlyCount[4,::,ref_lat,ref_lon],label=4)
        plt.legend()
        plt.xticks(np.arange(0,288,24),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
        plt.title('Mean frequency of each category by hour by month'.format(str(time)))
        plt.savefig('{}/cat_test_monthly_all_{}_{}.png'.format(Save_dir,str(ref_lat),str(ref_lon)))
        plt.show()
        plt.close()

        plt.figure()
        plt.plot(MonthlyTrend[0,time::24,ref_lat,ref_lon]*25,label=0)
        plt.plot(MonthlyTrend[1,time::24,ref_lat,ref_lon]*25,label=1)
        plt.plot(MonthlyTrend[2,time::24,ref_lat,ref_lon]*25,label=2)
        plt.plot(MonthlyTrend[3,time::24,ref_lat,ref_lon]*25,label=3)
        plt.plot(MonthlyTrend[4,time::24,ref_lat,ref_lon]*25,label=4)
        plt.legend()
        plt.xticks(np.arange(12),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
        plt.title('{} UTC trend of each category by month'.format(str(time)))
        plt.savefig('{}/cat_test_monthly_trend_{}z_{}_{}.png'.format(Save_dir,str(time),str(ref_lat),str(ref_lon)))
        plt.show()
        plt.close()

        MonthlyMeanCat=(MonthlyCount[4,::,::,::]*4)+(MonthlyCount[3,::,::,::]*3)+(MonthlyCount[2,::,::,::]*2)+(MonthlyCount[1,::,::,::]*1)

        del MonthlyCount
    del MonthlyTrend

    plt.figure()
    plt.plot(MonthlyMeanCat[::,ref_lat,ref_lon])
    plt.xticks(np.arange(0,288,24),['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
    if WBGT_type=='WBGT':
       plt.title('Mean WBGT by hour by month')
    else:
        plt.title('Mean WBGT category by hour by month')
    plt.savefig('{}/cat_test_monthly_mean_{}_{}.png'.format(Save_dir,str(time),str(ref_lat),str(ref_lon)))
    plt.show()
    plt.close()

    for i in range(4,9):
       print('This better be plotting. Why didn\'t it')
       plt.figure()
#       plt.plot(MonthlyMeanCat[i*24:(i+1)*24,ref_lat,ref_lon])
       plt.plot(MonthlyMeanCat[i*24:(i+1)*24,30,65])
       plt.xticks(np.arange(0,24,6),['0Z','6Z','12Z','18Z'],fontsize=16)
       if WBGT_type=='WBGT':
           plt.yticks(np.arange(56,84.1,4),fontsize=16)
       elif WBGT_type=='Cat':
           plt.yticks(np.arange(0,1.81,0.2),fontsize=16)
           plt.ylim([-0.04,1.84])
       elif WBGT_type=='ModCat': 
           plt.yticks(np.arange(0,1.81,0.2),fontsize=16)
           plt.ylim([-0.04,1.84])     
       plt.xlabel('Time (UTC)',fontsize=20)
       if WBGT_type=='WBGT':
           plt.ylabel('WBGT (F)',fontsize=20)
           plt.title('{} mean WBGT by hour'.format(month_labels[i]),fontsize=24)
       else:
           plt.ylabel('Mean WBGT Category',fontsize=20)
           plt.title('{} mean WBGT category by hour'.format(month_labels[i]),fontsize=24)
       plt.tight_layout()
       plt.savefig('{}/cat_test_monthly_mean_{}_{}_{}.png'.format(Save_dir,month_labels[i],str(30),str(65))) 
       plt.show()
       plt.close()

       plt.figure()
#       plt.plot(MonthlyMeanCat[i*24:(i+1)*24,ref_lat,ref_lon])
       plt.plot(MonthlyMeanCat[i*24:(i+1)*24,20,60])
       plt.xticks(np.arange(0,24,6),['0Z','6Z','12Z','18Z'],fontsize=16)
       if WBGT_type=='WBGT':
           plt.yticks(np.arange(60,84.1,4),fontsize=16)
       elif WBGT_type=='Cat':
           plt.yticks(np.arange(0,1.81,0.2),fontsize=16)
           plt.ylim([-0.04,1.84])
       elif WBGT_type=='ModCat':
           plt.yticks(np.arange(0,1.81,0.2),fontsize=16)
           plt.ylim([-0.04,1.84])
       plt.xlabel('Time (UTC)',fontsize=20)
       if WBGT_type=='WBGT':
           plt.ylabel('WBGT (F)',fontsize=20)
           plt.title('{} mean WBGT by hour'.format(month_labels[i]),fontsize=24)
       else:
           plt.ylabel('Mean WBGT Category',fontsize=20)
           plt.title('{} mean WBGT category by hour'.format(month_labels[i]),fontsize=24)
       plt.tight_layout()
       plt.savefig('{}/cat_test_monthly_mean_{}_{}_{}.png'.format(Save_dir,month_labels[i],str(20),str(65)))
       plt.show()
       plt.close()

       plt.figure()
#       plt.plot(MonthlyMeanCat[i*24:(i+1)*24,ref_lat,ref_lon])
       plt.plot(MonthlyMeanCat[i*24:(i+1)*24,8,60])
       plt.xticks(np.arange(0,24,6),['0Z','6Z','12Z','18Z'],fontsize=16)
       if WBGT_type=='WBGT':
           plt.yticks(np.arange(60,84.1,4),fontsize=16)
       elif WBGT_type=='Cat':
           plt.yticks(np.arange(0,1.81,0.2),fontsize=16)
           plt.ylim([-0.04,1.84])
       elif WBGT_type=='ModCat':
           plt.yticks(np.arange(0,1.81,0.2),fontsize=16)
           plt.ylim([-0.04,1.84])
       plt.xlabel('Time (UTC)',fontsize=20)
       if WBGT_type=='WBGT':
           plt.ylabel('WBGT (F)',fontsize=20)
           plt.title('{} mean WBGT by hour'.format(month_labels[i]),fontsize=24)
       else:
           plt.ylabel('Mean WBGT Category',fontsize=20)
           plt.title('{} mean WBGT category by hour'.format(month_labels[i]),fontsize=24)
       plt.tight_layout()
       plt.savefig('{}/cat_test_monthly_mean_{}_{}_{}.png'.format(Save_dir,month_labels[i],str(8),str(60)))
       plt.show()
       plt.close()

       plt.figure()
       plt.plot(MonthlyMeanCat[i*24:(i+1)*24,30,65],color='blue',label='South')
       plt.plot(MonthlyMeanCat[i*24:(i+1)*24,20,60],color='red',label='Central')
       plt.plot(MonthlyMeanCat[i*24:(i+1)*24,8,60],color='green',label='North')
       plt.xticks(np.arange(0,24,6),['0','6','12','18'],fontsize=16)
       if WBGT_type=='WBGT':
           plt.yticks(np.arange(56,88.1,4),fontsize=16)
           plt.ylim([55,89])
       elif WBGT_type=='Cat':
           plt.yticks(np.arange(0,2.01,0.2),fontsize=16)
           plt.ylim([-0.04,2.04])
       elif WBGT_type=='ModCat':
           plt.yticks(np.arange(0,2.01,0.2),fontsize=16)
           plt.ylim([-0.04,2.04])
       plt.xlabel('Time (UTC)',fontsize=20)
       if WBGT_type=='WBGT':
           plt.ylabel('WBGT (F)',fontsize=20)
           #plt.title('{} mean WBGT by hour'.format(month_labels[i]),fontsize=24)
       else:
           plt.ylabel('Mean WBGT Category',fontsize=20)
           #plt.title('{} mean WBGT category by hour'.format(month_labels[i]),fontsize=24)
       plt.legend(loc=2,fontsize=20)
       plt.tight_layout()
       plt.savefig('{}/cat_test_monthly_mean_{}_all.png'.format(Save_dir,month_labels[i]),dpi=figdpi)
       plt.show()
       plt.close()

       plt.figure()
       plt.plot(MonthlyMeanCat[i*24:(i+1)*24,30,65],color='blue',label='South')
       plt.plot(MonthlyMeanCat[i*24:(i+1)*24,20,60],color='red',label='Central')
       plt.plot(MonthlyMeanCat[i*24:(i+1)*24,8,60],color='green',label='North')
       plt.xticks(np.arange(0,24,6),['0','6','12','18'],fontsize=16)
       if WBGT_type=='WBGT':
           plt.yticks(np.arange(56,88.1,4),fontsize=16)
           plt.ylim([55,89])
       elif WBGT_type=='Cat':
           plt.yticks(np.arange(0,2.01,0.2),fontsize=16)
           plt.ylim([-0.04,2.04])
       elif WBGT_type=='ModCat':
           plt.yticks(np.arange(0,2.01,0.2),fontsize=16)
           plt.ylim([-0.04,2.04])
       plt.xlabel('Time (UTC)',fontsize=20)
       if WBGT_type=='WBGT':
           plt.ylabel('WBGT (F)',fontsize=20)
           #plt.title('{} mean WBGT by hour'.format(month_labels[i]),fontsize=24)
       else:
           plt.ylabel('Mean WBGT Category',fontsize=20)
           #plt.title('{} mean WBGT category by hour'.format(month_labels[i]),fontsize=24)
#       plt.legend(loc=2)
       plt.tight_layout()
       plt.savefig('{}/cat_test_monthly_mean_{}_all_nolegend.png'.format(Save_dir,month_labels[i]),dpi=figdpi)
       plt.show()
       plt.close()


       plt.figure()
       plt.plot(np.nanmean(MonthlyMeanCat[i*24:(i+1)*24,minlats[0]:maxlats[0]:,minlons[0]:maxlons[0]:],axis=(1,2)),color='blue',label='South')
       plt.plot(np.nanmean(MonthlyMeanCat[i*24:(i+1)*24,minlats[1]:maxlats[1]:,minlons[1]:maxlons[1]:],axis=(1,2)),color='red',label='Central')
       plt.plot(np.nanmean(MonthlyMeanCat[i*24:(i+1)*24,minlats[2]:maxlats[2]:,minlons[2]:maxlons[2]:],axis=(1,2)),color='green',label='North')
       plt.xticks(np.arange(0,24,6),['0','6','12','18'],fontsize=16)
       if WBGT_type=='WBGT':
           plt.yticks(np.arange(56,88.1,4),fontsize=16)
           plt.ylim([55,89])
       elif WBGT_type=='Cat':
           plt.yticks(np.arange(0,2.01,0.2),fontsize=16)
           plt.ylim([-0.04,2.04])
       elif WBGT_type=='ModCat':
           plt.yticks(np.arange(0,2.01,0.2),fontsize=16)
           plt.ylim([-0.04,2.04])
       plt.xlabel('Time (UTC)',fontsize=20)
       if WBGT_type=='WBGT':
           plt.ylabel('WBGT (F)',fontsize=20)
           #plt.title('{} mean WBGT by hour'.format(month_labels[i]),fontsize=24)
       else:
           plt.ylabel('Mean WBGT Category',fontsize=20)
           #plt.title('{} mean WBGT category by hour'.format(month_labels[i]),fontsize=24)
#       plt.legend(loc=2)
       plt.tight_layout()
       plt.savefig('{}/cat_test_monthly_mean_{}_all_Regional_Mean_nolegend.png'.format(Save_dir,month_labels[i]),dpi=figdpi)
       plt.legend(loc=2,fontsize=20)
       plt.savefig('{}/cat_test_monthly_mean_{}_all_Regional_Mean.png'.format(Save_dir,month_labels[i]),dpi=figdpi)
       plt.show()
       plt.close()

    del MonthlyMeanCat

#    del AnnMonthlyMeanCat
#    del MonthlyTrend
#    del indc
#    del indm
#    del E
#    del xhat
#    del AnnMonthlyCount
#    del MonthlyCount
#    del MonthlyMeanCat

if DailyClimo:
    print('##############')
    print('# DailyClimo #')
    print('##############')

    Cat_Trends=True

#sum data by day each year
    AnnDailyCount=np.zeros((num_cats,num_years,num_days,I,J),dtype=np.float32)
    AnnDailyExt=np.ones((2,num_years,num_days,I,J),dtype=np.float32)*np.nan
    AnnDailyInt=np.ones((num_years,num_days,I,J),dtype=np.float32)*np.nan
#calc trend (for each cat? time in each cat? max/min cat?) by day
    for y in range(num_years):
        for d in range(num_days):
            indmax=np.nanmax(Cats[y,d*num_hours:(d+1)*num_hours,::,::],axis=0)
            indmin=np.nanmin(Cats[y,d*num_hours:(d+1)*num_hours,::,::],axis=0)
            AnnDailyExt[0,y,d,::,::]=indmin.astype(np.float32)
            AnnDailyExt[1,y,d,::,::]=indmax.astype(np.float32)
            del indmax
            del indmin

            if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                AnnDailyInt[y,d,::,::]=np.nansum(Cats[y,d*num_hours:(d+1)*num_hours,::,::],axis=0)
            else:
                AnnDailyInt[y,d,::,::]=np.nanmean(Cats[y,d*num_hours:(d+1)*num_hours,::,::],axis=0)

            for c in range(num_cats):
                indc=np.nansum((Cats[y,d*num_hours:(d+1)*num_hours,::,::]==c),axis=0)
                AnnDailyCount[c,y,d,::,::]=indc
            del indc
#sum data by day
    DailyCount=np.nansum(AnnDailyCount,axis=1)
    DailyMinCount=np.zeros((num_cats,num_days,I,J),dtype=np.float32)
    DailyMaxCount=np.zeros((num_cats,num_days,I,J),dtype=np.float32)  
    print(DailyCount[0,::,ref_lat,ref_lon])

    print('Starting loop to count daily max and min cats')
    for d in range(num_days):
        for c in range(num_cats):
            indcmin=np.nansum((AnnDailyExt[0,::,d,::,::]==c),axis=0)
            indcmax=np.nansum((AnnDailyExt[1,::,d,::,::]==c),axis=0)
            DailyMinCount[c,d,::,::]=indcmin.astype(np.float32)
            DailyMaxCount[c,d,::,::]=indcmax.astype(np.float32)

    del indcmin
    del indcmax

    if Cat_Trends==True:
        #MonMinCount=np.zeros((num_cats,num_years,num_months,I,J),dtype=np.float32)
        MonMaxCount=np.zeros((num_cats,num_years,num_months,I,J),dtype=np.float32)
        for y in range(num_years):
            for m in range(num_months):
                for c in range(num_cats):
                    #indcmin=np.nansum((AnnDailyExt[0,y,month_start_days[m]:int(expected_ei[m]/24):,::,::]==c),axis=0)
                    indcmax=np.nansum((AnnDailyExt[1,y,month_start_days[m]:int(expected_ei[m]/24):,::,::]==c),axis=0)
                    #MonMinCount[c,y,m,::,::]=indcmin.astype(np.float32)
                    MonMaxCount[c,y,m,::,::]=indcmax.astype(np.float32)
        #AnnMinCount=np.nansum(MonMinCount,axis=2)
        
        MJJASMaxCount=np.nansum(MonMaxCount[::,::,4:9:,::],axis=2)
        AnnMaxCount=np.nansum(MonMaxCount,axis=2)

        #print('AnnMinCount shape 1',AnnMinCount.shape)
        if RegionalMean:
            IJ_reg=len(maxlats)
            AnnMaxCount_Region=np.zeros((num_cats,num_years,IJ_reg))
            for i in range(IJ_reg):
                AnnMaxCount_Region[::,::,i]=np.nanmean(MJJASMaxCount[::,::,minlats[i]:maxlats[i],minlons[i]:maxlons[i]],axis=(2,3))

        Full_MC=False
        if Full_MC:
            IJ=I*J
            #AnnMinCount.shape=(num_cats,num_years,IJ)
            AnnMaxCount.shape=(num_cats,num_years,IJ)
            MJJASMaxCount.shape=(num_cats,num_years,IJ)
            #MonMinCount.shape=(num_cats,num_years,num_months,IJ)
            MonMaxCount.shape=(num_cats,num_years,num_months,IJ)
            #print('AnnMinCount shape 2a',AnnMinCount.shape)

        else:
            ref_lat_list=[30,20,8]
            ref_lon_list=[65,60,60]
            IJ=len(ref_lat_list)
            #print('AnnMinCount shape 2b',AnnMinCount.shape)
            #AnnMinCount_trim=np.zeros((num_cats,num_years,IJ))
            AnnMaxCount_trim=np.zeros((num_cats,num_years,IJ))
            MJJASMaxCount_trim=np.zeros((num_cats,num_years,IJ))
            #MonMinCount_trim=np.zeros((num_cats,num_years,num_months,IJ))
            MonMaxCount_trim=np.zeros((num_cats,num_years,num_months,IJ))

            for i in range(len(ref_lat_list)):
                #print('AnnMinCount shape 3',AnnMinCount.shape)
                #AnnMinCount_trim[::,::,i]=AnnMinCount[::,::,ref_lat_list[i],ref_lon_list[i]]
                AnnMaxCount_trim[::,::,i]=AnnMaxCount[::,::,ref_lat_list[i],ref_lon_list[i]]
                MJJASMaxCount_trim[::,::,i]=MJJASMaxCount[::,::,ref_lat_list[i],ref_lon_list[i]]
                #MonMinCount_trim[::,::,::,i]=MonMinCount[::,::,::,ref_lat_list[i],ref_lon_list[i]]
                MonMaxCount_trim[::,::,::,i]=MonMaxCount[::,::,::,ref_lat_list[i],ref_lon_list[i]]


            #AnnMinCount=AnnMinCount_trim.copy()
            AnnMaxCount=AnnMaxCount_trim.copy()
            MJJASMaxCount=MJJASMaxCount_trim.copy()
            #MonMinCount=MonMinCount_trim.copy()
            MonMaxCount=MonMaxCount_trim.copy()

            #del AnnMinCount_trim
            del AnnMaxCount_trim
            del MJJASMaxCount_trim
            #del MonMinCount_trim
            del MonMaxCount_trim

        
        EAnn=np.zeros((num_cats,num_years,2))
        xhatAnn=np.zeros((num_cats,2,IJ))
        EMJJAS=np.zeros((num_cats,num_years,2))
        xhatMJJAS=np.zeros((num_cats,2,IJ))
        EMon=np.zeros((num_cats,num_months,num_years,2))
        xhatMon=np.zeros((num_cats,num_months,2,IJ))
        AnnMaxTrend=np.zeros((num_cats,num_years,IJ))
        MJJASMaxTrend=np.zeros((num_cats,num_years,IJ))
        MonMaxTrend=np.zeros((num_cats,num_years,num_months,IJ))
        PAnn=np.ones((num_cats,IJ))*np.nan
        PAnn_tf=np.zeros((num_cats,IJ))
        PMJJAS=np.ones((num_cats,IJ))*np.nan
        PMJJAS_tf=np.zeros((num_cats,IJ))
        PMon=np.ones((num_cats,num_months,IJ))*np.nan
        PMon_tf=np.zeros((num_cats,num_months,IJ))
        if RegionalMean:
            AnnMaxTrend_Region=np.zeros((num_cats,num_years,IJ_reg))
            PAnn_Reg=np.ones((num_cats,IJ_reg))*np.nan
            PAnn_tf_Reg=np.zeros((num_cats,IJ_reg))
            EAnn_Reg=np.zeros((num_cats,num_years,2))
            xhatAnn_Reg=np.zeros((num_cats,2,IJ_reg))

        for c in range(num_cats):
            print(c)
            E,xhat=LinRegress([np.arange(Start_year,End_year+1)],AnnMaxCount[c])
            print(E.shape)
            EAnn[c]=E
            xhatAnn[c]=xhat
            print('SigTest')
            P,P_tf=SigTrend(AnnMaxCount[c],xhatAnn[c][0],5000)
            PAnn[c]=P
            PAnn_tf[c]=P_tf
            AnnMaxTrend[c]=np.dot(E,xhat)

            E,xhat=LinRegress([np.arange(Start_year,End_year+1)],MJJASMaxCount[c])
            print(E.shape)
            EMJJAS[c]=E
            xhatMJJAS[c]=xhat
            print('SigTest')
            P,P_tf=SigTrend(MJJASMaxCount[c],xhatMJJAS[c][0],5000)
            PMJJAS[c]=P
            PMJJAS_tf[c]=P_tf
            MJJASMaxTrend[c]=np.dot(E,xhat)

            if RegionalMean:
                E,xhat=LinRegress([np.arange(Start_year,End_year+1)],AnnMaxCount_Region[c])
                print(E.shape)
                EAnn_Reg[c]=E
                xhatAnn_Reg[c]=xhat
                print('SigTest')
                P,P_tf=SigTrend(AnnMaxCount_Region[c],xhatAnn_Reg[c][0],5000)
                PAnn_Reg[c]=P
                PAnn_tf_Reg[c]=P_tf
                AnnMaxTrend_Region[c]=np.dot(E,xhat)
                
            for m in range(num_months):
                print(m)
                E,xhat=LinRegress([np.arange(Start_year,End_year+1)],MonMaxCount[c,::,m,::])
                EMon[c,m]=E
                xhatMon[c,m]=xhat
                print('SigTest')
                P,P_tf=SigTrend(MonMaxCount[c,::,m,::],xhatMon[c,m][0],5000)
                PMon[c,m]=P
                PMon_tf[c,m]=P_tf
                MonMaxTrend[c,::,m,::]=np.dot(E,xhat)
        
       
        if Full_MC:
            #AnnMinCount.shape=(num_cats,num_years,I,J)
            AnnMaxCount.shape=(num_cats,num_years,I,J)
            MJJASMaxCount.shape=(num_cats,num_years,I,J)
            AnnMaxTrend.shape=(num_cats,num_years,I,J)
            MJJASMaxTrend.shape=(num_cats,num_years,I,J)
            PAnn.shape=(num_cats,I,J)
            PAnn_tf.shape=(num_cats,I,J)
            PMJJAS.shape=(num_cats,I,J)
            PMJJAS_tf.shape=(num_cats,I,J)
        
            #MonMinCount.shape=(num_cats,num_years,num_months,I,J)
            MonMaxCount.shape=(num_cats,num_years,num_months,I,J)
            MonMaxTrend.shape=(num_cats,num_years,num_months,I,J)
            PMon.shape=(num_cats,num_months,I,J)
            PMon_tf.shape=(num_cats,num_months,I,J)

    #del indcmin
    del indcmax

    AnnDailyCount.shape=(num_cats,num_years,num_days,I*J)
    AnnDailyCount.shape=(num_cats,num_years,num_days*I*J)

    print('Calculating daily trend') #Note: I dont use the daily trend for some reason
#    DailyTrend=np.ones((num_cats,num_days*I*J))*np.nan

#    for c in range(num_cats):
#        E,xhat=LinRegress([np.arange(Start_year,End_year+1)],AnnDailyCount[c,::,::])
#        DailyTrend[c,::]=xhat[0]
    AnnDailyInt.shape=(num_years,num_days,I*J)
    AnnDailyInt.shape=(num_years,num_days*I*J)
    Eint,xhat_int=LinRegress([np.arange(Start_year,End_year+1)],AnnDailyInt)

    print('Reshaping daily trends')
#    AnnDailyCount.shape=(num_cats,num_years,num_days,I*J)
    AnnDailyInt.shape=(num_years,num_days,I*J)
#    DailyTrend.shape=(num_cats,num_days,I*J)
    xhat_int.shape=(2,num_days,I*J)
#    AnnDailyCount.shape=(num_cats,num_years,num_days,I,J)
    AnnDailyInt.shape=(num_years,num_days,I,J)
#    DailyTrend.shape=(num_cats,num_days,I,J)
    xhat_int.shape=(2,num_days,I,J)

    del AnnDailyCount
    del Eint

    print('Calculating daily thresholds')

    if Pre_calc:
        if WBGT_type=='Cat':
            f='{}/{}/{}CatClimo.nc'.format(ERA5_path,Var_name_path,Var_name_path)
            #Varn, Varnm, lat_pre, lon_pre = ReadNetCDF(f,['time','WBGT90MaxRM','WBGT90IntRM'])#,chunk_start,chunk_end)
            Varn=xr.load_dataset(f)
            Pre_vars=['WBGT90MaxRM','WBGT90IntRM']
        elif WBGT_type=='ModCat':
            f='{}/{}/{}ModCatClimo.nc'.format(ERA5_path,Var_name_path,Var_name_path,y)
            #Varn, Varnm, lat_pre, lon_pre = ReadNetCDF(f,['time','WBGT90MaxRM','WBGT90IntRM'])#,chunk_start,chunk_end)
            Varn=xr.load_dataset(f)
            Pre_vars=['WBGT90MaxRM','WBGT90IntRM']

        elif WBGT_type=='WBGT':
            f='{}/{}/{}Climo.nc'.format(ERA5_path,Var_name_path,Var_name_path)
            #Varn, Varnm, lat_pre, lon_pre = ReadNetCDF(f,['time','WBGT90RM','WBGT50RM'])#,chunk_start,chunk_end)
            Varn=xr.load_dataset(f)
            Pre_vars=['WBGT90RM','WBGT50RM']
        #time_pre=Varn[0]

        Varn_units=[]
        for Pre_var in Pre_vars:
            mv=Varn[Pre_var].encoding['missing_value']
            Varn_units.append(Varn[Pre_var].units)
            Varn[Pre_var]=xr.where(Varn[Pre_var]==mv,np.nan,Varn[Pre_var])

        if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
            Daily_max_thresh=np.rint(Varn[Pre_vars[0]].values).astype(np.float32)
            Daily_int_thresh=np.rint(Varn[Pre_vars[1]].values).astype(np.float32)
        elif WBGT_type=='WBGT':
            Daily_max_thresh=Varn[Pre_vars[0]].astype(np.float32)
            Daily_int_thresh=Varn[Pre_vars[1]].astype(np.float32)
            Daily_max_thresh=Convert_T_Unit(Daily_max_thresh,Varn_units[0],'F')
            Daily_int_thresh=Convert_T_Unit(Daily_int_thresh,Varn_units[1],'F')
        del Varn
        #del Varnm

    else:
        Daily_int_thresh=np.nanpercentile(AnnDailyInt,90,axis=0,interpolation='linear').astype(np.float32)
        Daily_max_thresh=np.nanpercentile(AnnDailyExt[1,::,::,::,::],90,axis=0,interpolation='linear').astype(np.float32)

    if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
        Daily_int_thresh[np.where(Daily_int_thresh<1)]=1
        Daily_max_thresh[np.where(Daily_max_thresh<1)]=1

    DailyCount=DailyCount/num_years
    DailyCount[::,59,::,::]=DailyCount[::,59,::,::]*num_years/num_leap_years

    PlotDailyClimo=True
    if PlotDailyClimo:
        if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
            plt.figure()
            Risk_labels=['No Risk','Low','Moderate','High','Extreme']
            Risk_colors=['tab:blue','tab:orange','tab:green','tab:red','tab:purple']
            if RegionalMean:
                print('Regional Mean')
                for p in range(IJ_reg):
                    plt.figure()
                    for c in range(num_cats):
                        print(PAnn_Reg[c,p],PAnn_tf_Reg[c,p])
                        if PAnn_tf_Reg[c,p]:
                            lw=4
                        else: lw=2
                        print(lw)
                        plt.plot(np.arange(Start_year,End_year+1),AnnMaxCount_Region[c,::,p],color=Risk_colors[c],label=Risk_labels[c])
                        plt.plot(np.arange(Start_year,End_year+1),AnnMaxTrend_Region[c,::,p],color=Risk_colors[c],linestyle='--',linewidth=lw)
                    plt.ylim([-10,160])
                    plt.xticks(fontsize=16)
                    plt.yticks(fontsize=16)
                    plt.savefig('{}/Ann_Max_Cat_Count_{}_Mean_no_legend.png'.format(Save_dir,Region_Labels[p]))
                    plt.legend()
                    plt.savefig('{}/Ann_Max_Cat_Count_{}_Mean.png'.format(Save_dir,Region_Labels[p]))
                    plt.show()
                    plt.close()

            if Full_MC:
                plt.figure()
                for c in range(num_cats):
                    plt.plot(np.arange(Start_year,End_year+1),AnnMaxCount[c,::,ref_lat,ref_lon],color=Risk_colors[c],label=Risk_labels[c])
                    plt.plot(np.arange(Start_year,End_year+1),AnnMaxTrend[c,::,ref_lat,ref_lon],color=Risk_colors[c],linestyle='--')
                plt.ylim([-10,365])
                plt.savefig('{}/Ann_Max_Cat_Count_{}_{}_no_legend.png'.format(Save_dir,str(ref_lat),str(ref_lon)))
                plt.legend()
                plt.savefig('{}/Ann_Max_Cat_Count_{}_{}.png'.format(Save_dir,str(ref_lat),str(ref_lon)))
                plt.show()
                plt.close()

                plt.figure()    
                for c in range(num_cats):
                    plt.plot(np.arange(Start_year,End_year+1),MJJASMaxCount[c,::,ref_lat,ref_lon],color=Risk_colors[c],label=Risk_labels[c])
                    plt.plot(np.arange(Start_year,End_year+1),MJJASMaxTrend[c,::,ref_lat,ref_lon],color=Risk_colors[c],linestyle='--')
                plt.ylim([-10,160]) #-10,153?
                plt.savefig('{}/MJJAS_Max_Cat_Count_{}_{}_no_legend.png'.format(Save_dir,str(ref_lat),str(ref_lon)))
                plt.legend()
                plt.savefig('{}/MJJAS_Max_Cat_Count_{}_{}.png'.format(Save_dir,str(ref_lat),str(ref_lon)))
                plt.show()
                plt.close()
    
                for m in range(num_months):
                    plt.figure()
                    for c in range(num_cats):
                        plt.plot(np.arange(Start_year,End_year+1),MonMaxCount[c,::,m,ref_lat,ref_lon],color=Risk_colors[c],label=Risk_labels[c])
                        plt.plot(np.arange(Start_year,End_year+1),MonMaxCount[c,::,m,ref_lat,ref_lon],color=Risk_colors[c],linestyle='--')
                    plt.ylim([0,35])
                    plt.xlabel(fontsize=16)
                    plt.ylabel(fontsize=16) 
                    plt.savefig('{}/{}_Max_Cat_Count_{}_{}_no_legend.png'.format(Save_dir,month_labels[m],str(ref_lat),str(ref_lon)))
                    plt.legend()
                    plt.savefig('{}/{}_Max_Cat_Count_{}_{}.png'.format(Save_dir,month_labels[m],str(ref_lat),str(ref_lon)))
                    plt.show()
                    plt.close()
            else:
                for p in range(len(ref_lat_list)):
                    
                    plt.figure()
                    for c in range(num_cats):
                        print(PAnn[c,p],PAnn_tf[c,p])
                        if PAnn_tf[c,p]:
                            lw=4
                        else: lw=2
                        print(lw)
                        plt.plot(np.arange(Start_year,End_year+1),AnnMaxCount[c,::,p],color=Risk_colors[c],label=Risk_labels[c])
                        plt.plot(np.arange(Start_year,End_year+1),AnnMaxTrend[c,::,p],color=Risk_colors[c],linestyle='--',linewidth=lw)
                    plt.ylim([-10,365])
                    plt.xticks(fontsize=16)
                    plt.yticks(fontsize=16)
                    plt.savefig('{}/Ann_Max_Cat_Count_{}_{}_no_legend.png'.format(Save_dir,str(ref_lat_list[p]),str(ref_lon_list[p])))
                    plt.legend()
                    plt.savefig('{}/Ann_Max_Cat_Count_{}_{}.png'.format(Save_dir,str(ref_lat_list[p]),str(ref_lon_list[p])))
                    plt.show()
                    plt.close()

                    plt.figure()
                    for c in range(num_cats):
                        print(PMJJAS[c,p],PMJJAS_tf[c,p])
                        if PMJJAS_tf[c,p]:
                            lw=4
                        else: lw=2
                        print(lw)
                        plt.plot(np.arange(Start_year,End_year+1),MJJASMaxCount[c,::,p],color=Risk_colors[c],label=Risk_labels[c])
                        plt.plot(np.arange(Start_year,End_year+1),MJJASMaxTrend[c,::,p],color=Risk_colors[c],linestyle='--',linewidth=lw)
                    plt.ylim([-10,160]) #-10,153?
                    plt.xticks(fontsize=16)
                    plt.yticks(fontsize=16)

                    plt.savefig('{}/MJJAS_Max_Cat_Count_{}_{}_no_legend.png'.format(Save_dir,str(ref_lat_list[p]),str(ref_lon_list[p])),dpi=figdpi)
                    plt.legend()
                    plt.savefig('{}/MJJAS_Max_Cat_Count_{}_{}.png'.format(Save_dir,str(ref_lat_list[p]),str(ref_lon_list[p])),dpi=figdpi)
                    plt.show()
                    plt.close()


                    """
                    for m in range(num_months):
                        plt.figure()
                        for c in range(num_cats):
                            print(PMon[c,m,p],PMon_tf[c,m,p])
                            if PMon_tf[c,m,p]:
                                lw=3
                            else: lw=2
                            print(lw)
                            plt.plot(np.arange(Start_year,End_year+1),MonMaxCount[c,::,m,p],color=Risk_colors[c],label=Risk_labels[c])
                            plt.plot(np.arange(Start_year,End_year+1),MonMaxTrend[c,::,m,p],color=Risk_colors[c],linestyle='--',linewidth=lw)
                        plt.ylim([-1,35])
                        plt.xticks(fontsize=16)
                        plt.yticks(fontsize=16)

                        plt.savefig('{}/{}_Max_Cat_Count_{}_{}_no_legend.png'.format(Save_dir,month_labels[m],str(ref_lat_list[p]),str(ref_lon_list[p])))
                        plt.legend()
                        plt.savefig('{}/{}_Max_Cat_Count_{}_{}.png'.format(Save_dir,month_labels[m],str(ref_lat_list[p]),str(ref_lon_list[p])))
                        plt.show()
                        plt.close()
                    """
    
            #del AnnMinCount
            del AnnMaxCount
            del MJJASMaxCount
            #del MonMinCount
            del MonMaxCount
            del AnnMaxTrend
            del MJJASMaxTrend
            del MonMaxTrend
            del PAnn
            del PAnn_tf
            del PMJJAS
            del PMJJAS_tf
            del PMon
            del PMon_tf
            #exit()           

            if RegionalMean:
                DailyMaxCount_Region=np.ones((num_cats,num_days,IJ_reg))*np.nan
                for r in range(IJ_reg):
                    DailyMaxCount_Region[::,::,r]=np.nanmean(DailyMaxCount[::,::,minlats[r]:maxlats[r],minlons[r]:maxlons[r]],axis=(2,3))

                DailyMaxCount_Region=DailyMaxCount_Region/num_years 
                DailyMaxCount_Region[::,59,::]=DailyMaxCount_Region[::,59,::]*num_years/num_leap_years
                for r in range(IJ_reg):
                    plt.figure()
                    plt.plot(DailyMaxCount_Region[0,::,r],label='No Risk')
                    plt.plot(DailyMaxCount_Region[1,::,r],label='Low')
                    plt.plot(DailyMaxCount_Region[2,::,r],label='Moderate')
                    plt.plot(DailyMaxCount_Region[3,::,r],label='High')
                    plt.plot(DailyMaxCount_Region[4,::,r],label='Extreme')
                    plt.legend()
                    plt.xticks(month_start_days,month_labels,fontsize=14)
                    plt.yticks(np.arange(0,1.01,0.2),fontsize=14)
                    #plt.title('Fraction of daily max\n category by day',fontsize=24)
                 #   plt.axhline(y=60)
                    plt.tight_layout()
                    plt.savefig('{}/cat_test_dailymax_mean_{}_Mean.png'.format(Save_dir,Region_Labels[r]),dpi=figdpi)
                    plt.show()
                    plt.close()


            plt.figure()
            plt.plot(DailyCount[0,::,ref_lat,ref_lon],label='0')
            plt.plot(DailyCount[1,::,ref_lat,ref_lon],label='1')
            plt.plot(DailyCount[2,::,ref_lat,ref_lon],label='2')
            plt.plot(DailyCount[3,::,ref_lat,ref_lon],label='3')
            plt.plot(DailyCount[4,::,ref_lat,ref_lon],label='4')
            plt.legend()
            plt.xticks(month_start_days,month_labels)
            plt.title('Mean hours in category by day')
    #        plt.axhline(y=60)
            plt.savefig('{}/cat_test_daily_mean_{}_{}.png'.format(Save_dir,str(ref_lat),str(ref_lon)))

            plt.show()
            plt.close()

            DailyMinCount=DailyMinCount/num_years
            DailyMinCount[::,59,::,::]=DailyMinCount[::,59,::,::]*num_years/num_leap_years
            plt.figure()
            plt.plot(DailyMinCount[0,::,ref_lat,ref_lon],label='0')
            plt.plot(DailyMinCount[1,::,ref_lat,ref_lon],label='1')
            plt.plot(DailyMinCount[2,::,ref_lat,ref_lon],label='2')
            plt.plot(DailyMinCount[3,::,ref_lat,ref_lon],label='3')
            plt.plot(DailyMinCount[4,::,ref_lat,ref_lon],label='4')
            plt.legend()
            plt.xticks(month_start_days,month_labels)
            plt.title('Fraction of daily min\n category by day')
            #    plt.axhline(y=60)
            plt.savefig('{}/cat_test_dailymin_mean_{}_{}.png'.format(Save_dir,str(ref_lat),str(ref_lon)))
            plt.show()
            plt.close()

            DailyMaxCount=DailyMaxCount/num_years
            DailyMaxCount[::,59,::,::]=DailyMaxCount[::,59,::,::]*num_years/num_leap_years
            plt.figure()
            plt.plot(DailyMaxCount[0,::,ref_lat,ref_lon],label='No Risk')
            plt.plot(DailyMaxCount[1,::,ref_lat,ref_lon],label='Low')
            plt.plot(DailyMaxCount[2,::,ref_lat,ref_lon],label='Moderate')
            plt.plot(DailyMaxCount[3,::,ref_lat,ref_lon],label='High')
            plt.plot(DailyMaxCount[4,::,ref_lat,ref_lon],label='Extreme')
            plt.legend()
            plt.xticks(month_start_days,month_labels,fontsize=14)
            plt.yticks(np.arange(0,1.01,0.2),fontsize=14)
            #plt.title('Fraction of daily max\n category by day',fontsize=24)
         #   plt.axhline(y=60)
            plt.tight_layout()
            plt.savefig('{}/cat_test_dailymax_mean_{}_{}.png'.format(Save_dir,str(ref_lat),str(ref_lon)),dpi=figdpi)
            plt.show()
            plt.close()

            plt.figure()
            plt.plot(xhat_int[0,::,ref_lat,ref_lon]*trend_scale)
            plt.xticks(month_start_days,month_labels)
            plt.tick_params(labelsize=16)
            plt.ylabel(r'$\Delta$Cat/{} years'.format(str(trend_scale)))
            plt.title('Trend in integrated category by day')
            plt.savefig('{}/cat_integral_daily_{}_{}.png'.format(Save_dir,str(ref_lat),str(ref_lon)))
            plt.show()
            plt.close()

    del DailyCount
    del DailyMinCount
    del DailyMaxCount
    del xhat_int

    print('Starting Monthly thresholds')
    MonthDailyMax_thresh=np.ones((num_days,I,J),dtype=np.float32)*np.nan
    MonthDailyInt_thresh=np.ones((num_days,I,J),dtype=np.float32)*np.nan

    if Pre_calc:
        if WBGT_type=='Cat':
            f='{}/{}/{}CatClimo.nc'.format(ERA5_path,Var_name_path,Var_name_path)
            #Varn, Varnm, lat_pre, lon_pre = ReadNetCDF(f,['time','WBGT90MaxM','WBGT90IntM'])#,chunk_start,chunk_end)
            Varn=xr.load_dataset(f)
            Pre_vars=['WBGT90MaxM','WBGT90IntM']

        elif WBGT_type=='ModCat':
            f='{}/{}/{}ModCatClimo.nc'.format(ERA5_path,Var_name_path,Var_name_path,y)
            #Varn, Varnm, lat_pre, lon_pre = ReadNetCDF(f,['time','WBGT90MaxM','WBGT90IntM'])#,chunk_start,chunk_end)
            Varn=xr.load_dataset(f)
            Pre_vars=['WBGT90MaxM','WBGT90IntM']

        elif WBGT_type=='WBGT':
            f='{}/{}/{}Climo.nc'.format(ERA5_path,Var_name_path,Var_name_path)
            #Varn, Varnm, lat_pre, lon_pre = ReadNetCDF(f,['time','WBGT90M','WBGT50M'])#,chunk_start,chunk_end)
            Varn=xr.load_dataset(f)
            Pre_vars=['WBGT90M','WBGT50M']

        Varn_units=[]
        for Pre_var in Pre_vars:
            mv=Varn[Pre_var].encoding['missing_value']
            Varn_units.append(Varn[Pre_var].units)
            Varn[Pre_var]=xr.where(Varn[Pre_var]==mv,np.nan,Varn[Pre_var])

        #time_pre=Varn[0]
        if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
            MonthDailyMax_thresh=np.rint(Varn[Pre_vars[0]]).astype(np.float32)
            MonthDailyInt_thresh=np.rint(Varn[Pre_vars[1]]).astype(np.float32)
        else:
            MonthDailyMax_thresh=Varn[Pre_vars[0]].astype(np.float32)
            MonthDailyInt_thresh=Varn[Pre_vars[1]].astype(np.float32)
            MonthDailyMax_thresh=Convert_T_Unit(MonthDailyMax_thresh,Varn_units[0],'F')
            MonthDailyInt_thresh=Convert_T_Unit(MonthDailyInt_thresh,Varn_units[1],'F')
        del Varn
        #del Varnm

    for i in range(6,num_months):
        print(month_labels[i])
        MonthAnnDailyExt=AnnDailyExt[1,::,month_start_days[i]:int(expected_ei[i]/24):,::,::]
        #if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
        MonthAnnDailyInt=AnnDailyInt[::,month_start_days[i]:int(expected_ei[i]/24):,::,::]

        print(MonthAnnDailyExt.shape)
        Y2,M,I2,J2=MonthAnnDailyExt.shape

        MonthAnnDailyExt=MonthAnnDailyExt.reshape(Y*M,I,J)
        #if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
        MonthAnnDailyInt=MonthAnnDailyInt.reshape(Y*M,I,J)
        if not Pre_calc:
            MonthDailyMax_thresh[month_start_days[i]:int(expected_ei[i]/24):,::,::]=np.nanpercentile(MonthAnnDailyExt,90,axis=0,interpolation='linear').astype(np.float32)
            MonthDailyInt_thresh[month_start_days[i]:int(expected_ei[i]/24):,::,::]=np.nanpercentile(MonthAnnDailyInt,90,axis=0,interpolation='linear').astype(np.float32)

#        MonthDailyMax_thresh[np.where(MonthDailyMax_thresh<1)]=1
#        MonthDailyInt_thresh[np.where(MonthDailyInt_thresh<1)]=1


        MonthAnnDailyExt.shape=(Y,M,I,J)
        #if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
        MonthAnnDailyInt.shape=(Y,M,I,J)
        
        print(MonthAnnDailyExt.shape)
        #MonthAnnDailyExt[np.where(MonthAnnDailyExt==0)]=np.nan
        #MonthAnnDailyCount.shape=((num_years*int(month_start_days[i+1]-month_start_days[i],)))
        #print(MonthAnnDailyCount.flatten())

        if i==6:#6:#6:#6:#6:#6:#0:
            lon_2D,lat_2D=np.meshgrid(lon_1D,lat_1D)

        if PlotDailyClimo:

            if RegionalMean:
                for r in range(IJ_reg):
                    plt.figure()
                    if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                        n,bins,patches=plt.hist(MonthAnnDailyExt[::,::,minlats[r]:maxlats[r]:,minlons[r]:maxlons[r]:].flatten(),bins=range(6),density=True)
                    else:
                        binsize=5
                        n,bins=np.histogram(MonthAnnDailyExt[::,::,minlats[r]:maxlats[r]:,minlons[r]:maxlons[r]:].flatten(),bins=range(40,101,binsize),density=True)
        #                 n,bins,patches=plt.hist(MonthAnnDailyExt[::,::,ref_lat,ref_lon].flatten(),bins=range(0,81,binsize),density=True)
                        n,bins,patches=plt.hist(bins[:-1],bins,weights=n*binsize)
        
                    plt.plot((bins[:-1]+bins[1:])/2,np.cumsum(n))
                        #n,bins,patches=plt.hist(MonthAnnDailyExt.flatten(),bins=range(6),density=True,cumulative=True,histtype='step')
                    plt.axhline(0.9)
                    if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                        plt.xlabel('Category',fontsize=16)
                    else:
                        plt.xlabel('WBGT (F)',fontsize=16)
                    plt.tick_params(labelsize=16)
                    plt.tight_layout()
                        #plt.title('{} max WBGT distribution'.format(month_labels[i]))
                    if i == 6:
                        plt.savefig('{}/cat_test_max_WBGT_dist_{}_{}_Mean'.format(Save_dir,month_labels[i],Region_Labels[r]),dpi=figdpi)
                    else:
                        plt.savefig('{}/cat_test_max_WBGT_dist_{}_{}_Mean'.format(Save_dir,month_labels[i],Region_Labels[r]))
                    plt.show()
                    plt.close()
                    print(n)
                    print(bins)


                    plt.figure()
                    if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                        n,bins,patches=plt.hist(MonthAnnDailyInt[::,::,minlats[r]:maxlats[r]:,minlons[r]:maxlons[r]:].flatten(),bins=range(40),density=True)
                    else:
                        binsize=5
                        n,bins=np.histogram(MonthAnnDailyInt[::,::,minlats[r]:maxlats[r]:,minlons[r]:maxlons[r]:].flatten(),bins=range(40,101,binsize),density=True)
                        n,bins,patches=plt.hist(bins[:-1],bins,weights=n*binsize)
        
                    plt.plot((bins[:-1]+bins[1:])/2,np.cumsum(n))
                    #n,bins,patches=plt.hist(MonthAnnDailyExt.flatten(),bins=range(6),density=True,cumulative=True,histtype='step')
                    plt.axhline(0.9)
                    if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                        plt.xlabel('Category-Hours',fontsize=16)
                    else:
                        plt.xlabel('WBGT (F)',fontsize=16)
                    plt.tick_params(labelsize=16)
        
                    plt.tight_layout()
                    #plt.title('{} integrated category distribution'.format(month_labels[i]))
                    if i==6:
                        plt.savefig('{}/cat_test_int_cat_dist_{}_{}_Mean'.format(Save_dir,month_labels[i],Region_Labels[r]),dpi=figdpi)
                    else:
                        plt.savefig('{}/cat_test_int_cat_dist_{}_{}_Mean'.format(Save_dir,month_labels[i],Region_Labels[r]))
                    plt.show()
                    plt.close()

                    #put maps here np.nanmean(MonthAnnDailyInt[::,::,minlats[r]:maxlats[r]:,minlons[r]:maxlons[r]:],axis=(0,1))
                    #              np.nanstd(MonthAnnDailyInt[::,::,minlats[r]:maxlats[r]:,minlons[r]:maxlons[r]:],axis=(0,1))

                    plt.figure(figsize=pltsize)
                    m=Basemap(resolution='i', projection='merc', llcrnrlat=minlat, urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon, lat_ts=10)
                    try:
                        m.drawcoastlines(linewidth=0.50)
                    except:
                        pass
                    m.drawcountries(linewidth=0.25)
                    m.drawstates(linewidth=0.15)
                    m.drawmeridians(np.arange(0,360,5),labels=[False,False,False,True])
                    m.drawparallels(np.arange(-90,90,5),labels=[True,False,False,False])
                    #plt.xticks(fontsize=16)
                    #plt.yticks(fontsize=16)
                    x,y=m(lon_2D,lat_2D)
                    if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                        levs=np.arange(-0.5,4.51,1)
                        lev_ticks=[0,1,2,3,4]
                        lev_labels=['No risk','Low','Moderate','High','Extreme']
                        exttype='neither'
                    else:
                        levs=np.arange(30.0,100.1,5)
                        exttype='both'
                    #xi,yi=m(lon[ref_lat,ref_lon],lat[ref_lat,ref_lon])
                    cmap=plt.get_cmap('RdYlGn_r')
                    m.contourf(x,y,np.nanmean(MonthAnnDailyExt[::,::,::,::],axis=(0,1)),cmap=cmap,extend=exttype,levels=levs)
                    if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                        cbar=plt.colorbar(ticks=lev_ticks)
                        cbar.ax.set_yticklabels(lev_labels,rotation=90,va='center',fontsize=16)
                        cbar.set_label("WBGT Category", fontsize=20)
                        #plt.title('{} 90% Max \nWBGT Category'.format(month_labels[i]),fontsize=24)
                    else:
                        cbar=plt.colorbar()
                        cbar.set_label("WBGT (F)",fontsize=20)
                        #plt.title('{} 90% Max \nWBGT'.format(month_labels[i]),fontsize=24)
                    cbar.ax.tick_params(labelsize=16)
                    CS=m.contour(x,y,np.nanstd(MonthAnnDailyExt[::,::,::,::],axis=(0,1)),colors='k')
                    plt.clabel(CS,CS.levels,inline=True,fontsize=12,colors='k')
                    plt.tight_layout()
                    if i==6:
                        plt.savefig('{}/cat_test_maxcat_map_meanstd_{}'.format(Save_dir,month_labels[i]),dpi=figdpi)
                    else:
                        plt.savefig('{}/cat_test_maxcat_map_meanstd_{}'.format(Save_dir,month_labels[i]))
                    plt.show()
                    plt.close()
    
                    plt.figure(figsize=pltsize)
                    m=Basemap(resolution='i', projection='merc', llcrnrlat=minlat, urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon, lat_ts=10)
                    try:
                        m.drawcoastlines(linewidth=0.50)
                    except:
                        pass
                    m.drawcountries(linewidth=0.25)
                    m.drawstates(linewidth=0.15)
                    m.drawmeridians(np.arange(0,360,5),labels=[False,False,False,True])
                    m.drawparallels(np.arange(-90,90,5),labels=[True,False,False,False])
                    #plt.xticks(fontsize=16)
                    #plt.yticks(fontsize=16)
                    x,y=m(lon_2D,lat_2D)
                    if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                        levs=np.arange(-0.5,4.51,1)
                        lev_ticks=[0,1,2,3,4]
                        lev_labels=['No risk','Low','Moderate','High','Extreme']
                        exttype='neither'
                    else:
                        levs=np.arange(30.0,100.1,5)
                        exttype='both'
                    #xi,yi=m(lon[ref_lat,ref_lon],lat[ref_lat,ref_lon])
                    cmap=plt.get_cmap('RdYlGn_r')
                    m.contourf(x,y,np.nanmean(MonthAnnDailyInt[::,::,::,::],axis=(0,1)),cmap=cmap,extend=exttype,levels=levs)
                    if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                        cbar=plt.colorbar(ticks=lev_ticks)
                        cbar.ax.set_yticklabels(lev_labels,rotation=90,va='center',fontsize=16)
                        cbar.set_label("WBGT Category", fontsize=20)
                        #plt.title('{} 90% Max \nWBGT Category'.format(month_labels[i]),fontsize=24)
                    else:
                        cbar=plt.colorbar()
                        cbar.set_label("WBGT (F)",fontsize=20)
                        #plt.title('{} 90% Max \nWBGT'.format(month_labels[i]),fontsize=24)
                    cbar.ax.tick_params(labelsize=16)
                    CS=m.contour(x,y,np.nanstd(MonthAnnDailyInt[::,::,::,::],axis=(0,1)),colors='k')
                    plt.clabel(CS,CS.levels,inline=True,fontsize=12,colors='k')
                    plt.tight_layout()
                    if i==6:
                        plt.savefig('{}/cat_test_intcat_map_meanstd_{}'.format(Save_dir,month_labels[i]),dpi=figdpi)
                    else:
                        plt.savefig('{}/cat_test_intcat_map_meanstd_{}'.format(Save_dir,month_labels[i]))
                    plt.show()
                    plt.close()

            #if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
            plt.figure()
            if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                n,bins,patches=plt.hist(MonthAnnDailyExt[::,::,ref_lat,ref_lon].flatten(),bins=range(6),density=True)
            else:
                binsize=5
                n,bins=np.histogram(MonthAnnDailyExt[::,::,ref_lat,ref_lon].flatten(),bins=range(40,101,binsize),density=True)
#                 n,bins,patches=plt.hist(MonthAnnDailyExt[::,::,ref_lat,ref_lon].flatten(),bins=range(0,81,binsize),density=True)
                n,bins,patches=plt.hist(bins[:-1],bins,weights=n*binsize)

            plt.plot((bins[:-1]+bins[1:])/2,np.cumsum(n))
                #n,bins,patches=plt.hist(MonthAnnDailyExt.flatten(),bins=range(6),density=True,cumulative=True,histtype='step')
            plt.axhline(0.9)
            if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                plt.xlabel('Category',fontsize=16)
            else:
                plt.xlabel('WBGT (F)',fontsize=16)
            plt.tick_params(labelsize=16)
            plt.tight_layout()
                #plt.title('{} max WBGT distribution'.format(month_labels[i]))
            if i == 6:
                plt.savefig('{}/cat_test_max_WBGT_dist_{}_{}_{}'.format(Save_dir,month_labels[i],str(ref_lat),str(ref_lon)),dpi=figdpi)
            else:
                plt.savefig('{}/cat_test_max_WBGT_dist_{}_{}_{}'.format(Save_dir,month_labels[i],str(ref_lat),str(ref_lon)))
            plt.show()
            plt.close()
            print(n)
            print(bins)

            #if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
            plt.figure()
            if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                n,bins,patches=plt.hist(MonthAnnDailyInt[::,::,ref_lat,ref_lon].flatten(),bins=range(40),density=True)
            else:
                binsize=5
                n,bins=np.histogram(MonthAnnDailyInt[::,::,ref_lat,ref_lon].flatten(),bins=range(40,101,binsize),density=True)
                n,bins,patches=plt.hist(bins[:-1],bins,weights=n*binsize)

            plt.plot((bins[:-1]+bins[1:])/2,np.cumsum(n))
            #n,bins,patches=plt.hist(MonthAnnDailyExt.flatten(),bins=range(6),density=True,cumulative=True,histtype='step')
            plt.axhline(0.9)
            if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                plt.xlabel('Category-Hours',fontsize=16)
            else:
                plt.xlabel('WBGT (F)',fontsize=16)
            plt.tick_params(labelsize=16)

            plt.tight_layout()
            #plt.title('{} integrated category distribution'.format(month_labels[i]))
            if i==6:
                plt.savefig('{}/cat_test_int_cat_dist_{}_{}_{}'.format(Save_dir,month_labels[i],str(ref_lat),str(ref_lon)),dpi=figdpi)
            else:
                plt.savefig('{}/cat_test_int_cat_dist_{}_{}_{}'.format(Save_dir,month_labels[i],str(ref_lat),str(ref_lon)))
            plt.show()
            plt.close()

            print(lon_1D.shape)
            print(lat_1D.shape)
      
#            if i==0:
#                lon_2D,lat_2D=np.meshgrid(lon_1D,lat_1D)

            print(lon_2D.shape)
            print(lat_2D.shape)
            print(MonthDailyMax_thresh[month_start_days[i],::,::].shape)


            plt.figure(figsize=pltsize)
            m=Basemap(resolution='i', projection='merc', llcrnrlat=minlat, urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon, lat_ts=10)
            try:
                m.drawcoastlines(linewidth=0.50)
            except:
                pass
            m.drawcountries(linewidth=0.25)
            m.drawstates(linewidth=0.15)
            m.drawmeridians(np.arange(0,360,5),labels=[False,False,False,True])
            m.drawparallels(np.arange(-90,90,5),labels=[True,False,False,False])
            #plt.xticks(fontsize=16)
            #plt.yticks(fontsize=16)
            x,y=m(lon_2D,lat_2D)
            if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                levs=np.arange(-0.5,4.51,1)
                lev_ticks=[0,1,2,3,4]
                lev_labels=['No risk','Low','Moderate','High','Extreme']
                exttype='neither'
            else:
                levs=np.arange(30.0,100.1,5)
                exttype='both'            
            #xi,yi=m(lon[ref_lat,ref_lon],lat[ref_lat,ref_lon])
            cmap=plt.get_cmap('RdYlGn_r')
            m.contourf(x,y,MonthDailyMax_thresh[month_start_days[i],::,::],cmap=cmap,extend=exttype,levels=levs)
            if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                cbar=plt.colorbar(ticks=lev_ticks)
                cbar.ax.set_yticklabels(lev_labels,rotation=90,va='center',fontsize=16)
                cbar.set_label("WBGT Category", fontsize=20)
                #plt.title('{} 90% Max \nWBGT Category'.format(month_labels[i]),fontsize=24)
            else:
                cbar=plt.colorbar()
                cbar.set_label("WBGT (F)",fontsize=20)
                #plt.title('{} 90% Max \nWBGT'.format(month_labels[i]),fontsize=24)
            cbar.ax.tick_params(labelsize=16)
            plt.tight_layout()
            if i==6:
                plt.savefig('{}/cat_test_maxcat_map_{}'.format(Save_dir,month_labels[i]),dpi=figdpi)
            else:
                plt.savefig('{}/cat_test_maxcat_map_{}'.format(Save_dir,month_labels[i]))
            plt.show()
            plt.close()

            #if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
            plt.figure(figsize=pltsize)
            m=Basemap(resolution='i', projection='merc', llcrnrlat=minlat, urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon, lat_ts=10)
            try:
                m.drawcoastlines(linewidth=0.50)
            except:
                pass
            m.drawcountries(linewidth=0.25)
            m.drawstates(linewidth=0.15)
            m.drawmeridians(np.arange(0,360,5),labels=[False,False,False,True])
            m.drawparallels(np.arange(-90,90,5),labels=[True,False,False,False])
            x,y=m(lon_2D,lat_2D)
            #xi,yi=m(lon[ref_lat,ref_lon],lat[ref_lat,ref_lon])
            cmap=plt.get_cmap('RdYlGn_r')
            if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                levs=np.arange(0,35.01,1)
                exttype='max'
            else:
                levs=np.arange(30.0,100.1,5)
                exttype='both'
            m.contourf(x,y,MonthDailyInt_thresh[month_start_days[i],::,::],cmap=cmap,extend=exttype,levels=levs)

            if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                cbar=plt.colorbar()
                cbar.set_label("WBGT Category", fontsize=20)
                #plt.title('{} 90% Max \nWBGT Category'.format(month_labels[i]),fontsize=24)
            else:
                cbar=plt.colorbar()
                cbar.set_label("WBGT (F)",fontsize=20)

            cbar.ax.tick_params(labelsize=16)
            plt.tight_layout()

#            if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
#                plt.title('{} 90th Percentile \nIntegrated WBGT Category'.format(month_labels[i]))
#            else:
#                plt.title('{} 90th Percentile Mean WBGT'.format(month_labels[i]))
            if i==6:
                plt.savefig('{}/cat_test_intcat_map_{}'.format(Save_dir,month_labels[i]),dpi=figdpi)
            else:
                plt.savefig('{}/cat_test_intcat_map_{}'.format(Save_dir,month_labels[i]))

            plt.show()
            plt.close()

            #Xi,Yi=m(Xi,Yi)
#    HWDaysInt=CalcHeatDays(AnnDailyInt,MonthDailyInt_thresh)
    if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
        MonthDailyMax_thresh[np.where(MonthDailyMax_thresh<1)]=1
        MonthDailyInt_thresh[np.where(MonthDailyInt_thresh<1)]=1

    print('Calculating heat wave days max')
#    HWDaysMax=CalcHeatDays(AnnDailyExt[1],MonthDailyMax_thresh,MaxBreak=HWMaxBreak)
    HWDaysMax=CalcHeatDays(AnnDailyExt[1],Daily_max_thresh,MaxBreak=HWMaxBreak)

    del AnnDailyExt

    print(HWDaysMax.shape)
    #for i in range(num_days):
    #    print(i,HWDaysMax[::,i,ref_lat,ref_lon])

    print('Calculating heat wave days int')
 #   print((HWDaysMax[::,::,ref_lat,ref_lon]==2).sum(axis=1))
#    HWDaysInt=CalcHeatDays(AnnDailyInt,MonthDailyInt_thresh,MaxBreak=HWMaxBreak)

    #if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):    
    HWDaysInt=CalcHeatDays(AnnDailyInt,Daily_int_thresh,MaxBreak=HWMaxBreak)

    del AnnDailyInt 

    print(HWDaysInt.shape)
    #for i in range(num_days):
    #    print(i,HWDaysInt[::,i,ref_lat,ref_lon])
#    print((HWDaysInt[::,::,ref_lat,ref_lon]==2).sum(axis=1))

    print('Thresholds: Max then Int')
#    print(MonthDailyMax_thresh[::,ref_lat,ref_lon])
#    print(MonthDailyInt_thresh[::,ref_lat,ref_lon])
    print(Daily_max_thresh[::,ref_lat,ref_lon])
    print(Daily_int_thresh[::,ref_lat,ref_lon])

    del Daily_max_thresh
    del Daily_int_thresh
    del MonthDailyMax_thresh
    del MonthDailyInt_thresh

    print('Heat Waves')
    Num_HW_Max=(HWDaysMax==2).sum(axis=1)
    #if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
    Num_HW_Int=(HWDaysInt==2).sum(axis=1)

    print('Max:',Num_HW_Max[::,ref_lat,ref_lon])
    #if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
    print('Int:',Num_HW_Int[::,ref_lat,ref_lon])

    print('Heat Wave Days')
#    print((HWDaysMax[::,::,ref_lat,ref_lon]>=1).sum(axis=1))
#    print((HWDaysInt[::,::,ref_lat,ref_lon]>=1).sum(axis=1))

    HW_Days_Max_year=(HWDaysMax>=1).sum(axis=1)
    #if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
    HW_Days_Int_year=(HWDaysInt>=1).sum(axis=1)
    print('Max:',HW_Days_Max_year[::,ref_lat,ref_lon])
    #if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
    print('Int:',HW_Days_Int_year[::,ref_lat,ref_lon])

#    print((HWDaysMax[::,::,ref_lat,ref_lon]>=1).sum(axis=0))
#    print((HWDaysInt[::,::,ref_lat,ref_lon]>=1).sum(axis=0))

    HW_Days_Max_day=(HWDaysMax>=1).sum(axis=0)
    #if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
    HW_Days_Int_day=(HWDaysInt>=1).sum(axis=0)

    print('Max:',HW_Days_Max_day[::,ref_lat,ref_lon])
    #if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
    print('Int:',HW_Days_Int_day[::,ref_lat,ref_lon])

    print('Heat Wave Lengths')
    
    print('Max:',HW_Days_Max_year[::,ref_lat,ref_lon]/Num_HW_Max[::,ref_lat,ref_lon])
    #if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
    print('Int:',HW_Days_Int_year[::,ref_lat,ref_lon]/Num_HW_Int[::,ref_lat,ref_lon])

    print('Max:',np.nansum(HW_Days_Max_year[::,ref_lat,ref_lon])/np.nansum(Num_HW_Max[::,ref_lat,ref_lon]))
    #if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
    print('Int:',np.nansum(HW_Days_Int_year[::,ref_lat,ref_lon])/np.nansum(Num_HW_Int[::,ref_lat,ref_lon]))

    #################
    # Plot Mean heat wave days per year
    # Plot Mean heat waves per year
    # Plot Mean heat wave length

    # Calculate and plot trend for each of the above

    #Insert Sigtest here

    # Heat Wave Days 
    HW_Days_Max_year.shape=(num_years,num_lat*num_lon)
    HW_Days_Int_year.shape=(num_years,num_lat*num_lon)

    E,xhat=LinRegress([np.arange(Start_year,End_year+1)],HW_Days_Max_year)
    HW_Days_Max_trend=xhat[0]

    P_max_days,P_tf_max_days=SigTrend(HW_Days_Max_year,HW_Days_Max_trend,5000)

    E,xhat=LinRegress([np.arange(Start_year,End_year+1)],HW_Days_Int_year)
    HW_Days_Int_trend=xhat[0]

    HW_trend=np.dot(E,xhat)

    P_int_days,P_tf_int_days=SigTrend(HW_Days_Int_year,HW_Days_Int_trend,5000)



    HW_Days_Max_year.shape=(num_years,num_lat,num_lon)
    HW_Days_Max_trend.shape=(num_lat,num_lon)
    P_max_days.shape=(num_lat,num_lon)
    P_tf_max_days.shape=(num_lat,num_lon)
    HW_Days_Int_year.shape=(num_years,num_lat,num_lon)
    HW_Days_Int_trend.shape=(num_lat,num_lon)
    P_int_days.shape=(num_lat,num_lon)
    P_tf_int_days.shape=(num_lat,num_lon)

    HW_trend.shape=(num_years,num_lat,num_lon) #??
         
    print('testing Monte Carlo heat wave days max trend')
    print(P_max_days)
    print(P_tf_max_days)    

    print(HW_Days_Int_trend[ref_lat,ref_lon])
    plt.figure()
    plt.plot(np.arange(Start_year,End_year+1),HW_Days_Max_year[::,ref_lat,ref_lon])
    plt.plot(np.arange(Start_year,End_year+1),HW_trend[::,ref_lat,ref_lon])
    plt.show()
    plt.close()

    #HW Days Max
    del E
    del xhat

    levs=np.arange(15,40.01,2.0)
    exttype='max'
    plttitle='Heat wave days'
    units='Days per year'
    pltsave='{}/HW_Days_Max.png'.format(Save_dir)
    Data=np.nanmean(HW_Days_Max_year,axis=0)
    PlotMap(Data,lon_2D,lat_2D,levs,exttype,plttitle,pltsave,colormap='YlOrRd',units=units)

    #HW Days Int

    levs=np.arange(15,40.01,2.0)
    exttype='max'
    plttitle='Heat wave days'
    units='Days per year'
    pltsave='{}/HW_Days_Int.png'.format(Save_dir)
    Data=np.nanmean(HW_Days_Int_year,axis=0)
    PlotMap(Data,lon_2D,lat_2D,levs,exttype,plttitle,pltsave,colormap='YlOrRd',units=units)

    #HW Days Max trend

    levs=np.arange(-0.75*trend_scale,0.7501*trend_scale,.05*trend_scale)
    exttype='both'
    plttitle='HW days trend'
    units=r'$\Delta$'+'Days/year per {} years'.format(str(trend_scale))
    pltsave='{}/HW_Days_Max_trend.png'.format(Save_dir)
    Data=HW_Days_Max_trend*trend_scale
    PlotMap(Data,lon_2D,lat_2D,levs,exttype,plttitle,pltsave,units=units,UseSig=True,Sig=P_tf_max_days)

    #HW Days Int

    levs=np.arange(-0.75*trend_scale,0.7501*trend_scale,0.05*trend_scale)
    exttype='both'
    plttitle='HW days trend'
    units=r'$\Delta$'+'Days/year per {} years'.format(str(trend_scale))
    pltsave='{}/HW_Days_Int_trend.png'.format(Save_dir)
    Data=HW_Days_Int_trend*trend_scale
    PlotMap(Data,lon_2D,lat_2D,levs,exttype,plttitle,pltsave,units=units,UseSig=True,Sig=P_tf_int_days)

    del HW_Days_Max_trend
    del HW_Days_Int_trend

   # Heat Waves

    Num_HW_Max.shape=(num_years,num_lat*num_lon)
    Num_HW_Int.shape=(num_years,num_lat*num_lon)


    E,xhat=LinRegress([np.arange(Start_year,End_year+1)],Num_HW_Max)
    Num_HW_Max_trend=xhat[0]

    P_max_num,P_tf_max_num=SigTrend(Num_HW_Max,Num_HW_Max_trend,5000)

    E,xhat=LinRegress([np.arange(Start_year,End_year+1)],Num_HW_Int)
    Num_HW_Int_trend=xhat[0]

    P_int_num,P_tf_int_num=SigTrend(Num_HW_Int,Num_HW_Int_trend,5000)

    Num_HW_Max.shape=(num_years,num_lat,num_lon)
    Num_HW_Max_trend.shape=(num_lat,num_lon)
    P_max_num.shape=(num_lat,num_lon)
    P_tf_max_num.shape=(num_lat,num_lon)
    Num_HW_Int.shape=(num_years,num_lat,num_lon)
    Num_HW_Int_trend.shape=(num_lat,num_lon)
    P_int_num.shape=(num_lat,num_lon)
    P_tf_int_num.shape=(num_lat,num_lon)

    del E
    del xhat

    #HW Max

    levs=np.arange(5,13.01,0.5)
    exttype='max'
    plttitle='Heat waves'
    units='Heat waves per year'
    pltsave='{}/HW_Max.png'.format(Save_dir)
    Data=np.nanmean(Num_HW_Max,axis=0)
    PlotMap(Data,lon_2D,lat_2D,levs,exttype,plttitle,pltsave,colormap='YlOrRd',units=units)

    #HW Int

    levs=np.arange(5,13.01,0.5)
    exttype='max'
    plttitle='Heat waves'
    units='Heat waves per year'
    pltsave='{}/HW_Int.png'.format(Save_dir)
    Data=np.nanmean(Num_HW_Int,axis=0)
    PlotMap(Data,lon_2D,lat_2D,levs,exttype,plttitle,pltsave,colormap='YlOrRd',units=units)

    #HW Max trend

    levs=np.arange(-0.20*trend_scale,0.2001*trend_scale,0.01*trend_scale)
    exttype='both'
    plttitle='HW trend'
    units=r'$\Delta$'+'HWs/year per {} years'.format(str(trend_scale))
    pltsave='{}/HW_Max_trend.png'.format(Save_dir)
    Data=Num_HW_Max_trend*trend_scale
    PlotMap(Data,lon_2D,lat_2D,levs,exttype,plttitle,pltsave,units=units,UseSig=True,Sig=P_tf_max_num)

    #HW Int trend

    levs=np.arange(-0.20*trend_scale,0.2001*trend_scale,0.01*trend_scale)
    exttype='both'
    plttitle='HW trend'
    units=r'$\Delta$'+'HWs/year per {} years'.format(str(trend_scale))
    pltsave='{}/HW_Int_trend.png'.format(Save_dir)
    Data=Num_HW_Int_trend*trend_scale
    PlotMap(Data,lon_2D,lat_2D,levs,exttype,plttitle,pltsave,units=units,UseSig=True,Sig=P_tf_int_num)

    del Num_HW_Max_trend
    del Num_HW_Int_trend

   # Heat Waves
    HW_Days_Max_year.shape=(num_years,num_lat*num_lon)
    HW_Days_Int_year.shape=(num_years,num_lat*num_lon)
    Num_HW_Max.shape=(num_years,num_lat*num_lon)
    Num_HW_Int.shape=(num_years,num_lat*num_lon)

    HW_len_Max=HW_Days_Max_year/Num_HW_Max
    HW_len_Int=HW_Days_Int_year/Num_HW_Int
    HW_len_Max[np.where(np.isnan(HW_len_Max))]=0
    HW_len_Int[np.where(np.isnan(HW_len_Int))]=0

    E,xhat=LinRegress([np.arange(Start_year,End_year+1)],HW_len_Max)
    HW_len_Max_trend=xhat[0]

    P_max_len,P_tf_max_len=SigTrend(HW_len_Max,HW_len_Max_trend,5000)

    E,xhat=LinRegress([np.arange(Start_year,End_year+1)],HW_len_Int)
    HW_len_Int_trend=xhat[0]

    P_int_len,P_tf_int_len=SigTrend(HW_len_Int,HW_len_Int_trend,5000)

    HW_len_Max_trend.shape=(num_lat,num_lon)
    HW_len_Int_trend.shape=(num_lat,num_lon)
    HW_Days_Max_year.shape=(num_years,num_lat,num_lon)
    HW_Days_Int_year.shape=(num_years,num_lat,num_lon)
    Num_HW_Max.shape=(num_years,num_lat,num_lon)
    Num_HW_Int.shape=(num_years,num_lat,num_lon)
    HW_len_Max.shape=(num_years,num_lat,num_lon)
    HW_len_Int.shape=(num_years,num_lat,num_lon)
    P_max_len.shape=(num_lat,num_lon)
    P_tf_max_len.shape=(num_lat,num_lon)
    P_int_len.shape=(num_lat,num_lon)
    P_tf_int_len.shape=(num_lat,num_lon)

    del E
    del xhat

    #HW length Max

    levs=np.arange(2,5.01,0.25)
    exttype='max'
    plttitle='Heat wave length'
    units='Days'
    pltsave='{}/HW_len_Max.png'.format(Save_dir)
    Data=np.nanmean(HW_Days_Max_year,axis=0)/np.nanmean(Num_HW_Max,axis=0)
    PlotMap(Data,lon_2D,lat_2D,levs,exttype,plttitle,pltsave,colormap='YlOrRd',units=units)

    #HW length Int

    levs=np.arange(2,5.01,0.25)
    exttype='max'
    plttitle='Heat wave length'
    units='Days'
    pltsave='{}/HW_len_Int.png'.format(Save_dir)
    Data=np.nanmean(HW_Days_Int_year,axis=0)/np.nanmean(Num_HW_Int,axis=0)
    PlotMap(Data,lon_2D,lat_2D,levs,exttype,plttitle,pltsave,colormap='YlOrRd',units=units)

    #HW length Max trend

    levs=np.arange(-0.03*trend_scale,0.0301*trend_scale,0.005*trend_scale)
    exttype='both'
    plttitle='HW length trend'
    units=r'$\Delta$'+'Days/HW per {} years'.format(str(trend_scale))
    pltsave='{}/HW_len_Max_trend.png'.format(Save_dir)
    Data=HW_len_Max_trend*trend_scale
    PlotMap(Data,lon_2D,lat_2D,levs,exttype,plttitle,pltsave,units=units,UseSig=True,Sig=P_tf_max_len)

    #HW length Int trend

    levs=np.arange(-0.03*trend_scale,0.0301*trend_scale,0.005*trend_scale)
    exttype='both'
    plttitle='HW length trend'
    units=r'$\Delta$'+'Days/HW per {} years'.format(str(trend_scale))
    pltsave='{}/HW_len_Int_trend.png'.format(Save_dir)
    Data=HW_len_Int_trend*trend_scale
    PlotMap(Data,lon_2D,lat_2D,levs,exttype,plttitle,pltsave,units=units,UseSig=True,Sig=P_tf_int_len)

    del HW_len_Max_trend
    del HW_len_Int_trend

    SaveHW=True
    if (SaveHW) and (WBGT_type=='WBGT'):
        print(HWDaysMax.shape)
        HWDaysMax.shape=(num_years*num_days,num_lat,num_lon)

        Input_path=ERA5_path
        lon=lon[:]
        lat=lat[:]
        #year=np.arange(Start_year,End_year+1)[:]
        #day=np.arange(num_years)[:] 
        time=np.arange(num_years*num_days)[:]
        Dims=[lon,lat,time]
        DimsName=['longitude','latitude','time']
        DimsUnits=['degrees_east','degrees_north','days_since_1-1-1960']
        DimsLN=['longitude','latitude','time']
        Vars=[HWDaysMax]
        VarsNames=['HWDaysMax']
        VarsUnits=['HW']
        VarsLN=['Heat Wave Days Max WBGT']
        VarsDims=[('time','latitude','longitude')]
        output_dir='{}/{}'.format(Input_path,'HWDays')
        filename='{}.nc'.format('HWDaysMax')

        WriteNetCDF(Dims,DimsName,DimsUnits,DimsLN,Vars,VarsNames,VarsUnits,VarsLN,VarsDims,output_dir,filename)

        print(HWDaysInt.shape)
        HWDaysInt.shape=(num_years*num_days,num_lat,num_lon)

        Input_path=ERA5_path
        lon=lon[:]
        lat=lat[:]
        #year=np.arange(Start_year,End_year+1)[:]
        #day=np.arange(num_years)[:] 
        time=np.arange(num_years*num_days)[:]
        Dims=[lon,lat,time]
        DimsName=['longitude','latitude','time']
        DimsUnits=['degrees_east','degrees_north','days_since_1-1-1960']
        DimsLN=['longitude','latitude','time']
        Vars=[HWDaysInt]
        VarsNames=['HWDaysInt']
        VarsUnits=['HW']
        VarsLN=['Heat Wave Days Int WBGT']
        VarsDims=[('time','latitude','longitude')]
        output_dir='{}/{}'.format(Input_path,'HWDays')
        filename='{}.nc'.format('HWDaysInt')

        WriteNetCDF(Dims,DimsName,DimsUnits,DimsLN,Vars,VarsNames,VarsUnits,VarsLN,VarsDims,output_dir,filename)

#    del AnnDailyCount
#    del AnnDailyExt
#    del AnnDailyInt
#    del indmax
#    del indmin
#    del indc
#    del indcmin
#    del indcmax
#    del DailyCount
#    del DailyMinCount
#    del DailyMaxCount
#    del DailyTrend
#    del E
#    del xhat
#    del Eint
#    del xhat_int
#    del MonthDailyMax_thresh
#    del MonthDailyInt_thresh
#    del Daily_max_thresh
#    del Daily_int_thresh
    del HWDaysMax
    del HWDaysInt
    del Num_HW_Max
    del Num_HW_Int
    del HW_Days_Max_day
    del HW_Days_Int_day
    del HW_len_Max
    del HW_len_Int

#    for i in range(month_start_days[6],month_start_days[9]):
#        print(AnnDailyExt[1,::,i,ref_lat,ref_lon])
#        print(HWDaysMax[::,i,ref_lat,ref_lon])
#        print(Daily_max_thresh[i,ref_lat,ref_lon])




#    print((HWDaysMax[::,::,ref_lat,ref_lon]>=1).sum()/~np.isnan(AnnDailyExt[1,::,::,ref_lat,ref_lon]).sum())
#    print((HWDaysInt[::,::,ref_lat,ref_lon]>=1).sum()/~np.isnan(AnnDailyInt[::,::,ref_lat,ref_lon]).sum())

#    E,xhat=LinRegress([np.arange(Start_year,End_year+1)],AnnMonthlyFullCount[c,::,::])
#    MonthlyFullTrend[c,::]=xhat[0]










#Same frequency as daily climo but shifted by 4 hours so day starts at 20z for purposes of identifying distribution of nighttime lows
if ShiftedDailyClimo:
#sum data by day each year
    print('#####################')
    print('# ShiftedDailyClimo #')
    print('#####################')

    print('creating cats_shifted')

    print(type(Cats))
 #   Catc.astype(np.float32)
    print(type(Cats))
    
    Cats_shifted=np.ones(Cats.shape,dtype=np.float32)*np.nan
    Cats_shifted[::,4::,::,::]=Cats[::,:-4:,::,::]

    del Cats

    print('Done creating cats_shifted')

    if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
        AnnDailyCount=np.zeros((num_cats,num_years,num_days,I,J),dtype=np.float32)
    else:
        AnnDailyExt=np.ones((2,num_years,num_days,I,J),dtype=np.float32)*np.nan

    for y in range(num_years):
        for d in range(num_days):
            if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                for c in range(num_cats):
                    indc=np.nansum((Cats_shifted[y,d*num_hours:(d+1)*num_hours,::,::]==c),axis=0,dtype=np.float32)
                    AnnDailyCount[c,y,d,::,::]=indc
            else:
                indc=np.nanmin(Cats_shifted[y,d*num_hours:(d+1)*num_hours,::,::],axis=0)
                AnnDailyExt[0,y,d,::,::]=indc.astype(np.float32)
        del indc

    del Cats_shifted

    if Pre_calc:
        if WBGT_type=='Cat':
            f='{}/{}/{}CatClimo.nc'.format(ERA5_path,Var_name_path,Var_name_path)
            #Varn, Varnm, lat_pre, lon_pre = ReadNetCDF(f,['time','WBGT10MinRM','WBGT10MinM'])#,chunk_start,chunk_end)
            Varn=xr.load_dataset(f)
            Pre_vars=['WBGT10MinRM','WBGT10MinM']
        elif WBGT_type=='ModCat':
            f='{}/{}/{}ModCatClimo.nc'.format(ERA5_path,Var_name_path,Var_name_path,y)
            #Varn, Varnm, lat_pre, lon_pre = ReadNetCDF(f,['time','WBGT10MinRM','WBGT10MinM'])#,chunk_start,chunk_end)
            Varn=xr.load_dataset(f)
            Pre_vars=['WBGT10MinRM','WBGT10MinM']

        elif WBGT_type=='WBGT':
            f='{}/{}/{}Climo.nc'.format(ERA5_path,Var_name_path,Var_name_path)
            #Varn, Varnm, lat_pre, lon_pre = ReadNetCDF(f,['time','WBGT10RM','WBGT10M'])#,chunk_start,chunk_end)
            Varn=xr.load_dataset(f)
            Pre_vars=['WBGT10RM','WBGT10M']

        Varn_units=[]
        for Pre_var in Pre_vars:
            mv=Varn[Pre_var].encoding['missing_value']
            Varn_units.append(Varn[Pre_var].units)
            Varn[Pre_var]=xr.where(Varn[Pre_var]==mv,np.nan,Varn[Pre_var])

        #time_pre=Varn[0]
        if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
            Daily_min_thresh=np.rint(Varn[Pre_vars[0]]).astype(np.float32)
            MonthDailyCount_thresh=np.rint(Varn[Pre_vars[1]]).astype(np.float32)
        elif WBGT_type=='WBGT':
            Daily_min_thresh=Varn[Pre_vars[0]].astype(np.float32)
            MonthDailyCount_thresh=Varn[Pre_vars[1]].astype(np.float32)
            Daily_min_thresh=Convert_T_Unit(Daily_min_thresh,Varn_units[0],'F')
            MonthDailyCount_thresh=Convert_T_Unit(MonthDailyCount_thresh,Varn_units[1],'F')
    
        del Varn
        #del Varnm
    
    else:
        if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
            MonthDailyCount_thresh=np.ones((num_days,I,J),dtype=np.float32)*np.nan



    lon_2D,lat_2D=np.meshgrid(lon_1D,lat_1D)
    PlotShiftedDailyClimo=True
    for i in range(num_months): #put before if statement so I can add low WBGT
        if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
            MonthAnnDailyCount=AnnDailyCount[0,::,month_start_days[i]:int(expected_ei[i]/24):,::,::]
            print(MonthAnnDailyCount.shape)
            if i==1:
                MonthAnnDailyCount[np.where(MonthAnnDailyCount==4)]=np.nan
            if (i==0) or (i==2):
                MonthAnnDailyCount[np.where(MonthAnnDailyCount==20)]=np.nan
            #MonthAnnDailyCount.shape=((num_years*int(month_start_days[i+1]-month_start_days[i],)))
            #print(MonthAnnDailyCount.flatten())
            if not Pre_calc:
                Y2,M,I2,J2=MonthAnnDailyCount.shape
                MonthAnnDailyCount=MonthAnnDailyCount.reshape(Y*M,I,J)
                MonthDailyCount_thresh[month_start_days[i]:int(expected_ei[i]/24):,::,::]=np.nanpercentile(MonthAnnDailyCount,10,axis=0,interpolation='linear').astype(np.float32)
                MonthDailyCount_thresh[np.where(MonthDailyCount_thresh>23)]=23

                MonthAnnDailyCount=MonthAnnDailyCount.reshape(Y,M,I,J)
        else:
            MonthAnnDailyExt=AnnDailyExt[0,::,month_start_days[i]:int(expected_ei[i]/24):,::,::]
            print(MonthAnnDailyExt.shape)
            

#        if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
        if PlotShiftedDailyClimo:

            if RegionalMean:
                for r in range(IJ_reg):
                    plt.figure()
                    if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                        n,bins,patches=plt.hist(MonthAnnDailyExt[::,::,minlats[r]:maxlats[r]:,minlons[r]:maxlons[r]:].flatten(),bins=range(6),density=True)
                    else:
                        binsize=5
                        n,bins=np.histogram(MonthAnnDailyExt[::,::,minlats[r]:maxlats[r]:,minlons[r]:maxlons[r]:].flatten(),bins=range(40,101,binsize),density=True)
        #                 n,bins,patches=plt.hist(MonthAnnDailyExt[::,::,ref_lat,ref_lon].flatten(),bins=range(0,81,binsize),density=True)
                        n,bins,patches=plt.hist(bins[:-1],bins,weights=n*binsize)

                    plt.plot((bins[:-1]+bins[1:])/2,np.cumsum(n))
                        #n,bins,patches=plt.hist(MonthAnnDailyExt.flatten(),bins=range(6),density=True,cumulative=True,histtype='step')
                    plt.axhline(0.9)
                    if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                        plt.xlabel('Category',fontsize=16)
                    else:
                        plt.xlabel('WBGT (F)',fontsize=16)
                    plt.tick_params(labelsize=16)
                    plt.tight_layout()
                        #plt.title('{} max WBGT distribution'.format(month_labels[i]))
                    if i == 6:
                        plt.savefig('{}/cat_test_min_WBGT_dist_{}_{}_Mean'.format(Save_dir,month_labels[i],Region_Labels[r]),dpi=figdpi)
                    else:
                        plt.savefig('{}/cat_test_min_WBGT_dist_{}_{}_Mean'.format(Save_dir,month_labels[i],Region_Labels[r]))
                    plt.show()
                    plt.close()
                    print(n)
                    print(bins)

                plt.figure(figsize=pltsize)
                m=Basemap(resolution='i', projection='merc', llcrnrlat=minlat, urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon, lat_ts=10)
                try:
                    m.drawcoastlines(linewidth=0.50)
                except:
                    pass
                m.drawcountries(linewidth=0.25)
                m.drawstates(linewidth=0.15)
                m.drawmeridians(np.arange(0,360,5),labels=[False,False,False,True])
                m.drawparallels(np.arange(-90,90,5),labels=[True,False,False,False])
                #plt.xticks(fontsize=16)
                #plt.yticks(fontsize=16)
                x,y=m(lon_2D,lat_2D)
                if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                    levs=np.arange(-0.5,4.51,1)
                    lev_ticks=[0,1,2,3,4]
                    lev_labels=['No risk','Low','Moderate','High','Extreme']
                    exttype='neither'
                else:
                    levs=np.arange(30.0,100.1,5)
                    exttype='both'
                #xi,yi=m(lon[ref_lat,ref_lon],lat[ref_lat,ref_lon])
                cmap=plt.get_cmap('RdYlGn_r')
                m.contourf(x,y,np.nanmean(MonthAnnDailyExt[::,::,::,::],axis=(0,1)),cmap=cmap,extend=exttype,levels=levs)
                if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                    cbar=plt.colorbar(ticks=lev_ticks)
                    cbar.ax.set_yticklabels(lev_labels,rotation=90,va='center',fontsize=16)
                    cbar.set_label("WBGT Category", fontsize=20)
                    #plt.title('{} 90% Max \nWBGT Category'.format(month_labels[i]),fontsize=24)
                else:
                    cbar=plt.colorbar()
                    cbar.set_label("WBGT (F)",fontsize=20)
                    #plt.title('{} 90% Max \nWBGT'.format(month_labels[i]),fontsize=24)
                cbar.ax.tick_params(labelsize=16)
                CS=m.contour(x,y,np.nanstd(MonthAnnDailyExt[::,::,::,::],axis=(0,1)),colors='k')
                plt.clabel(CS,CS.levels,inline=True,fontsize=12,colors='k')
                plt.tight_layout()
                if i==6:
                    plt.savefig('{}/cat_test_mincat_map_meanstd_{}'.format(Save_dir,month_labels[i]),dpi=figdpi)
                else:
                    plt.savefig('{}/cat_test_mincat_map_meanstd_{}'.format(Save_dir,month_labels[i]))
                plt.show()
                plt.close()




            plt.figure()
            if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):            
                n,bins,patches=plt.hist(MonthAnnDailyCount[::,::,ref_lat,ref_lon].flatten(),bins=range(26),density=True)
                plt.plot((bins[:-1]+bins[1:])/2,np.cumsum(n))
            else:
                binsize=5
                n,bins=np.histogram(MonthAnnDailyExt[::,::,ref_lat,ref_lon].flatten(),bins=range(40,101,binsize),density=True)
#                n,bins,patches=plt.hist(MonthAnnDailyExt[::,::,ref_lat,ref_lon].flatten(),bins=range(0,81,binsize),density=True)
                n,bins,patches=plt.hist(bins[:-1],bins,weights=n*binsize)
                plt.plot((bins[:-1]+bins[1:])/2,np.cumsum(n))
            #n,bins,patches=plt.hist(MonthAnnDailyCount.flatten(),bins=range(26),density=True,cumulative=True,histtype='step')

            if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                plt.xlabel('Hours',fontsize=16)
            else:
                plt.xlabel('WBGT (F)',fontsize=16)
            plt.tick_params(labelsize=16)
            if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                plt.axhline(0.1)
                #plt.title('{} Cat 0 distribution'.format(month_labels[i]))
            else: #This is for if I create a distribution of low WBGT
                plt.axhline(0.9)
                #plt.title('{} Min WBGT distribution'.format(month_labels[i]))
            plt.tight_layout()
            if i==6:
                if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                    plt.savefig('{}/cat_test_Cat0_dist_{}_{}_{}'.format(Save_dir,month_labels[i],str(ref_lat),str(ref_lon)),dpi=figdpi)
                else:
                    plt.savefig('{}/cat_test_Min_dist_{}_{}_{}'.format(Save_dir,month_labels[i],str(ref_lat),str(ref_lon)),dpi=figdpi)
            else:
                if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                    plt.savefig('{}/cat_test_Cat0_dist_{}_{}_{}'.format(Save_dir,month_labels[i],str(ref_lat),str(ref_lon)))
                else:
                    plt.savefig('{}/cat_test_Min_dist_{}_{}_{}'.format(Save_dir,month_labels[i],str(ref_lat),str(ref_lon)))
            plt.show()
            plt.close()
            print(n)
            print(bins)
            
            if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                del MonthAnnDailyCount

        if PlotShiftedDailyClimo:
            plt.figure(figsize=pltsize)
            m=Basemap(resolution='i', projection='merc', llcrnrlat=minlat, urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon, lat_ts=10)
            try:
                m.drawcoastlines(linewidth=0.50)
            except:
                pass
            m.drawcountries(linewidth=0.25)
            m.drawstates(linewidth=0.15)
            m.drawmeridians(np.arange(0,360,5),labels=[False,False,False,True])
            m.drawparallels(np.arange(-90,90,5),labels=[True,False,False,False])
            print(lon_2D.shape,lat_2D.shape)
            x,y=m(lon_2D,lat_2D)
            #xi,yi=m(lon[ref_lat,ref_lon],lat[ref_lat,ref_lon])
            if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                levs=np.arange(0,24.1,1)
                cmap=plt.get_cmap('RdYlGn')
            else:
                levs=np.arange(30,100.1,5)
                cmap=plt.get_cmap('RdYlGn_r')
            m.contourf(x,y,MonthDailyCount_thresh[month_start_days[i],::,::],cmap=cmap,extend='max',levels=levs)

            if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
                cbar=plt.colorbar()
                cbar.set_label("No Risk Hours", fontsize=20)
                #plt.title('{} 90% Min \nWBGT Category'.format(month_labels[i]),fontsize=24)
            else:
                cbar=plt.colorbar()
                cbar.set_label("WBGT (F)",fontsize=20)

#            if (WBGT_type=='Cat') or (WBGT_type=='ModCat'):
#                plt.title('10th Percentile Category 0 For Month {}'.format(str(i)))
#            else:
#                plt.title('90th Percentile Minimum WBGT For Month {}'.format(str(i)))
            cbar.ax.tick_params(labelsize=16)
            plt.tight_layout()

            if i==6:
                plt.savefig('{}/cat_test_mincat_map_{}'.format(Save_dir,month_labels[i]),dpi=1000)
            else:
                plt.savefig('{}/cat_test_mincat_map_{}'.format(Save_dir,month_labels[i]))
            plt.show()
            plt.close()

    if (WBGT_type=='Cat') or (WBGT_type=='ModCat'): 
        HWDays0=CalcHeatDays(AnnDailyCount[0],Daily_min_thresh,'Below',MaxBreak=HWMaxBreak)
        del AnnDailyCount
    else:
        HWDays0=CalcHeatDays(AnnDailyExt[0],Daily_min_thresh,'Above',MaxBreak=HWMaxBreak)
        del AnnDailyExt
    print(HWDays0.shape)
    #for i in range(num_days):
    #    print(i,HWDays0[::,i,ref_lat,ref_lon])
    print((HWDays0[::,::,ref_lat,ref_lon]==2).sum(axis=1))

    print('Daily Threshold')
    print(Daily_min_thresh[::,ref_lat,ref_lon])
    print('Monthly Threshold')
    print(MonthDailyCount_thresh[::,ref_lat,ref_lon])

    del Daily_min_thresh
    del MonthDailyCount_thresh

    print('Heat Waves')
    Num_HW_0=(HWDays0==2).sum(axis=1)
#    Num_HW_Int=(HWDaysInt[::,::,ref_lat,ref_lon]==2).sum(axis=1)

    print('Low:',Num_HW_0[::,ref_lat,ref_lon])
#    print('Int:',Num_HW_Int)

    print('Heat Wave Days')
#    print((HWDaysMax[::,::,ref_lat,ref_lon]>=1).sum(axis=1))
#    print((HWDaysInt[::,::,ref_lat,ref_lon]>=1).sum(axis=1))

    HW_Days_0_year=(HWDays0>=1).sum(axis=1)
#    HW_Days_Int_year=(HWDaysInt[::,::,ref_lat,ref_lon]>=1).sum(axis=1)
    print('Low:',HW_Days_0_year[::,ref_lat,ref_lon])
#    print('Int:',HW_Days_Int_year)

#    print((HWDaysMax[::,::,ref_lat,ref_lon]>=1).sum(axis=0))
#    print((HWDaysInt[::,::,ref_lat,ref_lon]>=1).sum(axis=0))

    HW_Days_0_day=(HWDays0>=1).sum(axis=0)
#    HW_Days_Int_day=(HWDaysInt[::,::,ref_lat,ref_lon]>=1).sum(axis=0)

    print('Low:',HW_Days_0_day[::,ref_lat,ref_lon])
#    print('Int:',HW_Days_Int_day)

#    print(MonthDailyInt_thresh[::,ref_lat,ref_lon])
#    print(Daily_max_thresh[::,ref_lat,ref_lon])
#    print(Daily_int_thresh[::,ref_lat,ref_lon])

    print('Heat Wave Lengths')

    print('Low:',HW_Days_0_year[::,ref_lat,ref_lon]/Num_HW_0[::,ref_lat,ref_lon])
#    print('Int:',HW_Days_Int_year/Num_HW_Int)

    print('Low:',np.nansum(HW_Days_0_year[::,ref_lat,ref_lon])/np.nansum(Num_HW_0[::,ref_lat,ref_lon]))
#    print('Int:',np.nansum(HW_Days_Int_year)/np.nansum(Num_HW_Int))

    # Heat Wave Days
    HW_Days_0_year.shape=(num_years,num_lat*num_lon)

    E,xhat=LinRegress([np.arange(Start_year,End_year+1)],HW_Days_0_year)
    HW_Days_0_trend=xhat[0]

    P_min_days,P_tf_min_days=SigTrend(HW_Days_0_year,HW_Days_0_trend,5000)

    HW_Days_0_year.shape=(num_years,num_lat,num_lon)
    HW_Days_0_trend.shape=(num_lat,num_lon)
    P_min_days.shape=(num_lat,num_lon)
    P_tf_min_days.shape=(num_lat,num_lon)
    #HW Days Max
    del E
    del xhat

    levs=np.arange(15.0,40.01,2.0)
    exttype='max'
    plttitle='Heat wave days'
    units='Days per year'
    pltsave='{}/HW_Days_Min.png'.format(Save_dir)
    Data=np.nanmean(HW_Days_0_year,axis=0)
    PlotMap(Data,lon_2D,lat_2D,levs,exttype,plttitle,pltsave,colormap='YlOrRd',units=units)

    #HW Days Max trend

    levs=np.arange(-0.75*trend_scale,0.7501*trend_scale,0.05*trend_scale)
    exttype='both'
    plttitle='HW days/year trend'
    units=r'$\Delta$'+'Days/year per {} years'.format(str(trend_scale))
    pltsave='{}/HW_Days_Min_trend.png'.format(Save_dir)
    Data=HW_Days_0_trend*trend_scale
    PlotMap(Data,lon_2D,lat_2D,levs,exttype,plttitle,pltsave,units=units,UseSig=True,Sig=P_tf_min_days)

    del HW_Days_0_trend

   # Heat Waves
    Num_HW_0.shape=(num_years,num_lat*num_lon)

    E,xhat=LinRegress([np.arange(Start_year,End_year+1)],Num_HW_0)
    Num_HW_0_trend=xhat[0]

    P_min_num,P_tf_min_num=SigTrend(Num_HW_0,Num_HW_0_trend,5000)

    Num_HW_0.shape=(num_years,num_lat,num_lon)
    Num_HW_0_trend.shape=(num_lat,num_lon)
    P_min_num.shape=(num_lat,num_lon)
    P_tf_min_num.shape=(num_lat,num_lon)
    del E
    del xhat

    #HW Max

    levs=np.arange(5,13.01,0.5)
    exttype='max'
    plttitle='Heat waves'
    units='Heat waves per year'
    pltsave='{}/HW_Min.png'.format(Save_dir)
    Data=np.nanmean(Num_HW_0,axis=0)
    PlotMap(Data,lon_2D,lat_2D,levs,exttype,plttitle,pltsave,colormap='YlOrRd',units=units)

    #HW Max trend

    levs=np.arange(-0.20*trend_scale,0.2001*trend_scale,0.01*trend_scale)
    exttype='both'
    plttitle='HW trend'
    units=r'$\Delta$'+'HWs/year per {} years'.format(str(trend_scale))
    pltsave='{}/HW_Min_trend.png'.format(Save_dir)
    Data=Num_HW_0_trend*trend_scale
    PlotMap(Data,lon_2D,lat_2D,levs,exttype,plttitle,pltsave,units=units,UseSig=True,Sig=P_tf_min_num)

    del Num_HW_0_trend

   # Heat Waves
    Num_HW_0.shape=(num_years,num_lat*num_lon)
    HW_Days_0_year.shape=(num_years,num_lat*num_lon)

    HW_len_0=HW_Days_0_year/Num_HW_0
    HW_len_0[np.where(np.isnan(HW_len_0))]=0
    HW_len_0[np.where(np.isinf(HW_len_0))]=0

#    HW_len_Int=HW_Days_Int_year/Num_HW_Int
#    HW_len_Max[np.where(np.isnan(HW_len_Max))]=0

    E,xhat=LinRegress([np.arange(Start_year,End_year+1)],HW_len_0)
    HW_len_0_trend=xhat[0]

    P_min_len,P_tf_min_len=SigTrend(HW_len_0,HW_len_0_trend,5000)

    print('debugging...')
    print(np.where(np.isnan(HW_len_0)))
    print(np.where(np.isnan(HW_len_0_trend)))
    badnans=np.where(np.isnan(HW_len_0_trend))[0]
    for n in badnans:
        print(HW_len_0[::,n])

    HW_len_0_trend.shape=(num_lat,num_lon)
    Num_HW_0.shape=(num_years,num_lat,num_lon)
    HW_Days_0_year.shape=(num_years,num_lat,num_lon)
    HW_len_0.shape=(num_years,num_lat,num_lon)
    P_min_len.shape=(num_lat,num_lon)
    P_tf_min_len.shape=(num_lat,num_lon)

    del E
    del xhat

    #HW length Max

    levs=np.arange(2,5.01,0.25)
    exttype='max'
    plttitle='Heat wave length'
    units='Days'
    pltsave='{}/HW_len_Min.png'.format(Save_dir)
    Data=np.nanmean(HW_Days_0_year,axis=0)/np.nanmean(Num_HW_0,axis=0)
    PlotMap(Data,lon_2D,lat_2D,levs,exttype,plttitle,pltsave,colormap='YlOrRd',units=units)

    #HW length Max trend

    levs=np.arange(-0.03*trend_scale,0.0301*trend_scale,0.005*trend_scale)
    exttype='both'
    plttitle='HW length trend'
    units=r'$\Delta$'+'Days/HW per {} years'.format(str(trend_scale))
    pltsave='{}/HW_len_Min_trend.png'.format(Save_dir)
    Data=HW_len_0_trend*trend_scale
    PlotMap(Data,lon_2D,lat_2D,levs,exttype,plttitle,pltsave,units=units,UseSig=True,Sig=P_tf_min_len)

    del HW_len_0_trend

    SaveHW=True
    if (SaveHW) & (WBGT_type=='WBGT'):
        print(HWDays0.shape)
        HWDays0.shape=(num_years*num_days,num_lat,num_lon)

        Input_path=ERA5_path
        lon=lon[:]
        lat=lat[:]
        #year=np.arange(Start_year,End_year+1)[:]
        #day=np.arange(num_years)[:] 
        time=np.arange(num_years*num_days)[:]
        Dims=[lon,lat,time]
        DimsName=['longitude','latitude','time']
        DimsUnits=['degrees_east','degrees_north','days_since_1-1-1960']
        DimsLN=['longitude','latitude','time']
        Vars=[HWDays0]
        VarsNames=['HWDaysMin']
        VarsUnits=['HW']
        VarsLN=['Heat Wave Days Min WBGT']
        VarsDims=[('time','latitude','longitude')]
        output_dir='{}/{}'.format(Input_path,'HWDays')
        filename='{}.nc'.format('HWDaysMin')

        WriteNetCDF(Dims,DimsName,DimsUnits,DimsLN,Vars,VarsNames,VarsUnits,VarsLN,VarsDims,output_dir,filename)

    del HWDays0
    del Num_HW_0
    del HW_Days_0_year
    del HW_Days_0_day
    del HW_len_0

#sum data by month each year

#calc trends for each month

#sum data for each month

if MonthlyFullClimo:
    AnnMonthlyFullCount=np.zeros((num_cats,num_years,num_months,I,J))
    for y in range(num_years):
        MonthStart=0
        for m in range(num_months):
            MonthEnd=expected_ei[m]
            for c in range(num_cats):
                indc=np.nansum((Cats[y,MonthStart:MonthEnd:,::,::]==c),axis=0)
                #print("indc shape",indc.shape)
                AnnMonthlyFullCount[c,y,m,::,::]=indc
                #MonthlyCount[c,m*num_hours+h,::,::]=np.nansum(HourlyCount[c,MonthStart+h:MonthEnd:24,::,::],axis=0)
            MonthStart=MonthEnd*1

    
    #print(AnnMonthlyCount)
    MonthlyFullCount=np.nansum(AnnMonthlyFullCount,axis=1)
    print(MonthlyFullCount.shape)
    print(MonthlyFullCount[0,::,ref_lat,ref_lon])

    AnnMonthlyFullCount.shape=(num_cats,num_years,num_months,I*J)
    AnnMonthlyFullCount.shape=(num_cats,num_years,num_months*I*J)
    MonthlyFullTrend=np.ones((num_cats,num_months*I*J))*np.nan
    for c in range(num_cats):
        E,xhat=LinRegress([np.arange(Start_year,End_year+1)],AnnMonthlyFullCount[c,::,::])
        MonthlyFullTrend[c,::]=xhat[0]

    AnnMonthlyFullCount.shape=(num_cats,num_years,num_months,I*J)
    MonthlyFullTrend.shape=(num_cats,num_months,I*J)
    AnnMonthlyFullCount.shape=(num_cats,num_years,num_months,I,J)
    MonthlyFullTrend.shape=(num_cats,num_months,I,J)





#sum data by season each year
#calc trends for each season

#sum data for each season


#write data to files for plotting
	#There should be a separate block for each above block
	#Each above block should write to a new file


#######Code from general/heat_class.py to identify heat days#######
#######Modify for this data format
#######HWtemp is threshold variables
#######Full_data is the data array
#######Might not need all the reshaping idk
#######Make sure that "default" category (24 hours, 0 int, 0 cat) isnt considered heat day

  







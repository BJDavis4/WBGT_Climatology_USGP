######################################################################
# This will be very similar to the Plot Climo and Calc climo scripts #
######################################################################
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
from netCDF_mods import ReadNetCDF, WriteNetCDF
import glob
from wet_bulb import Vapr, Vapr_to_T
from Python_mods import Convert_T_Unit

#set parameters
mesonet_path='/data/deluge/observations/OKMesonet/grid'
#mesonet_path='/data/deluge/scratch/ERA5-Ben/2D/hourly'
ERA5_path='/data/deluge/scratch/ERA5-Ben/2D/hourly'
Output_path=ERA5_path
mesonet_vars=['time','TWBG','RELH','TAIR','SRAD','WS2M','PRES','WSPD']
ERA5_vars=['time','WBGT','d2m','t2m','msdwswrf','ws2m','sp','ws10']
var_units_plt=['seconds',r'${\degree}F$',r'${\degree}C$',r'${\degree}C$',r'$W m^\mathregular{-2}$',r'$m s^\mathregular{-1}$',r'$hPa$',r'$m s^\mathregular{-1}$']
mesotype='obs'#['obs','WBGT']
#mesonet, used for plotting
fext='nc'
hours=np.arange(24)
#print(hours)
#exit()

start_month=1
end_month=12
start_year=1960
end_year=2020

Num_years=end_year-start_year+1
Num_days=366
Num_hours=24
Num_lat=11
Num_lon=22

if mesotype=='WBGT':
    mesonet_vars=ERA5_vars[:2]
    ERA5_vars=ERA5_vars[:2]

def gridRMSD(arrin):#,nans,nansplij,nanspli,Y,T,J,I): #This needs changed
    arrIn=np.copy(arrin)
    #Figure out exactly what to do with the next two rows
    #arrIn[::,::,nans]=np.nan
    #arrIn[::,::,nansplj,nanspli]=np.nan
    

    arrIn=np.square(arrIn)
    arrCount=np.zeros(arrIn.shape)
    arrCount[np.where(~np.isnan(arrIn))]=1
    ArrRMSD=np.nansum(arrIn,axis=0)
    ArrCount=np.nansum(arrCount,axis=0)
    ArrRMSD=np.divide(ArrRMSD,ArrCount)    
    ArrRMSD=np.sqrt(ArrRMSD)

    #print("RMSD array shape:", ArrRMSD.shape)
    #arrmean=np.reshape(arrmean,(T,J*I))
    #arrmean=np.nanmean(arrmean,axis=1)
    return ArrRMSD, ArrCount

def combineRMSD(arrin,CountIn,nansplj,nanspli):
    RMSDin=np.copy(arrin)
    RMSDin[::,nansplj,nanspli]=np.nan
    RMSDin=np.square(RMSDin)
    RMSDin=np.multiply(RMSDin,CountIn)
    hourRMSD=[]
    hourCount=[]
    for i in range(24):
        hourRMSD.append(np.nansum(RMSDin[i::24]))
        hourCount.append(np.nansum(CountIn[i::24]))
    hourRMSD=np.array(hourRMSD)
    hourCount=np.array(hourCount)
    RMSDout=np.divide(hourRMSD,hourCount)
    RMSDout=np.sqrt(RMSDout)
    return RMSDout,hourCount

#start_hour=0
#end_hour=0

for month in range(start_month,end_month+1):
    m='{:02d}'.format(month)
    if mesotype=='obs':
        print(mesonet_path,m,fext)
        files=sorted(glob.glob('{}/OKMesonet.????{}.{}'.format(mesonet_path,m,fext)))
        print(files)
    elif mesotype=='WBGT':
        files=sorted(glob.glob('{}/WBGT-uncorrected-SRADnew/WBGT.????{}.{}'.format(mesonet_path,m,fext)))
    #print(files)
    if month in [1,3,5,7,8,10,12]:
        M_days=31
    elif month in [4,6,9,11]:
        M_days=30
    elif month==2:
        M_days=29

    bWBGT=np.ones((Num_years,M_days*Num_hours,Num_lat,Num_lon),dtype=np.float32)*np.nan
    if mesotype=='obs':
        bTDEW=np.ones((Num_years,M_days*Num_hours,Num_lat,Num_lon),dtype=np.float32)*np.nan
        bTAIR=np.ones((Num_years,M_days*Num_hours,Num_lat,Num_lon),dtype=np.float32)*np.nan
        bSRAD=np.ones((Num_years,M_days*Num_hours,Num_lat,Num_lon),dtype=np.float32)*np.nan
        bWS2M=np.ones((Num_years,M_days*Num_hours,Num_lat,Num_lon),dtype=np.float32)*np.nan
        bPRES=np.ones((Num_years,M_days*Num_hours,Num_lat,Num_lon),dtype=np.float32)*np.nan
        bWSPD=np.ones((Num_years,M_days*Num_hours,Num_lat,Num_lon),dtype=np.float32)*np.nan

#    eWBGT=np.ones((M_days*Num_hours,Num_lat,Num_lon))*np.nan
#    eTDEW=np.ones((M_days*Num_hours,Num_lat,Num_lon))*np.nan
#    eTAIR=np.ones((M_days*Num_hours,Num_lat,Num_lon))*np.nan
#    eSRAD=np.ones((M_days*Num_hours,Num_lat,Num_lon))*np.nan
#    eWS2M=np.ones((M_days*Num_hours,Num_lat,Num_lon))*np.nan
#    ePRES=np.ones((M_days*Num_hours,Num_lat,Num_lon))*np.nan
#    eWSPD=np.ones((M_days*Num_hours,Num_lat,Num_lon))*np.nan
#Create 4 D array (year, hour, space (x2)). I need at least 1 for the RMSD and 1 for the number of data points
    #start_hour=end_hour
    #end_hour=end_hour+M_days*Num_hours

###############
#Read in files#
###############
    for y, f in enumerate(files): #Need to adjust file paths/names
        print(y+start_year)
        eWBGT=np.ones((M_days*Num_hours,Num_lat,Num_lon),dtype=np.float32)*np.nan
        if mesotype=='WBGT':
            mWBGT=np.ones((M_days*Num_hours,361,720),dtype=np.float32)*np.nan
        elif mesotype=='obs':
            eTDEW=np.ones((M_days*Num_hours,Num_lat,Num_lon),dtype=np.float32)*np.nan
            eTAIR=np.ones((M_days*Num_hours,Num_lat,Num_lon),dtype=np.float32)*np.nan
            eSRAD=np.ones((M_days*Num_hours,Num_lat,Num_lon),dtype=np.float32)*np.nan
            eWS2M=np.ones((M_days*Num_hours,Num_lat,Num_lon),dtype=np.float32)*np.nan
            ePRES=np.ones((M_days*Num_hours,Num_lat,Num_lon),dtype=np.float32)*np.nan
            eWSPD=np.ones((M_days*Num_hours,Num_lat,Num_lon),dtype=np.float32)*np.nan


        print(f)
        #Read in mesonet data
        mesof=f#'{}/OKMesonet.????{}.nc'.format(mesonet_path,m)

        Varn, Varnm, mlat, mlon = ReadNetCDF(mesof,mesonet_vars)
        if mesotype=='obs':
            mlat=mlat[:][::-1] #Check these for mesotype==WBGT
            mlon=mlon[:]+360 #Check these for mesotype==WBGT
        elif mesotype=='WBGT':
            mlat=mlat[:]
            mlon=mlon[:]

        #print('mlat,mlon',mlat,mlon)
        time=Varn[0]

        if mesotype=='WBGT':
            month_days=Varn[1].shape[0]
            mWBGT[:month_days]=Varn[1][::,::,::].astype(np.float32)
        elif mesotype=='obs':
            mWBGT=Varn[1][::,::-1,::].astype(np.float32)
            mRELH=Varn[2][::,::-1,::].astype(np.float32)
            mTAIR=Varn[3][::,::-1,::].astype(np.float32)
            mSRAD=Varn[4][::,::-1,::].astype(np.float32)
            mWS2M=Varn[5][::,::-1,::].astype(np.float32)
            mPRES=Varn[6][::,::-1,::].astype(np.float32)
            mWSPD=Varn[7][::,::-1,::].astype(np.float32)
            mTDEW=Vapr_to_T(Vapr(mTAIR)*mRELH*0.01).astype(np.float32)
            mTDEW[np.where(mTDEW.mask==True)]=np.nan
#        print(mRELH)
#        print(mTAIR)
#        print(mTDEW.data)
        
        if mesotype=='obs':
            maxlat=np.nanmax(mlat)
            minlat=np.nanmin(mlat)
            maxlon=np.nanmax(mlon)
            minlon=np.nanmin(mlon)
        elif mesotype=='WBGT':
            maxlat=38.0
            minlat=33.0
            maxlon=-93.5+360 #Check format
            minlon=-104.0+360 #Check format

            #Trim to desired grid
            latinds=np.where((mlat<=maxlat)&(mlat>=minlat))[0]
            loninds=np.where((mlon<=maxlon)&(mlon>=minlon))[0]
            #print(latinds,loninds)
            mWBGT_temp=mWBGT[:,latinds]
            mWBGT_temp=mWBGT_temp[::,::,loninds]
            #print(mWBGT_temp.shape)
            mWBGT=mWBGT_temp

            mlat=mlat[latinds]
            mlon=mlon[loninds]
        
        #Read in ERA5 data
        EVarn=[]
        for i, var in enumerate(mesonet_vars[1:]):
            EVar=ERA5_vars[i+1]
            if EVar=='WBGT':
                ERA5f='{}/{}-Liljegren/{}.{}{}.{}'.format(ERA5_path,EVar,EVar,str(start_year+y),m,fext)
            elif mesotype=='obs':
                ERA5f='{}/{}/{}.{}{}.{}'.format(ERA5_path,EVar,EVar,str(start_year+y),m,fext)
            Varn, Varnm, lat, lon = ReadNetCDF(ERA5f,['time',EVar])
            lat=lat[:]
            lon=lon[:]
            #print('lat,lon',lat,lon)
            #print(lon)
            #lon[np.where(lon>180.0)]=lon[np.where(lon>180.0)]-360.0
            #print(lon)

            #Trim to mesonet grid
            latinds=np.where((lat<=maxlat)&(lat>=minlat))[0]
            loninds=np.where((lon<=maxlon)&(lon>=minlon))[0]
            print(latinds,loninds)
            EVarn_temp=Varn[1][:,latinds]
            EVarn_temp=EVarn_temp[::,::,loninds]
            print(EVarn_temp.shape)
            EVarn.append(EVarn_temp)

        lat=lat[latinds]
        lon=lon[loninds]
        eWBGT[:EVarn[0].shape[0],::,::]=EVarn[0].astype(np.float32)
        if mesotype=='WBGT':
            eWBGT=Convert_T_Unit(eWBGT,'K','F')
        if mesotype=='obs':
            eTDEW[:EVarn[1].shape[0],::,::]=EVarn[1].astype(np.float32)
            #if Stat_type=='Mean':
            eTDEW=eTDEW-273.15
            print("************Dewpoint TEST*************")
            print("min eTDEW:",np.nanmin(eTDEW))
            print("max eTDEW:",np.nanmax(eTDEW))
            eTAIR[:EVarn[2].shape[0],::,::]=EVarn[2].astype(np.float32)
            eTAIR=eTAIR-273.15
            print("************Temp TEST*************")
            print("min eTAIR:",np.nanmin(eTAIR))
            print("max eTAIR:",np.nanmax(eTAIR))
            eSRAD[:EVarn[3].shape[0],::,::]=EVarn[3].astype(np.float32)
            eWS2M[:EVarn[4].shape[0],::,::]=EVarn[4].astype(np.float32)#*1.25
            ePRES[:EVarn[5].shape[0],::,::]=EVarn[5].astype(np.float32)/100
            eWSPD[:EVarn[6].shape[0],::,::]=EVarn[6].astype(np.float32)

        mlon,mlat=np.meshgrid(mlon,mlat)
        elon,elat=lon,lat
        lon,lat=np.meshgrid(lon,lat)

        print(eWBGT.shape,mWBGT.shape)
        bWBGT[y,::,::,::]=eWBGT-mWBGT
        if mesotype=='obs':
            bTDEW[y,::,::]=eTDEW-mTDEW
            bTAIR[y,::,::]=eTAIR-mTAIR
            bSRAD[y,::,::]=eSRAD-mSRAD
            bWS2M[y,::,::]=eWS2M-mWS2M
            bPRES[y,::,::]=ePRES-mPRES
            bWSPD[y,::,::]=eWSPD-mWSPD

#################
#Calculate Stats#
#################
#Trim to plains domain
    if mesotype=='obs':
        nanspl=np.where(lon>360.0-97.0)
    elif mesotype=='WBGT':
        nanspl=np.where(lon>maxlon)
    nansplj=nanspl[0]
    nanspli=nanspl[1]
#Calculate RMSD for each grid point

#Calculate number of points at each grid point

    rWBGT,cWBGT=gridRMSD(bWBGT)
    if mesotype=='obs':
        rTDEW,cTDEW=gridRMSD(bTDEW)
        rTAIR,cTAIR=gridRMSD(bTAIR)
        rSRAD,cSRAD=gridRMSD(bSRAD)
        rWS2M,cWS2M=gridRMSD(bWS2M)
        rPRES,cPRES=gridRMSD(bPRES)
        rWSPD,cWSPD=gridRMSD(bWSPD)

#Sum up for each day
    hrWBGT,hcWBGT=combineRMSD(rWBGT,cWBGT,nansplj,nanspli)
    if mesotype=='obs':
        hrTDEW,hcTDEW=combineRMSD(rTDEW,cTDEW,nansplj,nanspli)
        hrTAIR,hcTAIR=combineRMSD(rTAIR,cTAIR,nansplj,nanspli)
        hrSRAD,hcSRAD=combineRMSD(rSRAD,cSRAD,nansplj,nanspli)
        hrWS2M,hcWS2M=combineRMSD(rWS2M,cWS2M,nansplj,nanspli)
        hrPRES,hcPRES=combineRMSD(rPRES,cPRES,nansplj,nanspli)
        hrWSPD,hcWSPD=combineRMSD(rWSPD,cWSPD,nansplj,nanspli)
       
#################
# Write to file #
#################

    print('Printing stats')
    print('t,lat,lon',time.shape,elat.shape,elon.shape)
    print(rWBGT.shape,cWBGT.shape,hrWBGT.shape,hcWBGT.shape)

    Dims=[time,hours,elat,elon]
    DimsName=['time','hours','latitude','longitude']
    DimsUnits=['hours since 1900-01-01 00:00:00.0','hour UTC','degrees_north','degrees_east']
    DimsLN=DimsName#['time','latitude','longitude']
    if mesotype=='obs':
        Vars=[rWBGT,cWBGT,hrWBGT,hcWBGT,rTDEW,cTDEW,hrTDEW,hcTDEW,rTAIR,cTAIR,hrTAIR,hcTAIR,rSRAD,cSRAD,hrSRAD,hcSRAD,rWS2M,cWS2M,hrWS2M,hcWS2M,rPRES,cPRES,hrPRES,hcPRES,rWSPD,cWSPD,hrWSPD,hcWSPD]
        VarsNames=['rWBGT','cWBGT','hrWBGT','hcWBGT',\
               'rTDEW','cTDEW','hrTDEW','hcTDEW',\
               'rTAIR','cTAIR','hrTAIR','hcTAIR',\
               'rSRAD','cSRAD','hrSRAD','hcSRAD',\
               'rWS2M','cWS2M','hrWS2M','hcWS2M',\
               'rPRES','cPRES','hrPRES','hcPRES',\
               'rWSPD','cWSPD','hrWSPD','hcWSPD']
        VarsUnits=['F','F','F','F','K','K','K','K','K','K','K','K',"W m**-2","W m**-2","W m**-2","W m**-2","m s**-1","m s**-1","m s**-1","m s**-1",'Pa','Pa','Pa','Pa',"m s**-1","m s**-1","m s**-1","m s**-1"]
        VarsLN=['WBGT RMSD grid','WBGT count grid','WBGT RMSD hourly','WBGT count hourly',\
               'TDEW RMSD grid','TDEW count grid','TDEW RMSD hourly','TDEW count hourly',\
               'TAIR RMSD grid','TAIR count grid','TAIR RMSD hourly','TAIR count hourly',\
               'SRAD RMSD grid','SRAD count grid','SRAD RMSD hourly','SRAD count hourly',\
               'WS2M RMSD grid','WS2M count grid','WS2M RMSD hourly','WS2M count hourly',\
               'PRES RMSD grid','PRES count grid','PRES RMSD hourly','PRES count hourly',\
               'WSPD RMSD grid','WSPD count grid','WSPD RMSD hourly','WSPD count hourly']
        VarsDims=[('time','latitude','longitude'),('time','latitude','longitude'),('hours',),('hours',),('time','latitude','longitude'),('time','latitude','longitude'),('hours',),('hours',),('time','latitude','longitude'),('time','latitude','longitude'),('hours',),('hours',),('time','latitude','longitude'),('time','latitude','longitude'),('hours',),('hours',),('time','latitude','longitude'),('time','latitude','longitude'),('hours',),('hours',),('time','latitude','longitude'),('time','latitude','longitude'),('hours',),('hours',),('time','latitude','longitude'),('time','latitude','longitude'),('hours',),('hours',)]
    else:
        Vars=[rWBGT,cWBGT,hrWBGT,hcWBGT]
        VarsNames=['rWBGT','cWBGT','hrWBGT','hcWBGT']
        VarsUnits=['F','F','F','F']
        VarsLN=['WBGT RMSD grid','WBGT count grid','WBGT RMSD hourly','WBGT count hourly']
        VarsDims=[('time','latitude','longitude'),('time','latitude','longitude'),('hours',),('hours',)]

    output_dir='{}/{}'.format(Output_path,'RMSD')
    filename='{}-uncorrected-SRADinterp.{}.{}'.format('RMSDinfo',m,fext)

    WriteNetCDF(Dims,DimsName,DimsUnits,DimsLN,Vars,VarsNames,VarsUnits,VarsLN,VarsDims,output_dir,filename)


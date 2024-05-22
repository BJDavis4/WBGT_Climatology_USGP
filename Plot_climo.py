from netCDF4 import Dataset
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
from netCDF_mods import ReadNetCDF, WriteNetCDF
from Python_mods import Calc_ws, interp_ws_log, Convert_T_Unit

mesonet_path='/data/deluge/observations/OKMesonet/grid'
ERA5_path='/data/deluge/scratch/ERA5-Ben/2D/hourly'
mesotype='obs'#['WBGT','obs]

mesonet_vars=['time','TWBG','TDEW','TAIR','SRAD','WS2M','PRES','WSPD']
ERA5_vars=['time','WBGT','d2m','t2m','msdwswrf','ws2m','sp','ws10']
var_units_plt=['seconds',r'${\degree}F$',r'${\degree}C$',r'${\degree}C$',r'$W m^\mathregular{-2}$',r'$m s^\mathregular{-1}$',r'$hPa$',r'$m s^\mathregular{-1}$'] #mesonet, used for plotting

if mesotype=='WBGT':
    mesonet_path=ERA5_path
    mesonet_vars=ERA5_vars[:2]
    ERA5_vars=ERA5_vars[:2]
    var_units_plt=var_units_plt[:2]

#W m**-2","W m**-2","W m**-2","m s**-1","m s**-1","m s**-1",'Pa','Pa','Pa',"m s**-1","m s**-1","m s**-1"]


Output_path=ERA5_path
Stat_type='STDev' #Mean #STDev
Calc_bias=True
Plot_data=False
save_path='/home/bjdavis/python/figures/Climo-Liljegren'


#vars_dict={'TWBG':[eWBGT,mWBGT],'TDEW':[eTDEW,mSRAD],'TAIR':[eTAIR,mTAIR],'SRAD':[eSRAD,mSRAD],'WS2M':[eWS2M,mWS2M],'PRES':[ePRES,mPRES],'WSPD':[eWSPD,mWSPD]}

def gridmean(arrin,nans,nansplij,nanspli,T,J,I):
    arrmean=np.copy(arrin)
    arrmean[nans]=np.nan
    arrmean[::,nansplj,nanspli]=np.nan
    arrmean=np.reshape(arrmean,(T,J*I))
    arrmean=np.nanmean(arrmean,axis=1)
    return arrmean

if Stat_type=='Mean':
    fext='ltm'
elif Stat_type=="STDev":
    fext='ltstd'

for month in range(7,8):
    m='{:02d}'.format(month)
    gpmaxlat=50.0
    if mesotype=='obs':
        mesof='{}/OKMesonet.hourly.{}.{}'.format(mesonet_path,m,fext)
    elif mesotype=='WBGT':
        mesof='{}/WBGT-uncorrected-SRADnew/WBGT.hourly.{}.{}'.format(mesonet_path,m,fext)

    Varn, Varnm, mlat, mlon = ReadNetCDF(mesof,mesonet_vars)
    if mesotype=='obs':
        mlat=mlat[:][::-1]
        mlon=mlon[:]+360
    elif mesotype=='WBGT':
        mlat=mlat[:]
        mlon=mlon[:]
    time=Varn[0]
    if mesotype=='WBGT':
        if month==2:
            mWBGT=np.ones((29*24,361,720))*np.nan
            mWBGT[:Varn[1].shape[0],::,::]=Varn[1][::,::,::].astype(np.float32)
        else:
            mWBGT=Varn[1][::,::,::].astype(np.float32)
        maxlat=38.0
        minlat=33.0
        maxlon=-93.5+360
        minlon=-104.0+360

        latinds=np.where((mlat<=maxlat)&(mlat>=minlat))[0]
        loninds=np.where((mlon<=maxlon)&(mlon>=minlon))[0]
        gplatinds=np.where((mlat<=gpmaxlat)&(mlat>=minlat))[0]
        print(latinds,loninds)
        mWBGT_temp=mWBGT[:,latinds].copy().astype(np.float32)
        mWBGT_temp=mWBGT_temp[::,::,loninds].astype(np.float32)
        mWBGTplot_temp=mWBGT[:,gplatinds].astype(np.float32)
        mWBGTplot_temp=mWBGTplot_temp[::,::,loninds].astype(np.float32)
        mWBGT=mWBGT_temp
        mWBGTplot=mWBGTplot_temp

    elif mesotype=='obs':
        mWBGT=Varn[1][::,::-1,::].astype(np.float32)
        mTDEW=Varn[2][::,::-1,::].astype(np.float32)
        mTAIR=Varn[3][::,::-1,::].astype(np.float32)
        mSRAD=Varn[4][::,::-1,::].astype(np.float32)
        mWS2M=Varn[5][::,::-1,::].astype(np.float32)
        mPRES=Varn[6][::,::-1,::].astype(np.float32)
        mWSPD=Varn[7][::,::-1,::].astype(np.float32)

        maxlat=np.nanmax(mlat)
        minlat=np.nanmin(mlat)
        maxlon=np.nanmax(mlon)
        minlon=np.nanmin(mlon)

    print(maxlat,minlat,maxlon,minlon)    
    #exit()

    EVarn=[]
    EVarnplot=[]
    for i, var in enumerate(mesonet_vars[1:]):
        #if var=='WSPD':
        #    EVar=ERA5_vars[i+1]
        #    ERA5f='{}/{}/{}.hourly.{}.{}'.format(ERA5_path,EVar,EVar,m,fext)
        #    Varn, Varnm, lat, lon = ReadNetCDF(ERA5f,['time',EVar])
        #    u10=Varn[1]
        #    EVar=ERA5_vars[i+2]
        #    ERA5f='{}/{}/{}.hourly.{}.{}'.format(ERA5_path,EVar,EVar,m,fext)
        #    Varn, Varnm, lat, lon = ReadNetCDF(ERA5f,['time',EVar])
        #    v10=Varn[1]
        #    Varn[1]=Calc_ws(u10,v10)#np.sqrt((u10**2)+(v10**2))

        #else:
        EVar=ERA5_vars[i+1]
        if mesotype=='WBGT':
            ERA5f='{}/{}-Liljegren/{}.hourly.{}.{}'.format(ERA5_path,EVar,EVar,m,fext)
        elif mesotype=='obs':
            ERA5f='{}/{}/{}.hourly.{}.{}'.format(ERA5_path,EVar,EVar,m,fext)
        Varn, Varnm, lat, lon = ReadNetCDF(ERA5f,['time',EVar])
            #if var=='WS2M':
            #    Varn[1]=Varn[1][::,::-1,::]
            #    lat=lat[::-1]
        lat=lat[:]
        lon=lon[:]
        print(var)
        print(lat)
        #lon[np.where(lon>180.0)]=lon[np.where(lon>180.0)]-360.0
        #print(lon)
        latinds=np.where((lat<=maxlat)&(lat>=minlat))[0]
        loninds=np.where((lon<=maxlon)&(lon>=minlon))[0]
        gplatinds=np.where((lat<=gpmaxlat)&(lat>=minlat))[0]
        print(latinds,loninds)
        EVarn_temp=Varn[1][:,latinds].astype(np.float32)
        EVarn_temp=EVarn_temp[::,::,loninds].astype(np.float32)
        EVarnplot_temp=Varn[1][:,gplatinds].astype(np.float32)
        EVarnplot_temp=EVarnplot_temp[::,::,loninds].astype(np.float32)
        print(EVarn_temp.shape)
        EVarn.append(EVarn_temp)
        EVarnplot.append(EVarnplot_temp)

    #u10=u10[:,latinds]
    #u10=u10[::,::,loninds]
    #v10=v10[:,latinds]
    #v10=v10[::,::,loninds]

    if mesotype=='WBGT':
        if Stat_type=='Mean':
            eWBGT=Convert_T_Unit(EVarn[0],'K','F')
        elif Stat_type=='STDev':
            eWBGT=EVarn[0]*1.8
    elif mesotype=='obs':
        eWBGT=EVarn[0]
        eTDEW=EVarn[1]
        if Stat_type=='Mean':
            eTDEW=eTDEW-273.15
        print("************DEWPOINT TEST*************")
        print(np.nanmin(eTDEW))
        print(np.nanmax(eTDEW))
        eTAIR=EVarn[2]
        if Stat_type=='Mean':
            eTAIR=eTAIR-273.15
        print("************Temp TEST*************")
        print(np.nanmin(eTAIR))
        print(np.nanmax(eTAIR))
        eSRAD=EVarn[3]
        eWS2M=EVarn[4]
        ePRES=EVarn[5]/100
        eWSPD=EVarn[6]

    if mesotype=='WBGT':
        if Stat_type=='Mean':
            eWBGTp=Convert_T_Unit(EVarnplot[0],'K','F')
        elif Stat_type=='STDev':
            eWBGTp=EVarnplot[0]*1.8
    elif mesotype=='obs':
        eWBGTp=EVarnplot[0]
        eTDEWp=EVarnplot[1]
        if Stat_type=='Mean':
            eTDEWp=eTDEWp-273.15
        print("************DEWPOINT TEST*************")
        print(np.nanmin(eTDEWp))
        print(np.nanmax(eTDEWp))
        eTAIRp=EVarnplot[2]
        if Stat_type=='Mean':
            eTAIRp=eTAIRp-273.15
        print("************Temp TEST*************")
        print(np.nanmin(eTAIRp))
        print(np.nanmax(eTAIRp))
        eSRADp=EVarnplot[3]
        eWS2Mp=EVarnplot[4]
        ePRESp=EVarnplot[5]/100
        eWSPDp=EVarnplot[6]

#    print('u10:', u10[0][52][72])
#    print('v10:', v10[0][52][72])
#    print('ws10:', eWSPD[0][52][72])
#    print('ws10v2:', Calc_ws(u10[0][52][72],v10[0][52][72]))
#    print('ws2m:', eWS2M[0][52][72])
#    print('ws2mv2:', interp_ws_log(2.0,np.ones(eWSPD.shape)*0.307,10.0,eWSPD[0][52][72]))
#    exit()

    if mesotype=='WBGT':
        latinds=np.where((lat<=maxlat)&(lat>=minlat))[0]
        loninds=np.where((lon<=maxlon)&(lon>=minlon))[0]
        gplatinds=np.where((lat<=gpmaxlat)&(lat>=minlat))[0]
        print(latinds,loninds)
        EVarn_temp=Varn[1][:,latinds].astype(np.float32)
        EVarn_temp=EVarn_temp[::,::,loninds].astype(np.float32)
        EVarnplot_temp=Varn[1][:,gplatinds].astype(np.float32)
        EVarnplot_temp=EVarnplot_temp[::,::,loninds].astype(np.float32)

    if mesotype=='obs':
        mlon,mlat=np.meshgrid(mlon,mlat)
    elif mesotype=='WBGT':
        mlon,mlat=np.meshgrid(mlon[loninds],mlat[latinds])
    lonp,latp=np.meshgrid(lon[loninds],lat[gplatinds])
    lon,lat=np.meshgrid(lon[loninds],lat[latinds])

    dWBGT=eWBGT-mWBGT#[::,::-1,::]
    if mesotype=='obs':
        dTDEW=eTDEW-mTDEW#[::,::-1,::]
        dTAIR=eTAIR-mTAIR#[::,::-1,::]
        dSRAD=eSRAD-mSRAD#[::,::-1,::]
        dWS2M=eWS2M-mWS2M#[::,::-1,::]
        dPRES=ePRES-mPRES#[::,::-1,::]
        dWSPD=eWSPD-mWSPD#[::,::-1,::]

    print('check lons')
    print(mlon)
    print(lon)
    print('check lats')
    print(mlat)
    print(lat)

    #############################
    ### Calculate single bias ### probably put into an if statement so I can reuse rest of code for stdev
    #############################
    if Calc_bias:
        nans=np.where(np.isnan(dWBGT))
        if mesotype=='obs':
            nanspl=np.where(lon>360.0-97.0)
        else:
            nanspl=np.where(lon>360.0-90.0)
        nansplj=nanspl[0]
        nanspli=nanspl[1]
        print(nanspl)
        #print(nanspl.shape)
        T,J,I=dWBGT.shape
    
        eWBGTm=gridmean(eWBGT,nans,nansplj,nanspli,T,J,I)
        mWBGTm=gridmean(mWBGT,nans,nansplj,nanspli,T,J,I)
        print(eWBGTm-mWBGTm)
        #print(eWBGTm2)
    
        bWBGT=gridmean(dWBGT,nans,nansplj,nanspli,T,J,I)
        print('WBGT bias shape:', bWBGT.shape)
        print(eWBGTm)
        print(mWBGTm)
        print(bWBGT)
        #print(bWBGT2)
        #exit()
    
        if mesotype=='obs':
            eTDEWm=gridmean(eTDEW,nans,nansplj,nanspli,T,J,I)
            mTDEWm=gridmean(mTDEW,nans,nansplj,nanspli,T,J,I)
            bTDEW=gridmean(dTDEW,nans,nansplj,nanspli,T,J,I)
        
            eTAIRm=gridmean(eTAIR,nans,nansplj,nanspli,T,J,I)
            mTAIRm=gridmean(mTAIR,nans,nansplj,nanspli,T,J,I)
            bTAIR=gridmean(dTAIR,nans,nansplj,nanspli,T,J,I)
        
            eSRADm=gridmean(eSRAD,nans,nansplj,nanspli,T,J,I)
            mSRADm=gridmean(mSRAD,nans,nansplj,nanspli,T,J,I)
            bSRAD=gridmean(dSRAD,nans,nansplj,nanspli,T,J,I)
        
            eWS2Mm=gridmean(eWS2M,nans,nansplj,nanspli,T,J,I)
            mWS2Mm=gridmean(mWS2M,nans,nansplj,nanspli,T,J,I)
            bWS2M=gridmean(dWS2M,nans,nansplj,nanspli,T,J,I)
        
            ePRESm=gridmean(ePRES,nans,nansplj,nanspli,T,J,I)*100
            mPRESm=gridmean(mPRES,nans,nansplj,nanspli,T,J,I)*100
            bPRES=gridmean(dPRES,nans,nansplj,nanspli,T,J,I)*100
        
            eWSPDm=gridmean(eWSPD,nans,nansplj,nanspli,T,J,I)
            mWSPDm=gridmean(mWSPD,nans,nansplj,nanspli,T,J,I)
            bWSPD=gridmean(dWSPD,nans,nansplj,nanspli,T,J,I)
    
        ###########################
        ### Write stuff to file ###
        ###########################
            
        Dims=[time]
        DimsName=['time']
        DimsUnits=['hours since 1900-01-01 00:00:00.0']
        DimsLN=['time']
        if mesotype=='obs':
            Vars=[eWBGTm,mWBGTm,bWBGT,eTDEWm,mTDEWm,bTDEW,eTAIRm,mTAIRm,bTAIR,eSRADm,mSRADm,bSRAD,eWS2Mm,mWS2Mm,bWS2M,ePRESm,mPRESm,bPRES,eWSPDm,mWSPDm,bWSPD]
            VarsNames=['eWBGTmean','mWBGTmean','WBGTbias',\
                       'eTDEWmean','mTDEWmean','TDEWbias',\
                       'eTAIRmean','mTAIRmean','TAIRbias',\
                       'eSRADmean','mSRADmean','SRADbias',\
                       'eWS2Mmean','mWS2Mmean','WS2Mbias',\
                       'ePRESmean','mPRESmean','PRESbias',\
                       'eWSPDmean','mWSPDmean','WSPDbias']
            VarsUnits=['F','F','F','K','K','K','K','K','K',"W m**-2","W m**-2","W m**-2","m s**-1","m s**-1","m s**-1",'Pa','Pa','Pa',"m s**-1","m s**-1","m s**-1"]
            VarsLN=['ERA5 WBGT mean','OKMesonet WBGT mean','ERA5 WBGT bias',\
                       'ERA5 TDEW mean','OKMesonet TDEW mean','ERA5 TDEW bias',\
                       'ERA5 TAIR mean','OKMesonet TAIR mean','ERA5 TAIR bias',\
                       'ERA5 SRAD mean','OKMesonet SRAD mean','ERA5 SRAD bias',\
                       'ERA5 WS2M mean','OKMesonet WS2M mean','ERA5 WS2M bias',\
                       'ERA5 PRES mean','OKMesonet PRES mean','ERA5 PRES bias',\
                       'ERA5 WSPD mean','OKMesonet WSPD mean','ERA5 WSPD bias']
            VarsDims=[('time',),('time',),('time',),('time',),('time',),('time',),('time',),('time',),('time',),('time',),('time',),('time',),('time',),('time',),('time',),('time',),('time',),('time',),('time',),('time',),('time',)]
        if mesotype=='WBGT':
            Vars=[eWBGTm,mWBGTm,bWBGT]
            VarsNames=['eWBGTmean','mWBGTmean','WBGTbias']
            VarsUnits=['F','F','F']    
            VarsLN=['ERA5 WBGT mean','OKMesonet WBGT mean','ERA5 WBGT bias']
            VarsDims=[('time',),('time',),('time',)]

        output_dir='{}/{}'.format(Output_path,'bias')
        if mesotype=='WBGT':
            filename='{}-Liljegren.{}.{}.nc'.format('biasinfo',m,fext)
        elif mesotype=='obs':
            filename='{}-uncorrected-SRADinterp.{}.{}.nc'.format('biasinfo',m,fext)
    
        WriteNetCDF(Dims,DimsName,DimsUnits,DimsLN,Vars,VarsNames,VarsUnits,VarsLN,VarsDims,output_dir,filename)
    
    
    ##################
    ### Plot Stuff ###
    ##################
    if Plot_data:
        vars_dict={'TWBG':[eWBGTp,mWBGT,dWBGT],'TDEW':[eTDEWp,mTDEW,dTDEW],'TAIR':[eTAIRp,mTAIR,dTAIR],'SRAD':[eSRADp,mSRAD,dSRAD],'WS2M':[eWS2Mp,mWS2M,dWS2M],'PRES':[ePRESp,mPRES,dPRES],'WSPD':[eWSPDp,mWSPD,dWSPD]}

        if Stat_type=='Mean':
            lev_start=[20, -10,-10, 0,   0,  880, 0]
            lev_end=  [101,27, 37,1001,7.1,1011,10.1]
            lev_inc=  [5,  2,  1, 50,  0.5,10,  0.5]
    
            dif_lev_start=[-15,-10,-10,-100,-5,-10,-5]
            dif_lev_end=[16,11,11,101,5.1,11,5.1]
            dif_lev_inc=[1,1,1,10,0.5,1,0.5]
    
        elif Stat_type=='STDev':
            lev_start=[0,   0,   0,   0,  0,  0, 0]
            lev_end=  [15.1,15.1,15.1,351,5.1,15.1,10.1]
            lev_inc=  [5,   0.5, 0.5, 25, 0.5,1,  0.5]

            dif_lev_start=[-15,-5,-5,  -100,-5, -10,-5]
            dif_lev_end=  [16, 5.1,5.1,101, 5.1,11,5.1]
            dif_lev_inc=  [1,  0.5,0.5,10,  0.5,1,0.5]
        mo=m
        for h in range(0,24):
            hh='{:02d}'.format(h)
            for i, var in enumerate(mesonet_vars[1:]):
                pltsave='{}/{}_{}_{}{}.png'.format(save_path,var,Stat_type,mo,hh)
                print(var)
                eVar=vars_dict[var][0][h]
                mVar=vars_dict[var][1][h]
                dVar=vars_dict[var][2][h]
                print('OK Mesonet')
                print('data min:',np.nanmin(mVar))
                print('data max:',np.nanmax(mVar))
                print('ERA5')
                print('data min:',np.nanmin(eVar))
                print('data max:',np.nanmax(eVar))
                
                fig = plt.figure(figsize=(9,5))
#                gs = gridspec.GridSpec(3,7)
                gs = gridspec.GridSpec(2,2)
#                ax=plt.subplot(gs[0,2:5])
                ax=plt.subplot(gs[:,0])
                m=Basemap(resolution='i', projection='merc', llcrnrlat=minlat, urcrnrlat=gpmaxlat,llcrnrlon=minlon, urcrnrlon=maxlon, lat_ts=35)
                m.drawcountries(linewidth=0.25)
                m.drawstates(linewidth=0.15)
                x,y=m(lonp,latp)
                cmap=plt.get_cmap('RdYlGn_r')
                m.contourf(x,y,eVar,cmap=cmap,extend='both',levels=np.arange(lev_start[i],lev_end[i],lev_inc[i]))
                plt.colorbar(label=var_units_plt[i+1])
                plt.title('ERA5')
    
#                ax=plt.subplot(gs[1,2:5])
                plt.subplot(gs[0,1])
                m=Basemap(resolution='i', projection='merc', llcrnrlat=minlat, urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon, lat_ts=35)
                m.drawcountries(linewidth=0.25)
                m.drawstates(linewidth=0.15)
                x,y=m(mlon,mlat)
                cmap=plt.get_cmap('RdYlGn_r')
                m.contourf(x,y,mVar,cmap=cmap,extend='both',levels=np.arange(lev_start[i],lev_end[i],lev_inc[i]))
                plt.colorbar(label=var_units_plt[i+1])
                plt.title('OK Mesonet')
    
#                ax=plt.subplot(gs[2,2:5])
                ax=plt.subplot(gs[1,1])
                m=Basemap(resolution='i', projection='merc', llcrnrlat=minlat, urcrnrlat=maxlat,llcrnrlon=minlon, urcrnrlon=maxlon, lat_ts=35)
                m.drawcountries(linewidth=0.25)
                m.drawstates(linewidth=0.15)
                x,y=m(lon,lat)
                cmap=plt.get_cmap('RdYlGn_r')
                m.contourf(x,y,dVar,cmap=cmap,extend='both',levels=np.arange(dif_lev_start[i],dif_lev_end[i],dif_lev_inc[i]))
                plt.colorbar(label=var_units_plt[i+1])
                plt.title('ERA5 minus OK Mesonet')
    
                plt.suptitle('{} for month {} at {} UTC'.format(var,mo,hh))
                plt.tight_layout()

                plt.savefig(pltsave)
                #plt.show()
                plt.close()

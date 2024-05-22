import numpy as np
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from netCDF_mods import ReadNetCDF

##################
# Set parameters #
##################
start_month=7#1
end_month=7#12
Stat_type='Mean'#'STDev'
rb='both' #rmsd bias both
mesotype='WBGT' #obs WBGT
if mesotype=='WBGT':
    formula='Liljegren'
elif mesotype=='obs':
    formula="Dimiceli" #Dimiceli Liljegren
inpath_bias='/data/deluge/scratch/ERA5-Ben/2D/hourly/bias'
inpath_rmsd='/data/deluge/scratch/ERA5-Ben/2D/hourly/RMSD'
Vars=['eWBGTmean','mWBGTmean','WBGTbias',\
                   'eTDEWmean','mTDEWmean','TDEWbias',\
                   'eTAIRmean','mTAIRmean','TAIRbias',\
                   'eSRADmean','mSRADmean','SRADbias',\
                   'eWS2Mmean','mWS2Mmean','WS2Mbias',\
                   'ePRESmean','mPRESmean','PRESbias',\
                   'eWSPDmean','mWSPDmean','WSPDbias']

RVars=['hrWBGT','hrTDEW','hrTAIR','hrSRAD','hrWS2M','hrPRES','hrWSPD']
var_units_plt=[r'${\degree}F$',r'${\degree}C$',r'${\degree}C$',r'$W m^\mathregular{-2}$',r'$m s^\mathregular{-1}$',r'$hPa$',r'$m s^\mathregular{-1}$']
#Vars=['eWS2Mmean','mWS2Mmean','WS2Mbias']
figdpi=1000#1000

if mesotype=='WBGT':
    Vars=Vars[:3]
    RVars=RVars[:1]
    var_units_plt[:1]

def hourly_mean(inarray):
    outarray=np.array([np.nanmean(inarray[i::24]) for i in range(24)])
    return outarray

def plot_ts(ERA5,Meso,Bias,unit,var,Stat_type,month,rmsd=None):
    if var in ['WBGT','TAIR','WS2M','WSPD']:
        legloc=4
    else: legloc=2
    hours=np.arange(len(ERA5))
    fig,ax1 = plt.subplots()
    ax2=ax1.twinx()
    #if var=='WS2M' and Stat_type=='Mean':
    #    ratio=Meso/ERA5
    #    print(ratio)
    #    print(np.mean(ratio))
    #    ax1.plot(hours,ratio,'m',label='ratio')
    if formula=="Liljegren":
        ax1.plot(hours,ERA5,'m-',label='ERA5 Liljegren')
        ax1.plot(hours,Meso,'b-',label='ERA5 Dimiceli')
    elif formula=="Dimiceli":
        ax1.plot(hours,ERA5,'b-',label='ERA5 Dimiceli')
        ax1.plot(hours,Meso,'r-',label='OK Mesonet')

    #if var=='PRES' and Stat_type=='Mean':
    #    plt.plot(hours,Bias*-1000,'g-',label='Diff*-1000')
    #else:
    if rb=='rmsd':
        rblabel="RMSD"
    else:
        rblabel="Difference"
    ax2.plot(hours,Bias,'g-',label=rblabel)
    if rb=='both':
        ax2.plot(hours,rmsd,'orange',label='RMSD')
    ax1.set_xlabel('Time (UTC)',fontsize=16)
    ax1.set_ylabel(unit,fontsize=16)
    ax2.set_ylabel('{} ({})'.format(rblabel,unit),fontsize=16)
    ax1.tick_params(axis='both', labelsize=12)
    ax2.tick_params(axis='both', labelsize=12)
    if Stat_type=="Mean":
        ax1.set_ylim(66,85)
        ax1.set_yticks(np.arange(66,85,2))
        if rb=='bias':
            ax2.set_ylim(-2.1,1.6)
            ax2.set_yticks(np.arange(-2,1.6,0.5))
        elif rb=='both':
            ax2.set_ylim(-2.1,2.2)
            ax2.set_yticks(np.arange(-2,2.2,0.5))
        else:
            ax2.set_ylim(0.7,2.2)
            ax2.set_yticks(np.arange(0.8,2.2,0.2))
    elif Stat_type=="STDev":
        ax1.set_ylim(2.8,4.1)
        ax1.set_yticks(np.arange(2.8,4.1,0.2))
        ax2.set_ylim(-0.26,0.31)
        ax2.set_yticks(np.arange(-0.25,0.31,0.05))
        legloc=2
    lines, labels=ax1.get_legend_handles_labels()
    lines2,labels2=ax2.get_legend_handles_labels()
    ax2.legend(lines+lines2,labels+labels2,loc=legloc,fontsize=12)
    #plt.legend()
    plt.title('{} hourly {} for month {}'.format(var,Stat_type,month),fontsize=20)
    plt.tight_layout()
    if formula=="Dimiceli": climopath="Climo_ws2m_SRADinterp_uncorrected"
    elif formula=="Liljegren": climopath="Climo-Liljegren"
    if rb=='rmsd':
        plt.savefig('/home/bjdavis/python/figures/{}/{}_{}_ts{}-rmsd_dpi{}.png'.format(climopath,var,Stat_type,month,str(figdpi)),dpi=figdpi)
    if rb=='both':
        plt.savefig('/home/bjdavis/python/figures/{}/{}_{}_ts{}-rmsd_bias_dpi{}.png'.format(climopath,var,Stat_type,month,str(figdpi)),dpi=figdpi)
    else:
        plt.savefig('/home/bjdavis/python/figures/{}/{}_{}_ts{}-bias_dpi{}.png'.format(climopath,var,Stat_type,month,str(figdpi)),dpi=figdpi)
    plt.show()

for month in range(start_month,end_month+1):
    m='{:02d}'.format(month)
    if formula=="Liljegren":
        infile_mean='{}/biasinfo-Liljegren.{}.ltm.nc'.format(inpath_bias,m)
        infile_std='{}/biasinfo-Liljegren.{}.ltstd.nc'.format(inpath_bias,m)
        infile_rmsd='{}/RMSDinfo-Liljegren.{}.nc'.format(inpath_rmsd,m)
    elif formula=="Dimiceli":
        infile_mean='{}/biasinfo_uncorrected_SRADinterp.{}.ltm.nc'.format(inpath_bias,m)
        infile_std='{}/biasinfo_uncorrected_SRADinterp.{}.ltstd.nc'.format(inpath_bias,m)
        infile_rmsd='{}/RMSDinfo-uncorrected-SRADinterp.{}.nc'.format(inpath_rmsd,m)

    if Stat_type=='Mean':
        Varn, Varnm, lat, lon=ReadNetCDF(infile_mean,Vars)
        print(len(Varn))
    elif Stat_type=='STDev':
        Varn, Varnm, lat, lon=ReadNetCDF(infile_std,Vars)
    eWBGTm,mWBGTm,bWBGT=Varn[0],Varn[1],Varn[2]

    print(eWBGTm)
    print(mWBGTm)
    print(bWBGT)
    #exit()

    if mesotype=='obs':
        eTDEWm,mTDEWm,bTDEW=Varn[3],Varn[4],Varn[5]
        eTAIRm,mTAIRm,bTAIR=Varn[6],Varn[7],Varn[8]
        eSRADm,mSRADm,bSRAD=Varn[9],Varn[10],Varn[11]
        eWS2Mm,mWS2Mm,bWS2M=Varn[12],Varn[13],Varn[14]
        ePRESm,mPRESm,bPRES=Varn[15]/100,Varn[16]/100,Varn[17]/100
        eWSPDm,mWSPDm,bWSPD=Varn[18],Varn[19],Varn[20]
        #eWS2Mm=eWS2Mm*1.25
        bWS2M=eWS2Mm-mWS2Mm

    eWBGTmh,mWBGTmh,bWBGTh=hourly_mean(eWBGTm),hourly_mean(mWBGTm),hourly_mean(bWBGT)
    if mesotype=='obs':
        eTDEWmh,mTDEWmh,bTDEWh=hourly_mean(eTDEWm),hourly_mean(mTDEWm),hourly_mean(bTDEW)
        eTAIRmh,mTAIRmh,bTAIRh=hourly_mean(eTAIRm),hourly_mean(mTAIRm),hourly_mean(bTAIR)
        eSRADmh,mSRADmh,bSRADh=hourly_mean(eSRADm),hourly_mean(mSRADm),hourly_mean(bSRAD)
        eWS2Mmh,mWS2Mmh,bWS2Mh=hourly_mean(eWS2Mm),hourly_mean(mWS2Mm),hourly_mean(bWS2M)
        ePRESmh,mPRESmh,bPRESh=hourly_mean(ePRESm),hourly_mean(mPRESm),hourly_mean(bPRES)
        eWSPDmh,mWSPDmh,bWSPDh=hourly_mean(eWSPDm),hourly_mean(mWSPDm),hourly_mean(bWSPD)
        print(ePRESmh)
        print(mPRESmh)
        print(bPRESh)

    if rb=='rmsd':
        Varn, Varnm, lat, lon=ReadNetCDF(infile_rmsd,RVars)
        bWBGTh=Varn[0]
        if mesotype=='obs':
            bTDEWh=Varn[1]
            bTAIRh=Varn[2]
            bSRADh=Varn[3]
            bWS2Mh=Varn[4]
            bPRESh=Varn[5]
            bWSPDh=Varn[6]
            print('RMSD:',Varn)
    elif rb=='both':
        Varn, Varnm, lat, lon=ReadNetCDF(infile_rmsd,RVars)
        bWBGThr=Varn[0]
        if mesotype=='obs':
            bTDEWhr=Varn[1]
            bTAIRhr=Varn[2]
            bSRADhr=Varn[3]
            bWS2Mhr=Varn[4]
            bPREShr=Varn[5]
            bWSPDhr=Varn[6]
            print('RMSD:',Varn)
        

    print(eWBGTmh)
    print(mWBGTmh)
    print(bWBGTh)
    if rb=='both':
        plot_ts(eWBGTmh,mWBGTmh,bWBGTh,var_units_plt[0],'WBGT',Stat_type,m,bWBGThr)
    else:
        plot_ts(eWBGTmh,mWBGTmh,bWBGTh,var_units_plt[0],'WBGT',Stat_type,m)
    exit()
    if mesotype=='obs':
        plot_ts(eTDEWmh,mTDEWmh,bTDEWh,var_units_plt[1],'TDEW',Stat_type,m)
        plot_ts(eTAIRmh,mTAIRmh,bTAIRh,var_units_plt[2],'TAIR',Stat_type,m)
        plot_ts(eSRADmh,mSRADmh,bSRADh,var_units_plt[3],'SRAD',Stat_type,m)
        plot_ts(eWS2Mmh,mWS2Mmh,bWS2Mh,var_units_plt[4],'WS2M',Stat_type,m)
        plot_ts(ePRESmh,mPRESmh,bPRESh,var_units_plt[5],'PRES',Stat_type,m)
        plot_ts(eWSPDmh,mWSPDmh,bWSPDh,var_units_plt[6],'WSPD',Stat_type,m)







import numpy as np
from read_mdf import *
from WBGT_func import Calc_WBGT_from_obs as Calc_WBGT
import glob


#mdf_files=["./mesonet/202109202100.mdf"]
year_start=1996
year_end=2020
month_start=1
month_end=12

input_path='/data/deluge/observations/OKMesonet'
output_path='/data/deluge/observations/OKMesonet/latlon'


geofile="./mesonet/geoinfo.csv"
geoinfo=np.loadtxt(geofile,dtype='str',delimiter=',')
StID_geo=geoinfo[::,1]
Lat_geo=geoinfo[::,7]
Lon_geo=geoinfo[::,8]

print('starting loop')
for y in range(year_start,year_end+1):
    #if y == 1996:
    #    month_start=4
    #else: month_start=1
    for month in range(month_start,month_end+1):
        m='{:02d}'.format(month)
        input_dir='{}/{}/{}'.format(input_path,y,m)
        print(input_dir)
        mdf_files=sorted(glob.glob(input_dir+'/*'))
        print(mdf_files)
        for mdf_file in mdf_files:
            file_name=mdf_file.split('/')[-1]
            print(file_name) 
            header,Var_names,data=read_mdf(mdf_file,"All")
            #missing_val=np.where(data<=-900)
            #data[missing_val]=np.nan
            #print(header)
            StID_mdf=data[::,0]

            ###################################
            ############ Calc WBGT ############
            ###################################
            print('Calculating WBGT')
            Tind=np.where(Var_names=='TAIR')[0]
            T=data[::,Tind].astype('float')

            RHind=np.where(Var_names=='RELH')[0]
            RH=data[::,RHind].astype('float')

            inSradind=np.where(Var_names=='SRAD')[0]
            inSrad=data[::,inSradind].astype('float')

            WS_2mind=np.where(Var_names=='WS2M')[0]
            WS_2m=data[::,WS_2mind].astype('float')

            pind=np.where(Var_names=='PRES')[0]
            p=data[::,pind].astype('float')

            print('STID,TAIR,RH,SRAD,WS2M,P')
            for row in range(len(StID_mdf)):
                print(StID_mdf[row],T[row],RH[row],inSrad[row],WS_2m[row],p[row])
            WBGT=Calc_WBGT(T,RH,inSrad,WS_2m,p,Moist_type='RH')
            #exit()
            #print(StID_mdf)
        
            #geofile="./mesonet/geoinfo.csv"
            #geoinfo=np.loadtxt(geofile,dtype='str',delimiter=',')
        
            #StID_geo=geoinfo[::,1]
            #Lat_geo=geoinfo[::,7]
            #Lon_geo=geoinfo[::,8]
        
            #for i in range(geoinfo.shape[0]):
            #    print(StID_geo[i],Lat_geo[i],Lon_geo[i])

            ##################################
            ########## Find lat/lon ##########
            ##################################
            print('finding lat/lon')
        
            mdf_lats=np.ones(StID_mdf.shape)*np.nan
            mdf_lons=np.ones(StID_mdf.shape)*np.nan
        
            for i in range(StID_mdf.shape[0]):
                StID=StID_mdf[i]
                StID_geo_ind=np.where(StID_geo==StID)[0]
                mdf_lats[i]=Lat_geo[StID_geo_ind]
                mdf_lons[i]=Lon_geo[StID_geo_ind]

            ##################################
            ######### Write new file #########
            ##################################
            print('writing file')
        
            mdf_lats.shape=(mdf_lats.shape[0],1)
            mdf_lons.shape=(mdf_lons.shape[0],1)
            WBGT.shape=(WBGT.shape[0],1)
            New_data=[WBGT,mdf_lats,mdf_lons]

            Var_names.shape=(1,Var_names.shape[0])
            New_vars=np.array([['TWBG','NLAT','ELON']])
        
            outfile='{}/{}/{}/{}'.format(output_path,y,m,file_name)
        
            Var_names,data=append_mdf(Var_names,data,New_vars,New_data,True,outfile,header)
        
            #outfile='./mesomod/202109202100.mdf'
            #write_mdf(outfile,header,Var_names,data)
        
            #for i in range(data.shape[0]):
            #    print(data[i,::])

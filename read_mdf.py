import numpy as np

def read_mdf(datafile,return_vars):#,vars_names)
    f=open(datafile,'r')
    header=f.readlines()
    header=header[0:2]
    f.close()

    data=np.loadtxt(datafile,skiprows=2,dtype='str')
    #Header=data[0:2,::]
    Var_names=data[0,::]
    data=data[1::,::]
    Vars_list=[]
    nan_inds=np.where((data=='-994') | (data=='-995') | (data=='-996') | (data=='-997') | (data=='-998') | (data=='-999'))
    data[nan_inds]=np.nan
    if return_vars=='All':
        return header,Var_names,data
    for var in return_vars:
        Var_loc=np.where(Var_names==var)[0]
        try:
            append_data=data[::,Var_loc].astype(float)
        except:
            append_data=data[::,Var_loc]
        Vars_list.append(append_data)
    print(len(Vars_list))
    return header,return_vars,Vars_list

def write_mdf(outfile, header, Var_names, data):
    data=np.append(Var_names,data,axis=0)
    data=data.astype("str")
#    for i in range(data.shape[0]):
#        print(data[i,::])
    #outfile='./mesomod/202109202100.mdf'
    for i in range(len(header)):
        header[i]=header[i].replace('\n','')
    np.savetxt(outfile,np.array(header),fmt='%s')
    f=open(outfile,'ba')
    np.savetxt(f,data,fmt='%s',delimiter='	')
    f.close()

def append_mdf(Var_names,data,New_Var_names,New_data,write=False,outfile=None,header=None):
    #Var_names.shape=(1,Var_names.shape[0])
    New_Var_names=np.array(New_Var_names)
    Var_names=np.append(Var_names,New_Var_names,axis=1)

#    print(type(New_data))
#    if data.shape[0] == New_data.shape[0]:
#        data=np.append(data,New_data,axis=1) 
#    print(data.shape)
#    print(New_data[0].shape)
    if isinstance(New_data,list):
        for i in range(len(New_data)):
#            print(New_data[i].shape)
            data=np.append(data,New_data[i],axis=1)       

    elif isinstance(New_data,np.ndarray):
        data=np.append(data,New_data,axis=1)

    if write:
        write_mdf(outfile,header,Var_names,data)

    return Var_names,data

#datafile='./mesonet/202109202005.mdf'
#return_vars=['STID','RELH','TAIR','PRES','SRAD','WS2M']
#[StID,RH,T,p,Srad,WS_2m]=read_mdf(datafile,return_vars)

#for i in range(len(StID)):
#    print(StID[i],RH[i],T[i],p[i],Srad[i],WS_2m[i])
#print(data.shape)
#print(Var_names)
#print(data)

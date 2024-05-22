from netCDF4 import Dataset
import numpy as np

def ReadNetCDF(infile,Vars,start=None,end=None):
    print(infile)
    nc = Dataset(infile,'r')
    Varn=[]
    Varnm=[]
    lat=[]
    lon=[]
    for i in range(len(Vars)):
        Variable=nc.variables[Vars[i]]
        if start==None or end==None:
            Varn.append(nc.variables[Vars[i]][:])
            try:
                missing_value=Variable.missing_value
                missing_inds=np.where(Varn[i].astype(float)==missing_value.astype(float))
                Varn[i][missing_inds]=np.nan
            except:                  
                print('no missing value')
            #Varn.append(nc.variables[Vars[i]][:])
        else:
            Varn.append(nc.variables[Vars[i]][start:end])
            try:
                missing_value=Variable.missing_value
                missing_inds=np.where(Varn[i].astype(float)==missing_value.astype(float))
                Varn[i][missing_inds]=np.nan
            except:
                pass
#            Varn.append(nc.variables[Vars[i]][start:end])
        try:
            Varnm.append=Varn[i].mask
        except:
            Varnm.append(np.zeros(Varn[i].shape, dtype=bool))
    try:
        lat=nc.variables['latitude']
    except:
        lat=np.nan
    try:
        lon=nc.variables['longitude']
    except:
        lon=np.nan
    #lat=np.multiply(lat,(180.0/np.pi))
    #lon=np.multiply(lon,(180.0/np.pi))
    return Varn, Varnm, lat, lon

def WriteNetCDF(Dims,DimsName,DimsUnits,DimsLN,Vars,VarsNames,VarsUnits,VarsLN,VarsDims,output_dir,filename):

   '''
    This function saves t/r wind and found tc lats and lons into a netcdf file

    Parameters
    ----------
    Dims : list
       List containing all of the values of the dimensions to write
    DimsName : list
       List containing all of the dimension names
    DimsUnits : list
       List containing all of the dimension units 
    DimsLN : list
       List containing all of the long names of the dimensions
    Vars : list
       List containing all of the values of output variables to write to file
    VarsNames: list
       List containing all of the names of the output variables
    VarsUnits : list
       List containing all of the units of the output variables
    VarsLN : list
       List containing all of the long names of the output variables
    VarsDims: list
       List containing the dimensions of all of the output variables
    output_dir : str
       Path to directory to write nc file to
    filename : str
       Name of nc file to write to

   Returns
   -------
   nothing

   Notes
   -----
   all missing values for each variable stored in the nc file will be np.nan
   for Dims lists the same index must correspond to data for the same variable
   for Vars lists the same index must correspond to data for the same variable

  '''

   #Put filepath together
   filepath = output_dir+"/"+filename
   ncfile = Dataset(filepath,mode="w",format='NETCDF4')

   #Create dimensions and corresponding variables
   for i in range(len(Dims)):
       dimi_dim = ncfile.createDimension(DimsName[i],len(Dims[i]))
       dimi_var = ncfile.createVariable(DimsName[i],np.float32,(DimsName[i],))
       dimi_var.units = DimsUnits[i]
       dimi_var.long_name = DimsLN[i]
       if DimsLN[i]=="time":
           dimi_var.calendar = "gregorian"
       print(Dims[i].shape)
       dimi_var[:] = Dims[i]

   #Write variables
   for j in range(len(Vars)):
       #pack data
       Vars[j][np.where(np.isnan(Vars[j]))]=-999

       dataMin=np.nanmin(Vars[j])
       dataMax=np.nanmax(Vars[j])
       offset=dataMin
       scale=(dataMax-dataMin)/((2**15)-1)

       print('Min',dataMin)
       print('Max',dataMax)
       print('offset',offset)
       print('scale',scale)
      
       #Vars[j][np.where(np.isnan(Vars[j]))]=-999
       packed_var=np.around((Vars[j]-offset)/scale)
       #packed_var[np.where(np.isnan(Vars[j]))]=np.nan
       print('packed calculated')
       packed_var=packed_var.astype(int)
       print('packed converted to int')
       
       data_var = ncfile.createVariable(VarsNames[j],np.int16,VarsDims[j])
       print(1)
       print(packed_var.shape)
       print(VarsNames[j])
       if len(VarsDims[j])==1:
           data_var[:] = packed_var
       elif len(VarsDims[j])==2:
           data_var[:,:] = packed_var
       elif len(VarsDims[j])==3:
           data_var[:,:,:] = packed_var#Vars[j]
       else:
           exit('need to add code for {} dimensions in WriteNetCDF'.format(len(VarsDims[j])))
       print(2)
       data_var.units = VarsUnits[j]
       print(3)
       data_var.long_name = VarsLN[j]
       print(4)
       data_var.add_offset = offset
       print(5)
       data_var.scale_factor = scale
       print(6)
       data_var.missing_value = -999
   #Close nc file
   ncfile.close()

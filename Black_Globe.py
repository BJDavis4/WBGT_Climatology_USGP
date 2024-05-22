import numpy as np
from Python_mods import *

#Lv=2.5e6 #J/kg
Lv=2.5e6
Rv=461.5 #J/kg/K
t0=273.15 #K
e0=6.11 #hPa
Rd=287 #hPa
Cp=1005 #J/kg/K
epsilon=Rd/Rv #unitlessimport numpy as np
from Python_mods import *

def Calc_T_Black_Globe(inTempC, inSrad, WS_2m):
    kCp=0.385
    kM=132.0
    kSurfaceArea=0.729659
    kEm=0.95
    kGrassAlbedo=0.5
    kDBeam=0.05
    diffBeam=0.95
    kSp=0.25
    ws2mMetersPerHour=WS_2m*3600

    calc1 = inSrad*kEm*kSurfaceArea
    q1=(calc1*kDBeam*kSp)+(calc1*diffBeam*(1.0+kGrassAlbedo))
    
    calc2= (q1/(kCp*kM))
    q2=q1-(0.115*(ws2mMetersPerHour**0.58)*calc2)

    Tbg=(q2/(kCp*kM))+inTempC
    return Tbg


#Tbg=Calc_T_Black_Globe(31.4,867,1.9)
#print(Tbg)

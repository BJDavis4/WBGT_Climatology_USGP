import numpy as np
from Python_mods import *
from Black_Globe import Calc_T_Black_Globe
from wet_bulb import *
import matplotlib
import matplotlib.pyplot as plt

def Calc_WBGT(T, Tw_nat, Tbg):
    Twbgt=(0.7*Tw_nat)+(0.2*Tbg)+(0.1*T)
    return Twbgt

def Calc_WBGT_from_obs(T,Moist,inSrad,WS_2m,p,Moist_type='RH'):
    Lv=2.5e6
    Rv=461.5 #J/kg/K
    t0=273.15 #K
    e0=6.11 #hPa
    Rd=287 #hPa
    Cp=1005 #J/kg/K
    epsilon=Rd/Rv #unitless

    if Moist_type=='RH':
        Td=Vapr_to_T(Vapr(T)*Moist*0.01) #Convert RH to Td
    elif Moist_type=='Td':
        Td=Moist
    
    #Calculate Natural Wet Bulb Temperature
    print('begin wet bulb')
#    print(inSrad,WS_2m,Td,T,p)
    Tw_nat=Calc_Tw_nat(inSrad=inSrad,WS_2m=WS_2m,Td=Td,T=T,p=p)
    print(Tw_nat)
    Tw_nat=Convert_T_Unit(Tw_nat,'C','F')
#    print('Wet Bulb:',Tw_nat)
    
    #Calculate Black Globe Temperature
    print('begin black globe')
    Tbg=Calc_T_Black_Globe(T,inSrad,WS_2m)
    Tbg=Convert_T_Unit(Tbg,'C','F')
#    print('Black Globe:',Tbg)
    
    #Calculate WBGT
    print('begin WBGT')
    T=Convert_T_Unit(T,'C','F')
#    print('T:',T,'Tw:',Tw_nat,'Tbg:',Tbg)
    Twbgt=Calc_WBGT(T,Tw_nat,Tbg)
#    print("WBGT",Twbgt)
    return Twbgt

#exit()
#plot relationship between Twbgt and WS2m
#plt.figure()
#plt.scatter(inSrad,Twbgt,color='b',label='WBGT')
#plt.scatter(inSrad,Tw_nat,color='r',label='Nat. Wet Bulb')
#plt.scatter(inSrad,Tbg,color='k',label='Black Globe')
#plt.legend(loc=0)
#plt.xlabel('Solar Radiation (W m^-2)')
#plt.ylabel('Temperature (F)')
#plt.ylim([70,120])
#plt.title('Sensitivity of WBGT to solar radiation')#im wind speed')
#plt.show()
##plt.savefig('srad-WBGT-25.png')
#plt.close()

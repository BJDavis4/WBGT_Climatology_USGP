import numpy as np
from Python_mods import Convert_moisture, Convert_T_Unit

#Lv=2.5e6 #J/kg
Lv=2.7e6
Rv=461.5 #J/kg/K
t0=273.15 #K
e0=6.11 #hPa
Rd=287 #hPa
Cp=1005 #J/kg/K
epsilon=Rd/Rv #unitless

T_univ=30.4
RH=53
Td_univ=Convert_moisture(RH,'RH','Percent','Td','C',Temp=T_univ,TempUnits='C')
#Tw=72
#Tw_nat=75
#WBGT=81
inSrad=709 #W/m^2
WS_2m=3.9 #m/s
p=966.75 #hPa




T=np.ones((1,1,1))*T_univ #K
Td=np.ones((1,1,1))*Td_univ #K
#p=974.66 #hPa

#def CC(T):
#    e=e0*np.exp((Lv/Rv)*((1.0/t0)-(1.0/T)))
#    return e

#def calc_Td(e):
#    T=((1.0/t0)-((Rv/Lv)*np.log(e/e0)))**(-1)
#    return T

def Vapr(T):
#    Tk=T+273.15
#    expo=23.832241-5.02808(np.log(Tk)/np.log(10.0))*(10**(11.334-(0.0303998*Tk)))
#    expo=expo+(8.1328e-3)*(10**(3.49149-(1302.8844/Tk)))-(2949.076/Tk)
#    ret=10**expo
#    print(T)
    try:
        denom=T.astype('float')+237.3
    except:
        denom=T+237.3
#    print(denom)
    expo=(7.5*T)/denom
#    print(expo,10**expo)
    vapr=6.11*(10.0**expo)
#    print(vapr)
    return vapr

def Vapr_to_T(e):
    d=np.log10(e/e0)
    b=237.3
    a=7.5
    T=(d*b)/(a-d)
    return T

def NWS_Theta_e(T,p,r,klcl):
    Rd=287.04
    cp=1005
    Theta_e=T*((1000.0/p)**((Rd/cp)*(1.0-(0.28*r))))*np.exp(((3376.0/klcl)-2.54)*r*(1.0+0.81*r))
    return Theta_e

def Calc_Tw_meso(Td,T,p):
    #note, Td and T must be entered in Celcius
    tguess=np.copy(T)
    shape=Td.shape
    currsign, prevsign= np.ones(shape),np.ones(shape)
    incr=np.ones(shape)*10.0
    thetadiff=np.ones(shape)
    estthetae=np.ones(shape)
#    print(Td)
    esurf=Vapr(Td)
#    print('esurf:',esurf)
    klcl = Td+273.15-(0.212+(0.001571*Td)-(0.000436*T))*(T-Td)
    mixr=621.97*esurf/(p-esurf)
    thetae=NWS_Theta_e(T+273.15,p,mixr/1000.0,klcl)
    debugcount=0
    while True:
        ktguess=tguess+273.15
        eguess=Vapr(tguess)
        mixrguess=621.97*eguess/(p-eguess)
        estthetae=NWS_Theta_e(ktguess,p,mixrguess/1000.0,klcl)
        thetadiff=thetae-estthetae
        if np.nanmax(abs(thetadiff))>0.05:
            currsign[np.where(thetadiff<0.0)]=-1
            currsign[np.where(thetadiff>=0.0)]=1
            #if thetadiff<0:
            #    currsign=-1
            #else: currsign=1
            #currsign=thetadiff/abs(thetadiff)
            signchange=np.where(currsign != prevsign)
            prevsign[signchange]=currsign[signchange]
            incr[signchange]=incr[signchange]*0.5
            #if (currsign != prevsign):
            #    prevsign=currsign
            #    incr=incr*0.5
            tguess=tguess+incr*prevsign

        else: break
#    print(tguess)
    naninds=np.where((p<=-900)|(np.isnan(p))|(T<=-900)|(np.isnan(p))|(Td<=-900)|(np.isnan(Td)))
    tguess[naninds]=np.nan
    return tguess

def Calc_Tw_psych(Td,T,p):
    e = T_to_e_2(Td)
    Tw=(T+Td)/2
    while True:
        eTw=T_to_e_2(Tw)
        Tw_new=T-((Lv*epsilon)/(Cp*p))*(eTw-e)
#        Tw_new=T-((Lv*epsilon)/((Cp*p)*(1+(0.0015*Tw))))*(eTw-e)
#        print(Cp/(Lv*epsilon)) 
#        print(Tw, Tw_new)
        if abs(Tw-Tw_new) < 0.001:
            Tw=Tw_new
            break
        else: Tw=(Tw+Tw_new)/2

    return Tw

def Calc_Tw_nat(Tw=None,inSrad=None,WS_2m=None,Td=-999,T=-999,p=-999):
#    print(inSrad,WS_2m,Td,T,p)
    if Tw==None:
        print('no Tw privided')
        Tw=Calc_Tw_meso(Td,T,p)
#    print(Tw)
#    print(Convert_T_Unit(Tw,'C','F'))
    try:
        if inSrad==None: exit('Must provide inSrad in Calc_TW_nat')
        if WS_2m==None: exit('Must provide WS_2m in Calc_TW_nat')
    except:
        pass
#    print(inSrad,WS_2m)
    Correction=(0.0021*inSrad)-(0.43*WS_2m)+1.93
#    print(Correction)
    try:
        Correction[np.where(Correction<0)]=0
    except:
        if Correction < 0: Correction = 0
    Tw_nat=Tw+Correction

    return Tw_nat

#T=Convert_T_Unit(T,'C','K')
#Td=Convert_T_Unit(Td,'C','K')

#Tw_nat=Calc_Tw_nat(inSrad=inSrad,WS_2m=WS_2m,Td=Td,T=T)
#Tw_nat=Convert_T_Unit(Tw_nat,'K','F')
#print(Tw_nat)











import numpy as np

Lv=2.5e6 #J/kg
Rv=461.5 #J/kg/K
t0=273.15 #K
e0=6.11 #hPa
Rd=287 #hPa
Cp=1005 #J/kg/K
epsilon=Rd/Rv #unitless

def calc_uk(z0,z1,u1):
    denom=np.log(z1/z0)
    uk=u1/denom
    return uk

def interp_ws_log(z,z0,z1,u1):
    uk=calc_uk(z0,z1,u1)
    u=uk*np.log(z/z0)
    return u

def interp_ws_pow(z,z1,u1,p):
    u=u1*((z/z1)**p)
    return u



def Calc_ws(u1,u2,u3=None):
    #C2=(u1**2)+(u2**2)
    C2=np.square(u1)+np.square(u2)
    if u3!=None:
        #C2=C2+(u3**2)
        C2=C2+np.square(u3)
    C=np.sqrt(C2)
    return C

def T_to_e(T):
    e=e0*np.exp((Lv/Rv)*((1.0/t0)-(1.0/T))) 
    return e

def T_to_e_2(T):
    e=e0*np.exp((17.67*T)/(T+243.5))
    return e

def e_to_T(e):
    T=((1.0/t0)-((Rv/Lv)*np.log(e/e0)))**(-1)
    return T

def e_to_T_2(e):
    a=17.67
    b=243.5
    d=np.log(e/e0)
    T=(b*d)/(a-d)
    return T

def q_to_w(q):
    w=q/(1-q)
    return w

def w_to_q(w):
    q=w/(1+w)
    return q

def w_to_e(w,p):
    e=(w*p)/(epsilon+w)
    return e

def e_to_w(e,p):
    w=((epsilon*e)/(p-e))
    return w

def q_to_e(q,p):
    e=(q*p)/epsilon
    return e

def e_to_q(e,p):
    q=(epsilon*e)/p
    return q

def RH_to_e(RH,es):
    e=RH*es
    return e

def e_to_RH(e,es):
    RH=e/es
    return RH

def Calc_Theta_e(T,p):
    p0=1000.00
    Theta_e=T*((p0/p)**(Rd/Cp))
    return Theta_e

def Convert_T_Unit(Temp,UnitIn,UnitOut):
    if UnitIn==UnitOut:
        return Temp
    elif UnitIn=='K':
        if UnitOut=='C':
            Temp=Temp-273.15
        elif UnitOut=='F':
            Temp=(1.8*Temp)-459.67
        else: exit('Invalid output unit in Convert_T_Units')
    elif UnitIn=='C':
        if UnitOut=='K':
            Temp=Temp+273.15
        elif UnitOut=='F':
            Temp=(1.8*Temp)+32.0
        else: exit('Invalid output unit in Convert_T_Units')
    elif UnitIn=='F':
        if UnitOut=='C':
            Temp=(Temp-32.0)/1.8
        elif UnitOut=='K':
            Temp=(Temp+459.67)/1.8
        else: exit('Invalid output unit in Convert_T_Units')
    else: exit('Invalid input units in Convert_T_Units')
    return Temp

def Convert_Pres_Units(Pres,UnitIn,UnitOut):
    if UnitIn==UnitOut:
        return Pres
    elif UnitIn=='hPa':
        if UnitOut=='Pa':
            Pres=Pres*100
        elif UnitOut=='kPa':
            Pres=Pres*0.1
        else: exit('Invalid output units in Convert_Pres_Units')
    elif UnitIn=='Pa':
        if UnitOut=='hPa':
            Pres=Pres*0.01
        elif UnitOut=='kPa':
            Pres=Pres*0.001
        else: exit('Invalid output units in Convert_Pres_Units')
    elif UnitIn=='kPa':
        if UnitOut=='Pa':
            Pres=Pres*1000
        elif UnitOut=='hPa':
            Pres=Pres*10
        else: exit('Invalid output units in Convert_Pres_Units')
    else: exit('Invalid input units in Convert_Pres_Units')
    return Pres

def Convert_qw_Unit(qw,UnitIn,UnitOut):
    if UnitIn=='gg':
        UnitIn='kgkg'
    if UnitOut=='gg':
        UnitOut='kgkg'
    if UnitIn==UnitOut:
        return qw
    elif UnitIn=='gkg':
        if UnitOut=='kgkg':
            qw*0.001
        else: exit('Invalid output units in Convert_qw_Unit')
    elif UnitIn=='kgkg':
        if unitOut=='gkg':
            qw=qw*1000
        else: exit('Invalid output units in Convert_qw_Unit')
    else: exit('Invalid input units in Convert_qw_Unit')    
    return qw

def Convert_Moist_Unit(Moist,MoistVar,UnitIn,UnitOut):
    if UnitIn==UnitOut:
        return Moist
    elif MoistVar=='Td':
        Moist = Convert_T_Unit(Moist,UnitIn,UnitOut)
    elif MoistVar=='e':
        Moist = Convert_Pres_Unit(Moist,UnitIn,UnitOut)
    elif MoistVar=='RH':
        if UnitIn=="Percent":
            if UnitOut=="Fraction":
                Moist=Moist*0.01
            else: exit('Invalid RH output units in Convert_Moist_Unit')
        elif UnitIn=="Fraction":
            if UnitOut=="Percent":
                Moist=Moist*100
            else: exit('Invalid RH output units in Convert_Moist_Unit')
        else: exit('Invalid RH input units in Convert_Moist_Unit')
    elif MostVar=='q':
        Moist = Convert_qw_Unit(Moist,UnitIn,UnitOut)
    elif MoistVar=='w':
        Moist = Convert_qw_Unit(Moist,UnitIn,UnitOut)
    else: exit('Invalid Input Variable in Convert_Moist_Unit')
    return Moist

def Convert_moisture(Moist, MoistVarIn, MoistUnitsIn, MoistVarOut, MoistUnitsOut, Temp=None, TempUnits=None, Pres=None, PresUnits=None):
    ValidMoist=['Td','e','RH','q','w']
    if MoistVarIn not in ValidMoist:
        exit('Invalid input moisture units')
    if MoistVarOut not in ValidMoist:
        exit('Invalid output moisture units')
    if MoistVarIn==MoistVarOut:
        if MoistUnitIn==MoistUnitOut:
            return Moist
        else:
            Moist=Convert_Moist_unit(Moist,MoistVarIn,MoistUnitsIn,MoistUnitsOut)
            return Moist

    if MoistVarIn=='q':
        Moist=Convert_Moist_Unit(Moist,MoistVarIn,MoistUnitsIn,'kgkg')
        if MoistUnitsOut=='w':
            Moist=q_to_w(Moist)
        else:
            if Pres==None:
                exit('Pressure required')
            elif PresUnits==None:
                exit('pressure units required')
            else:
                Pres=Convert_Pres_Unit(Pres,PresUnits,'hPa')
                PresUnits='hPa'
            Moist=q_to_e(Moist,Pres)
            MoistVarIn='e'
            MoistUnitsIn='hPa'
            if MoistVarOut=='e':
                Moist=Convert_Pres_Unit(Moist,MoistUnitsIn,MoistUnitsOut)
                return Moist

    if MoistVarIn=='w':
        Moist=Convert_Moist_Unit(Moist,MoistVarIn,MoistUnitsIn,'kgkg')
        if MoistUnitsOut=='q':
            Moist=w_to_q(Moist)
        else:
            if Pres==None:
                exit('Pressure required')
            elif PresUnits==None:
                exit('pressure units required')
            else:
                Pres=Convert_Pres_Unit(Pres,PresUnitsIn,'hPa')
                PresUnits='hPa'
            Moist=w_to_e(Moist,Pres)
            MoistVarIn='e'
            MoistUnitsIn='hPa'
            if MoistVarOut=='e':
                Moist=Convert_Pres_Unit(Moist,MoistUnitsIn,MoistUnitsOut)
                return Moist

    if MoistVarIn=='RH':
        Moist=Convert_Moist_Unit(Moist,MoistVarIn,MoistUnitsIn,'Fraction')
        if TempUnits==None:
            exit('Temperature units required')
        try:
            if Temp==None:
                exit('Temperature required for RH conversion')
            else:
                Temp=Convert_T_Unit(Temp,TempUnits,'C')
                TempUnits='C'
        except:
            Temp=Convert_T_Unit(Temp,TempUnits,'C')
            TempUnits='C'

        es=T_to_e_2(Temp)
        Moist=RH_to_e(Moist,es)
        MoistVarIn='e'
        MoistUnitsIn='hPa'
        if MoistVarOut=='e':
            Moist=Convert_Pres_Unit(Moist,MoistUnitsIn,MoistUnitsOut)
            return Moist
        print('e: ', Moist)
        print('es: ', es)

    if MoistVarIn=='Td':
        Moist=Convert_Moist_Unit(Moist,MoistVarIn,MoistUnitsIn,'K')
        Moist=T_to_e(Moist)
        MoistVarIn='e'
        MoistUnitsIn='hPa'                
        if MoistVarOut=='e':
            Moist=Convert_Pres_Unit(Moist,MoistUnitsIn,MoistUnitsOut)
            return Moist
        #else: e=Moist

    if MoistVarIn=='e':
        Moist=Convert_Moist_Unit(Moist,MoistVarIn,MoistUnitsIn,'hPa')
        if MoistVarOut=='q':
            Moist=e_to_q(Moist,Pres)
            Moist=Convert_Moist_Unit(Moist,'q','kgkg',MoistUnitsOut)
        
        if MoistVarOut=='w':
            Moist=e_to_w(Moist,Pres)
            Moist=Convert_Moist_Unit(Moist,'w','kgkg',MoistUnitsOut)

        if MoistVarOut=='Td':
            Moist=e_to_T_2(Moist)
            #Moist=Convert_Moist_Unit(Moist,'Td','K',MoistUnitsOut)
            Moist=Convert_Moist_Unit(Moist,'Td','C',MoistUnitsOut)

        if MoistVarOut=='RH':
            try:
                if Temp==None:
                    exit('Temperature required')
                elif TempUnits==None:
                    exit('Temperature units required')
                else:
                    Temp=Convert_T_Unit(Temp,TempUnits,'K')
                    TempUnits='K'    
            except:
                Temp=Convert_T_Unit(Temp,TempUnits,'K')
                TempUnits='K'
            es=T_to_e(Temp)
            Moist=e_to_RH(Moist,es)
            Moist=Convert_Moist_Unit(Moist,'RH','Fraction',MoistUnitsOut)

    return Moist


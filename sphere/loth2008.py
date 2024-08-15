import numpy as np
import math

def C08(M):
    if M <= 1.45:
        Cm = 5/3 + (2/3)*np.tanh(3*np.log(M+0.1))
    else:
        Cm = 2.044 + 0.2*np.exp(-1.8*(np.log(M/1.5))**2)
    return Cm

def G08(M):
    if M < 0.89:
        Gm = 1 - 1.525*M**4
    else:
        Gm = 0.0002 + 0.0008*np.tanh(12.77*(M-2.02))
    return Gm

def H08(M):
    Hm = 1 - (0.258*C08(M))/(1+514*G08(M))
    return Hm

def f08(Kn):
    fKn = 1/(1+Kn*(2.514+0.8*np.exp(-0.55/Kn)))
    return fKn

def CKR08(Kn,Re):
    CdKnRe = (24/Re)*(1+0.15*Re**0.687)*f08(Kn)
    return CdKnRe

def Cf08(M,Re,T_p,T_inf,gam,Kn):
    s = M*np.sqrt(gam/2)
    
    Cdfm = (1+2*s**2)*np.exp(-(s**2))/(s**3*np.sqrt(np.pi)) + \
            (4*s**4 + 4*s**2 - 1)*math.erf(s) / (2*s**4) + \
            2/(3*s)*np.sqrt(np.pi*T_p/T_inf)
    return Cdfm

def CfR08(M,Re,T_p,T_inf,gam,Kn):
    Cdfm = Cf08(M,Re,T_p,T_inf,gam,Kn)
    CdfmRe = Cdfm / (1+(Cdfm/1.63-1)*np.sqrt(Re/45))
    return CdfmRe
    
def Loth_2008_Cd(M,Re,T_p,T_inf,gam,Kn):
    
    if not isinstance(Re, np.ndarray):
        Re = np.asarray([Re])
    
    Cd=np.zeros(np.shape(Re))
    for rei in np.arange(len(Re)):
        if Re[rei] > 45:
            Cd[rei] = (24/Re[rei])*(1.0+0.15*Re[rei]**0.687)*H08(M) + (0.42*C08(M))/(1+(42500*G08(M))/Re[rei]**1.16)
        else:
            Cd[rei] = CKR08(Kn,Re[rei])/(1+M**4) + M**4*CfR08(M,Re[rei],T_p,T_inf,gam,Kn)/(1+M**4)
    return Cd
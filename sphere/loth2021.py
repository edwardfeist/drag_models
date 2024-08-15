import numpy as np
import math

def C21(M):
    if (M < 1.5):
        Cm = 1.65 + 0.65*np.tanh(4.0*M - 3.4)
    else:
        Cm = 2.18 - 0.13*np.tanh(0.9*M - 2.7)
    return Cm

def G21(M):
    if (M < 0.8):
        Gm = 166*M**3.0 + 3.29*M**2.0 - 10.9*M + 20.0
    else:
        Gm = 5.0 + 40.0*M**(-3.0)
    return Gm

def H21(M):
    if (M < 1.0):
        Hm = 0.0239*M**3.0 + 0.212*M**2.0 - 0.074*M + 1.0
    else:
        Hm = 0.93 + 1/(3.5+M**5.0)
    return Hm

def f21(Kn):
    fKn = 1.0/(1+Kn*(2.514+0.8*np.exp(-0.55/Kn)))
    return fKn

def CKR21(Kn,Re):
    CdKnRe = (24.0/Re)*(1.0 + 0.15*Re**0.687)*f21(Kn)
    return CdKnRe

def Cf21(M,Re,T_p,T_inf,gam,Kn):
    s = M*np.sqrt(gam/2.0)  
    Cdfm = (1.0 + 2.0*s**2.0)*np.exp(-(s**2.0))/(s**3.0*np.sqrt(np.pi)) + \
            (4.0*s**4.0 + 4.0*s**2.0 - 1.0)*math.erf(s) / (2.0*s**4.0) + \
            2.0/(3.0*s)*np.sqrt(np.pi)
    return Cdfm

def J21(M):
    if (M <= 1.0):
       Jm = 2.26 - 0.10/M + 0.14/M**3.0
    else:
       Jm = 1.60 + 0.25/M + 0.11/M**2.0 + 0.44/M**3.0
    return Jm

def CfR21(M,Re,T_p,T_inf,gam,Kn):
   Cdfm = Cf21(M,Re,T_p,T_inf,gam,Kn)
   CdfmRe = Cdfm / (1.0+(Cdfm/J21(M) - 1.0)*np.sqrt(Re/45.0))
   return CdfmRe

def Loth_2021_Cd(M,Re,T_p,T_inf,gam,Kn): 
    if not isinstance(Re, np.ndarray):
        Re = np.asarray([Re])
    
    Cd=np.zeros(np.shape(Re))
    for rei in np.arange(len(Re)):
        if (Re[rei] > 45.0):
            Cd[rei] = (24.0/Re[rei])*(1.0+0.15*Re[rei]**0.687)*H21(M) + \
                   (0.42*C21(M))/(1.0 + (42500.0/Re[rei]**(1.16*C21(M))) + G21(M)/Re[rei]**0.5)
        else:
            Cd[rei] = CKR21(Kn,Re[rei])/(1.0+M**4.0) + M**4.0*CfR21(M,Re[rei],T_p,T_inf,gam,Kn)/(1.0+M**4.0)
    return Cd
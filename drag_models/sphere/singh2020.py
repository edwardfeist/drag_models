import numpy as np
import math

delta_0   = 9.4
alpha_0   = 0.356
A1        = 2.514
A2        = 0.8
A3        = 0.55
eta       = 1.8
alpha_hoc = 1.27
Cdminf    = 0.9
#C0        = 24 / 9.06**2
C0        = 24 / delta_0**2
omega     = 0.74
e         = 0.9

def theta_func(M, gam):
    theta = (1 + (gam-1)*(M**2/2))**(gam/(gam-1))
    return theta

def C1_func(M, gam):
    C1 = (Cdminf - C0 * (1 + (gam-1)**2 / (4*gam)) ** (gam / (gam-1))) / (1 - (1/(alpha_0*M))*(gam-1)/(gam+1))
    return C1

def Ccd_func(M, Re, gam, R, T_inf, Ts):
    alpha = 1 / (alpha_0 * M + 1 - alpha_0)
    C1 = C1_func(M, gam)
    Uinf = M * np.sqrt(gam*R*T_inf)
    
    if (M > 1):
        M1s = M**2
        gm1 = gam - 1
        M2s = (gm1*M1s + 2) / (2*gam*M1s - gm1)
        M2 = np.sqrt(M2s)
        Us = M2 * np.sqrt(gam * R * Ts)
        
        theta_s = theta_func(M2, gam)
        Re_s = Re * (1 /alpha**2 * T_inf / Ts)**omega * theta_s**((gam+1/(2*gam)) - ((gam - 1)/gam)*omega)
        
        Ccd = C1 * (1-alpha * Us / Uinf) + C0 * theta_s * (1 + delta_0 / np.sqrt(Re_s))**2
        
    else:
        theta_s = theta_func(M, gam)
        Re_s    = Re
        
        Ccd = C0 * theta_s * (1 + delta_0 / np.sqrt(Re_s))**2
        
    return Ccd

def Cdfm_func(M,Re,T_p,T_inf,gam,Kn):
    s = M*np.sqrt(gam/2.0)  
    Cdfm = (1.0 + 2.0*s**2.0)*np.exp(-(s**2.0))/(s**3.0*np.sqrt(np.pi)) + \
            (4.0*s**4.0 + 4.0*s**2.0 - 1.0)*math.erf(s) / (2.0*s**4.0) + \
            (2.0*(1-e))/(3.0*s)*np.sqrt(np.pi * T_p / T_inf)
    return Cdfm

def fKnW_func(M, Re, Kn, T_p, Ts):
    Wr = M**(2*omega) / Re
    WTr = Wr * (1+T_p / Ts)**omega
    
    fKn = 1.0/(1+Kn*(A1+A2*np.exp(-A3/Kn))) * 1 / (1 + alpha_hoc * WTr)
    return fKn   
    
def Br_func(M, Re, gam, T_p, T_inf, Ts):
    
    Wr = M**(2*omega) / Re
    WTr = Wr * (1+T_p / Ts)**omega
    
    Br = WTr * (M**(2*omega-1) + 1) / M**(2*omega-1)
    return Br

def Singh_2020_Cd(M,Re,T_p,T_inf,gam,Kn,R): 
    if not isinstance(Re, np.ndarray):
        Re = np.asarray([Re])
        
    if (M > 1):
        M1s = M**2
        gm1 = gam - 1
        gp1 = gam + 1
        Ts = T_inf * ((2*gam*M1s - gm1)*(gm1*M1s+2)) / (gp1**2 * M1s)
    else:
        Ts = T_inf
    
    Cd=np.zeros(np.shape(Re))
    for rei in np.arange(len(Re)):
        Ccd  = Ccd_func(M, Re[rei], gam, R, T_inf, Ts)
        Br   = Br_func(M, Re[rei], gam, T_p, T_inf, Ts)
        fKnW = fKnW_func(M, Re[rei], Kn, T_p, Ts)
        Cdfm = Cdfm_func(M, Re[rei], T_p, T_inf, gam, Kn)
        
        Cd[rei] = Ccd * fKnW / (1+Br**eta) + Cdfm * Br**eta / (1+Br**eta)
    return Cd
import numpy as np

def low_mach(M,Re,T_rat,gam):
    S = M*np.sqrt(0.5*gam)

    # first term
    r1 = 4.33 + (3.65 - 1.53*T_rat)/(1.0 + 0.353*T_rat)
    r1 = 24.0/(Re + S*r1*np.exp(-0.247*Re/S))

    # second term
    r2 = 0.03*Re + 0.48*np.sqrt(Re)
    r2 = (4.5 + 0.38*r2)/(1.0 + r2) + 0.1*M**2.0 + 0.2*M**8.0
    r2 = np.exp(-0.5*M/Re)*r2

    # third term
    r3 = 0.6*S*(1.0 - np.exp(-M/Re))

    Cd = r1 + r2 + r3
    
    return Cd

def high_mach(M,Re,T_rat,gam):
    S = M*np.sqrt(0.5*gam)

    r1 = 1.86*np.sqrt(M/Re)
    r2 = 2.0 + 2.0/(S**2.0) + (1.058/S)*np.sqrt(T_rat) - 1.0/(S**4.0)

    Cd = (0.9 + 0.34/M**2.0 + r1*r2)/(1.0 + r1) # note: 0.34/M^2
    
    return Cd
    
def Henderson_Cd(M,Re,T_p,T_inf,gam,Kn):

    if (M>1.75):
       Cd = high_mach(M,Re,T_p/T_inf,gam)
    elif (M<1.0):
       Cd = low_mach(M,Re,T_p/T_inf,gam)
    else: # linear interpolation
       Cd_low = low_mach(1.0,Re,T_p/T_inf,gam) 
       Cd_high = high_mach(1.75,Re,T_p/T_inf,gam)
       Cd = Cd_low + (4.0/3.0)*(M - 1.0)*(Cd_high - Cd_low)

    return Cd
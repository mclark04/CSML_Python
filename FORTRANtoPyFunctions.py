#The uppercase MAX function is deinfed here to mimic FORTRAN MAX which max() doesn't capture
import numpy as np

def MAX(var,maxVar): 
    tmp=np.max(var)    
    outVar=np.max((tmp,maxVar))
    return outVar

def MIN(var,maxVar): 
    tmp=np.min(var)
    outVar=np.min((tmp,maxVar))
    return outVar

#The FORTRAN SIGN function assignes the sign of B on abs(A)
def SIGN(A,B):
    if B>=0:
        Aout=abs(A)
    else:
        Aout=-abs(A)
    return Aout
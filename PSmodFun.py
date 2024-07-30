import numpy as np
import scipy as sp
import math
from scipy import stats
import emcee
import camb
from camb import model 
from PSwincon import *



 




 




#==============================================================================
#                             Integration of linear Power 
# 2. Calculate integrations:---------------------------------------------------    
def Fmn(m,n,r,x,y):
    if(m==0) and (n==0):
        FMN=(7.0*x + 3.0*r - 10.0*r*x*x)**2/(14.0*14.0*r*r*y*y*y*y);
    if(m==0) and (n==1):
        FMN= (7.0*x + 3.0*r - 10.0*r*x*x)*(7.0*x - r - 6.0*r*x*x)/(14.0*14.0*r*r*y*y*y*y);
    if(m==0) and (n==2):
        FMN= (x*x - 1.0)*(7.0*x + 3.0*r - 10.0*r*x*x)/(14.0*r*y*y*y*y);
    if(m==0) and (n==3):
        FMN= (1.0 - x*x)*(3.0*r*x - 1.0)/(r*r*y*y);       
    if(m==1) and (n==0):
        FMN= x*(7.0*x + 3.0*r - 10.0*r*x*x)/(14.0*r*r*y*y);
    if(m==1) and (n==1):
        FMN= (7.0*x - r - 6.0*r*x*x)**2/(14.0*14.0*r*r*y*y*y*y);
    if(m==1) and (n==2):
        FMN= (x*x - 1.0)*(7.0*x - r - 6.0*r*x*x)/(14.0*r*y*y*y*y);
    if(m==1) and (n==3):
        FMN=( 4.0*r*x + 3.0*x*x - 6.0*r*x*x*x - 1.0)/(2.0*r*r*y*y);             
    if(m==2) and (n==0):
        FMN= (2.0*x + r - 3.0*r*x*x)*(7.0*x + 3.0*r - 10.0*r*x*x)/(14.0*r*r*y*y*y*y);
    if(m==2) and (n==1):
        FMN= (2.0*x + r - 3.0*r*x*x)*(7.0*x - r - 6.0*r*x*x)/(14.0*r*r*y*y*y*y);
    if(m==2) and (n==2):
        FMN= x*(7.0*x - r - 6.0*r*x*x)/(14.0*r*r*y*y);
    if(m==2) and (n==3):
        FMN= 3.0*(1.0-x*x)*(1.0-x*x)/(y*y*y*y); 
    if(m==3) and (n==0):
        FMN= (1.0 - 3.0*x*x - 3.0*r*x + 5.0*r*x*x*x)/(r*r*y*y);
    if(m==3) and (n==1):
        FMN=  (1.0 - 2*r*x)*(1.0 - x*x)/(2.0*r*r*y*y);
    if(m==3) and (n==2):
        FMN=  (1.0 - x*x)*(2.0 - 12.0*r*x - 3.0*r*r + 15.0*r*r*x*x)/(r*r*y*y*y*y);
    if(m==3) and (n==3):
        FMN=  (-4.0 + 12.0*x*x + 24.0*r*x - 40.0*r*x*x*x + 3.0*r*r - 30.0*r*r*x*x + 35.0*r*r*x*x*x*x)/(r*r*y*y*y*y);        
    return FMN

def Imn_inner_integ(x,m,n,k,r,Pl_spline):
    y    = np.sqrt(1.0+r*r-2.0*r*x)#x=cos theta
    Plkq = sp.interpolate.splev(np.log10(k*y), Pl_spline, der=0)
    return Fmn(m,n,r,x,y)*10.0**(Plkq)

def Imn_outer_integ(r,m,n,k,Pl_spline):
    integ,errin= sp.integrate.quad(Imn_inner_integ,-1.0,1.0,epsabs=0.0,epsrel=1.0e-4,args=(m,n,k,r,Pl_spline))
    Pl   = sp.interpolate.splev(np.log10(k*r),Pl_spline, der=0)
    return r*r*integ*10.0**(Pl)

def I_mn(kmod,Pl_spline,ks):
    Imn=np.zeros((len(ks),4,4))
    for i in range(len(ks)):
        I00,errin= sp.integrate.quad(Imn_outer_integ,kmod[0],kmod[-1],args=(0,0,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        Imn[i,0,0]= I00*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)
        I01,errin= sp.integrate.quad(Imn_outer_integ,kmod[0],kmod[-1],args=(0,1,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        Imn[i,0,1]= I01*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)
        I02,errin= sp.integrate.quad(Imn_outer_integ,kmod[0],kmod[-1],args=(0,2,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        Imn[i,0,2]= I02*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)
        I03,errin= sp.integrate.quad(Imn_outer_integ,kmod[0],kmod[-1],args=(0,3,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        Imn[i,0,3]= I03*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)

        I10,errin= sp.integrate.quad(Imn_outer_integ,kmod[0],kmod[-1],args=(1,0,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        Imn[i,1,0]= I10*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)
        I11,errin= sp.integrate.quad(Imn_outer_integ,kmod[0],kmod[-1],args=(1,1,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        Imn[i,1,1]= I11*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)
        I12,errin= sp.integrate.quad(Imn_outer_integ,kmod[0],kmod[-1],args=(1,2,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        Imn[i,1,2]= I12*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)
        I13,errin= sp.integrate.quad(Imn_outer_integ,kmod[0],kmod[-1],args=(1,3,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        Imn[i,1,3]= I13*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)

        I20,errin= sp.integrate.quad(Imn_outer_integ,kmod[0],kmod[-1],args=(2,0,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        Imn[i,2,0]= I20*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)
        I21,errin= sp.integrate.quad(Imn_outer_integ,kmod[0],kmod[-1],args=(2,1,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        Imn[i,2,1]= I21*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)
        I22,errin= sp.integrate.quad(Imn_outer_integ,kmod[0],kmod[-1],args=(2,2,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3)
        Imn[i,2,2]= I22*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)
        I23,errin= sp.integrate.quad(Imn_outer_integ,kmod[0],kmod[-1],args=(2,3,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        Imn[i,2,3]= I23*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)

        I30,errin= sp.integrate.quad(Imn_outer_integ,kmod[0],kmod[-1],args=(3,0,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        Imn[i,3,0]= I30*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)
        I31,errin= sp.integrate.quad(Imn_outer_integ,kmod[0],kmod[-1],args=(3,1,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3)  
        Imn[i,3,1]= I31*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)
        I32,errin= sp.integrate.quad(Imn_outer_integ,kmod[0],kmod[-1],args=(3,2,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        Imn[i,3,2]= I32*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)
        I33,errin= sp.integrate.quad(Imn_outer_integ,kmod[0],kmod[-1],args=(3,3,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        Imn[i,3,3]= I33*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)
    return Imn

#--------------------------------------------------------------------------
def Gmn(m,n,r):    
    if(m==0) and (n==0):
        GMN= (12.0/(r*r) - 158.0 + 100.0*r*r - 42.0*r*r*r*r + (3.0/(r*r*r))*(r*r - 1.0)*(r*r - 1.0)*(r*r - 1.0)*(7.0*r*r + 2.0)*np.log((r + 1.0)/np.abs(r - 1.0)))/3024.0;
    if(m==0) and (n==1):
        GMN= (24.0/(r*r) - 202.0 + 56.0*r*r - 30.0*r*r*r*r + (3.0/(r*r*r))*(r*r - 1.0)*(r*r - 1.0)*(r*r - 1.0)*(5.0*r*r + 4.0)*np.log((r + 1.0)/np.abs(r - 1.0)))/3024.0;
    if(m==0) and (n==2):
        GMN= (2.0*(r*r + 1.0)*(3.0*r*r*r*r - 14.0*r*r + 3.0)/(r*r) - (3.0/(r*r*r))*(r*r - 1.0)*(r*r - 1.0)*(r*r - 1.0)*(r*r - 1.0)*np.log((r + 1.0)/np.abs(r - 1.0)))/224.0;
    if(m==1) and (n==0):
        GMN= (-38.0 +48.0*r*r - 18.0*r*r*r*r + (9.0/r)*(r*r - 1.0)*(r*r - 1.0)*(r*r - 1.0)*np.log((r + 1.0)/np.abs(r - 1.0)))/1008.0;
    if(m==1) and (n==1):
        GMN= (12.0/(r*r) - 82.0 + 4.0*r*r - 6.0*r*r*r*r + (3.0/(r*r*r))*(r*r - 1.0)*(r*r - 1.0)*(r*r - 1.0)*(r*r + 2.0)*np.log((r + 1.0)/np.abs(r - 1.0)))/1008.0;
    if(m==2) and (n==0):
        GMN= (2.0*(9.0 - 109.0*r*r + 63.0*r*r*r*r - 27.0*r*r*r*r*r*r)/(r*r) + (9.0/(r*r*r))*(r*r - 1.0)*(r*r - 1.0)*(r*r - 1.0)*(3.0*r*r + 1.0)*np.log((r + 1.0)/np.abs(r - 1.0)))/672.0;
    return GMN

def Jmn_integ(q,m,n,k,Pl_spline):
    r    = q/k
    if(r==1.0):
        r=(q+0.00001*q)/k
    Pl   = sp.interpolate.splev(np.log10(q),Pl_spline, der=0)
    return Gmn(m,n,r)*10.0**(Pl)

def J_mn(kmod,Pl_spline,ks):
    Jmn=np.zeros((len(ks),2,3))
    for i in range(len(ks)):
        J00,errin= sp.integrate.quad(Jmn_integ,kmod[0],kmod[-1],args=(0,0,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3)
        Jmn[i,0,0] = J00/(2.0*math.pi*math.pi)
        J01,errin= sp.integrate.quad(Jmn_integ,kmod[0],kmod[-1],args=(0,1,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3)
        Jmn[i,0,1] = J01/(2.0*math.pi*math.pi)
        J02,errin= sp.integrate.quad(Jmn_integ,kmod[0],kmod[-1],args=(0,2,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        Jmn[i,0,2] = J02/(2.0*math.pi*math.pi)

        J10,errin= sp.integrate.quad(Jmn_integ,kmod[0],kmod[-1],args=(1,0,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        Jmn[i,1,0] = J10/(2.0*math.pi*math.pi)
        J11,errin= sp.integrate.quad(Jmn_integ,kmod[0],kmod[-1],args=(1,1,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        Jmn[i,1,1] = J11/(2.0*math.pi*math.pi)
        J20,errin= sp.integrate.quad(Jmn_integ,kmod[0],kmod[-1],args=(2,0,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3)
        Jmn[i,1,2] = J20/(2.0*math.pi*math.pi)
    return Jmn

#--------------------------------------------------------------------------
def hmn(m,n,s,r,x,y):
    if (m == 0) and (n == 0):
        if (s == 0):
            numer = 7.0*x + 3.0*r - 10.0*r*x*x;
            return numer/(14.0*r*y*y);
        if (s == 1):
            numer = (7.0*x + 3.0*r - 10.0*r*x*x)*(3.0*x*x - 4.0*r*x + 2.0*r*r - 1.0);
            return numer/(14.0*3.0*r*y*y*y*y);
    if (m == 0) and (n == 1):
        if (s == 0):
            return 1.0;
        if (s == 1):
            numer = (3.0*x*x - 4.0*r*x + 2.0*r*r - 1.0);
            return (numer*numer)/(3.0*3.0*y*y*y*y);
    if (m == 0) and (n == 2):
        if (s == 1):
            numer = (3.0*x*x - 4.0*r*x + 2.0*r*r - 1.0);
            return numer/(3.0*y*y);
    if (m == 1) and (n == 0):
        if (s == 0):
            numer = 7.0*x - r - 6.0*r*x*x;
            return numer/(14.0*r*y*y);
        if (s == 1):
            numer = (7.0*x - r - 6.0*r*x*x)*(3.0*x*x - 4.0*r*x + 2.0*r*r - 1.0);
            return numer/(14.0*3.0*r*y*y*y*y);
    if (m == 1) and (n == 1):
        if (s == 0):
            return x/r;
        if (s == 1):
            numer = (3.0*x*x - 4.0*r*x + 2.0*r*r - 1.0)*x/r;
            return numer/(3.0*y*y);
    if (m == 2) and (n == 0):
       if (s == 0):
           numer = (x*x - 1.0)/2.0
           return numer/(y*y);
       if (s == 1):
           numer = (x*x - 1.0)*(3.0*x*x - 4.0*r*x + 2.0*r*r - 1.0);
           return numer/(6.0*y*y*y*y); 
    if (m == 3) and (n == 0):
       if (s == 0):
           numer = (2.0*x + r - 3.0*r*x*x)/2.0;
           return numer/(r*y*y);
       if (s == 1):
           numer = (2.0*x + r - 3.0*r*x*x)*(3.0*x*x - 4.0*r*x + 2.0*r*r - 1.0);
           return numer/(6.0*r*y*y*y*y); 

def Kmn_inner_integ(x,m,n,s,k,r,Pl_spline):
    y    = np.sqrt(1.0+r*r-2.0*r*x)
    Plkq = sp.interpolate.splev(np.log10(k*y), Pl_spline, der=0)
    return hmn(m,n,s,r,x,y)*10.0**(Plkq)

def Kmn_outer_integ(r,m,n,s,k,Pl_spline):
    integ,errin= sp.integrate.quad(Kmn_inner_integ,-1.0,1.0,epsabs=0.0,epsrel=1.0e-3,args=(m,n,s,k,r,Pl_spline))
    Pl   = sp.interpolate.splev(np.log10(k*r),Pl_spline, der=0)
    return r*r*integ*10.0**(Pl)

def K_mn(kmod,Pl_spline,ks):
    KMN=np.zeros((len(ks),3,2))
    KSmn=np.zeros((len(ks),7))
    for i in range(len(ks)):
        K00,errin = sp.integrate.quad(Kmn_outer_integ,kmod[0],kmod[-1],args=(0,0,0,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        KMN[i,0,0]= K00*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi) 
        Ks00,errin= sp.integrate.quad(Kmn_outer_integ,kmod[0],kmod[-1],args=(0,0,1,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        KSmn[i,0] = Ks00*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi) 
        
        K01,errin = sp.integrate.quad(Kmn_outer_integ,kmod[0],kmod[-1],args=(0,1,0,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        KMN[i,0,1]= K01*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi) 
        Ks01,errin= sp.integrate.quad(Kmn_outer_integ,kmod[0],kmod[-1],args=(0,1,1,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        KSmn[i,1] = Ks01*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi) 
        Ks02,errin= sp.integrate.quad(Kmn_outer_integ,kmod[0],kmod[-1],args=(0,2,1,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        KSmn[i,2] = Ks02*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi) 

        K10,errin = sp.integrate.quad(Kmn_outer_integ,kmod[0],kmod[-1],args=(1,0,0,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        KMN[i,1,0]= K10*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)
        Ks10,errin= sp.integrate.quad(Kmn_outer_integ,kmod[0],kmod[-1],args=(1,0,1,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        KSmn[i,3] = Ks10*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)
        
        K11,errin = sp.integrate.quad(Kmn_outer_integ,kmod[0],kmod[-1],args=(1,1,0,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        KMN[i,1,1]= K11*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi) 
        Ks11,errin= sp.integrate.quad(Kmn_outer_integ,kmod[0],kmod[-1],args=(1,1,1,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        KSmn[i,4] = Ks11*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)

        K20,errin = sp.integrate.quad(Kmn_outer_integ,kmod[0],kmod[-1],args=(2,0,0,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        KMN[i,2,0]= K20*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)
        Ks20,errin= sp.integrate.quad(Kmn_outer_integ,kmod[0],kmod[-1],args=(2,0,1,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        KSmn[i,5] = Ks20*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)
        
        K30,errin = sp.integrate.quad(Kmn_outer_integ,kmod[0],kmod[-1],args=(3,0,0,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        KMN[i,2,1]= K30*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)
        Ks30,errin= sp.integrate.quad(Kmn_outer_integ,kmod[0],kmod[-1],args=(3,0,1,ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        KSmn[i,6] = Ks30*ks[i]*ks[i]*ks[i]/(4.0*math.pi*math.pi)
    return KMN ,KSmn

#--------------------------------------------------------------------------
def pk_integ(r,kval,Pl_spline):
    P_lin = sp.interpolate.splev(np.log10(r*kval),Pl_spline, der=0)
    return r*r*10.0**(P_lin)*10.0**(P_lin)

def KNORm(kmod,Pl_spline,ks): 
    xsa=np.zeros(len(ks))
    for i in range(len(ks)):
        integ,errin = sp.integrate.quad(pk_integ,kmod[0],kmod[len(kmod)-1],args=(ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3)
        xsa[i]=ks[i]*ks[i]*ks[i]*integ/(2.0*math.pi*math.pi)
    return xsa

def SIGMA3_inner_integ(x,r):
    y     = np.sqrt(1.0+r*r-2.0*r*x)
    numer = (2.0/7.0)*(x*x - 1.0)*(3.0*x*x - 4.0*r*x + 2.0*r*r - 1.0)
    return numer/(3.0*y*y) + 8.0/63.0
    
def SIGMA3_outter_integ(r,k,Pl_spline):  
    integ,errin= sp.integrate.quad(SIGMA3_inner_integ,-1.0,1.0,epsabs=0.0,epsrel=1.0e-4,args=(r))
    Pl   = sp.interpolate.splev(np.log10(k*r),Pl_spline, der=0)
    return integ*r*r*10.0**(Pl)

def SIGMA3(kmod,Pl_spline,ks):
    SIGMAs=np.zeros(len(ks))
    for i in range(len(ks)):
        integ,errin= sp.integrate.quad(SIGMA3_outter_integ,kmod[0],kmod[-1],args=(ks[i],Pl_spline),epsabs=0.0,epsrel=1.0e-3) 
        SIGMAs[i]=integ*(105.0*ks[i]*ks[i]*ks[i])/(64.0*math.pi*math.pi)
    return SIGMAs

#   the parameter estimation main code:------------------------------------
def PSlin_Integrations( KmodLim,hub,ombh2,omch2,ns, As,out_dir ):  
    #Linear spectra
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=hub,ombh2=ombh2, omch2=omch2)
    pars.set_dark_energy()  
    pars.InitPower.set_params(As=As,ns=ns)
    pars.set_matter_power(redshifts=[0.0], kmax=100.0)
    pars.NonLinear = model.NonLinear_none
    results = camb.get_results(pars)
    kh, z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=100.0, npoints = 692)
    pk=pk[0]; ks=kh[kh<=KmodLim]
    print( 'sigma8_fid_integ = ',results.get_sigma8(),'\n')
    Pl_spline = sp.interpolate.splrep(np.log10(kh),np.log10(pk ), s=0)
    # integrations:
    print( '\n Integration-Imn \n')
    Imn=I_mn(kh,Pl_spline,ks)
    print( '\n Integration-Jmn \n')
    Jmn=J_mn(kh,Pl_spline,ks)
    print( '\n Integration-Kmn \n')
    Kmn,Ksmn=K_mn(kh,Pl_spline,ks)
    print( '\n Integration-Sigma3 \n')
    SIGMA_3=SIGMA3(kh,Pl_spline,ks)
    print( '\n Integration-norm \n')
    Intg_norm=KNORm(kh,Pl_spline,ks)
    # outputs:
    I00=Imn[:,0,0]   ; I01=Imn[:,0,1]   ; I02=Imn[:,0,2]   ; I03=Imn[:,0,3]
    I10=Imn[:,1,0]   ; I11=Imn[:,1,1]   ; I12=Imn[:,1,2]   ; I13=Imn[:,1,3]
    I20=Imn[:,2,0]   ; I21=Imn[:,2,1]   ; I22=Imn[:,2,2]   ; I23=Imn[:,2,3]
    I30=Imn[:,3,0]   ; I31=Imn[:,3,1]   ; I32=Imn[:,3,2]   ; I33=Imn[:,3,3]
    J00=Jmn[:,0,0]   ; J01=Jmn[:,0,1]   ; J02=Jmn[:,0,2]   ; J10=Jmn[:,1,0] 
    J11=Jmn[:,1,1]   ; J20=Jmn[:,1,2]       
    K00=Kmn[:,0,0]   ; K01=Kmn[:,0,1]   
    K10=Kmn[:,1,0]   ; K11=Kmn[:,1,1]
    K20=Kmn[:,2,0]   ; K30=Kmn[:,2,1]       
    Ks00=Ksmn[:,0]   ; Ks01=Ksmn[:,1]   ; Ks02=Ksmn[:,2]   ; Ks10=Ksmn[:,3]        
    Ks11=Ksmn[:,4]   ; Ks20=Ksmn[:,5]   ; Ks30=Ksmn[:,6]
    outfile    = open(out_dir+'INTEG_CH.npy', 'w')
    for i in range(len(ks)):
       outfile.write("%12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf  %12.6lf\n"% (kh[i],pk[i],   I00[i],  I01[i],  I02[i],  I03[i],  I10[i],  I11[i],  I12[i],  I13[i],  I20[i],  I21[i],  I22[i],  I23[i],  I30[i],  I31[i],  I32[i],  I33[i],  J00[i],  J01[i],  J02[i],  J10[i],  J11[i],  J20[i],  K00[i],  Ks00[i], K01[i],  Ks01[i], Ks02[i], K10[i],  Ks10[i], K11[i],  Ks11[i], K20[i],  Ks20[i], K30[i],  Ks30[i], SIGMA_3[i],Intg_norm[i]))
    outfile.close()
    PT=np.loadtxt(out_dir+'INTEG_CH.npy')
    PT=PT.T
    PT[26,0:] -= PT[38,0:];PT[27,0:] -= 4.0/9.0*PT[38,0:];PT[28,0:] -= 2.0/3.0*PT[38,0:]
    np.save(out_dir+'INTEG_CH',np.array([PT,results.get_sigma8()],dtype=object),allow_pickle=True)    
    return PT,results.get_sigma8()





#==============================================================================
#                               Model power
# 3. calculate model mom-PS:---------------------------------------------------
def PijFun(parm,PT,MU,pstyp):
    if(pstyp=='auto'):
        f, b1,b2,bs,b3nl, sigv1_squre,sigv2_squre,sigvf_squre,sigma8rat = parm
    if(pstyp=='crs'):  
        f, b1,b2,bs,b3nl, b1v,b2v,bsv,b3nlv, sigv1_squre,sigv2_squre,sigvf_squre,  sigvv1_squre,sigvv2_squre,sigvvf_squre,  sigma8rat = parm  
    k  =PT[0,0:] ; Pl =PT[1,0:]   
    I00=PT[2,0:]  ;I01=PT[3,0:]  ;I02=PT[4,0:]  ;I03=PT[5,0:]  ;I10=PT[6,0:]  ;I11=PT[7,0:]  ;I12=PT[8,0:]  ;I13=PT[9,0:] ;
    I20=PT[10,0:] ;I21=PT[11,0:] ;I22=PT[12,0:] ;I23=PT[13,0:] ;I30=PT[14,0:] ;I31=PT[15,0:] ;I32=PT[16,0:] ;I33=PT[17,0:] ;
    J00=PT[18,0:] ;J01=PT[19,0:];J02=PT[20,0:]; J10=PT[21,0:]; J11=PT[22,0:];  J20=PT[23,0:]; K00=PT[24,0:]; Ks00=PT[25,0:];
    K01=PT[26,0:]; Ks01=PT[27,0:];Ks02=PT[28,0:];K10=PT[29,0:];Ks10=PT[30,0:]; K11=PT[31,0:];
    Ks11=PT[32,0:];K20=PT[33,0:];Ks20=PT[34,0:];K30=PT[35,0:];Ks30=PT[36,0:];  sig3_squre=PT[37,0:];#Intg_norm=PT[38,0:];   
    sig4_squre=1.0/(24.0*math.pi*math.pi)*(sigma8rat**2*sp.integrate.simps( (I23+2./3.*I32+I33/5.) /k**2,k)) 
    
    #sigvL_squre=1.0/(6.0*math.pi**2)*sp.integrate.simps(sigma8rat*k*k*(Pl + sigma8rat*(2.0*I11 + 6.*k*k*Pl*J11)), k) 
    sigvL_squre=sigv1_squre +0-0     ;    sigv1_squre =0; sigv2_squre=0;    sigvv1_squre=0; sigvv2_squre=0;  
    
    P_00=np.zeros((len(k),len(MU)));P_01=np.zeros((len(k),len(MU)));P_02=np.zeros((len(k),len(MU)))
    P_03=np.zeros((len(k),len(MU)));P_04=np.zeros((len(k),len(MU)));P_11=np.zeros((len(k),len(MU)))
    P_12=np.zeros((len(k),len(MU)));P_13=np.zeros((len(k),len(MU)));P_22=np.zeros((len(k),len(MU)))
    for ss in range(len(MU)):
        mu=MU[ss]
        if(pstyp=='auto'):# Note: sig3_squre is not sigv3_squre, and sig4_squre is not sigv4_squre, they are different parameters, do not confused with them.
            P00=(b1+b1)*sigma8rat**2*(b2*K00+bs*Ks00)+b1*b1*sigma8rat*(Pl+2.*sigma8rat*(I00+3.*k*k*Pl*J00))+1./2.*sigma8rat**2*(b2*b2*K01+bs*bs*Ks01+2.*b2*bs*Ks02+4.*b3nl*sig3_squre*Pl)
            P01=f*b1*sigma8rat*(Pl+2.*sigma8rat*(I01+b1*I10+3.*k*k*Pl*(J01+b1*J10))-b2*sigma8rat*K11-bs*sigma8rat*Ks11)-f*sigma8rat**2*(b2*K10+bs*Ks10+b3nl*sig3_squre*Pl)
            P02=f*f*b1*sigma8rat**2*(I02+mu**2*I20+2.*k*k*Pl*(J02+mu**2*J20))-f*f*k*k*(sigvL_squre+sigv1_squre/f**4 )* P00+f*f*sigma8rat**2*(b2*(K20+mu**2*K30)+bs*(Ks20+mu**2*Ks30))
            P03=-f*f*k*k*(sigvL_squre+sigv2_squre/f**4)*P01
            P04=-1./2.*f*f*f*f*b1*k*k*(sigvL_squre+sigv1_squre/f**4)*sigma8rat**2*(I02+mu**2*I20+2.*k*k*Pl*(J02+mu**2*J20)) +1./4.*f*f*f*f*b1*b1*k**4*P00*((sigvL_squre+sigv1_squre/f**4)**2+sig4_squre/3. )
            P11=f*f*sigma8rat*(mu**2*(Pl+sigma8rat*(2.*I11+(b1+b1+b1+b1)*I22+b1*b1*I13+6.*k**2*Pl*(J11+(b1+b1)*J10)))  +b1*b1*sigma8rat*I31)
            P12=f*f*f*sigma8rat**2*(I12+mu**2*I21-b1*(I03+mu**2*I30)+2.*k**2*Pl*(J02+mu**2*J20))   -f*f*k*k*(sigvL_squre+sigv1_squre/f**4)* P01+2.*f*f*f*k*k*sigma8rat**2*(I01+I10+3.*k*k*Pl*(J01+J10))   *sigvf_squre 
            P13=-f*f*k*k*f*f*sigma8rat*((sigvL_squre+sigv2_squre/f**4)*mu**2*(Pl+sigma8rat*(2.*I11+(b1+b1+b1+b1)*I22 +6.*k*k*Pl*(J11+(b1+b1)*J10))) +(sigvL_squre+sigv1_squre/f**4)*b1*b1*sigma8rat*(mu**2*I13+I31) )
            P22=1./4.*f*f*f*f*sigma8rat**2*(I23+2.*mu**2*I32+mu**4*I33)+f*f*f*f*k**4*(sigvL_squre+sigv1_squre/f**4 )**2 *P00-f*f*k*k*(sigvL_squre+sigv1_squre/f**4 )*(2.*P02-f*f*sigma8rat**2*(b2*(K20+mu**2*K30)+bs*(Ks20+mu**2*Ks30)))          
        if(pstyp=='crs'): 
            P00=(b1+b1)*sigma8rat**2*(b2*K00+bs*Ks00)+b1*b1*sigma8rat*(Pl+2.*sigma8rat*(I00+3.*k*k*Pl*J00))+1./2.*sigma8rat**2*(b2*b2*K01+bs*bs*Ks01+2.*b2*bs*Ks02+4.*b3nl*sig3_squre*Pl)
            P01=f*b1*sigma8rat*(Pl+2.*sigma8rat*(I01+b1*I10+3.*k*k*Pl*(J01+b1*J10))-b2*sigma8rat*K11-bs*sigma8rat*Ks11)-f*sigma8rat**2*(b2*K10+bs*Ks10+b3nl*sig3_squre*Pl)
            P02=f*f*b1*sigma8rat**2*(I02+mu**2*I20+2.*k*k*Pl*(J02+mu**2*J20))-f*f*k*k*(sigvL_squre+sigv1_squre/f**4 )* P00+f*f*sigma8rat**2*(b2*(K20+mu**2*K30)+bs*(Ks20+mu**2*Ks30))
            P03=-f*f*k*k*(sigvL_squre+sigv2_squre/f**4)*P01
            P04=-1./2.*f*f*f*f*b1*k*k*(sigvL_squre+sigv1_squre/f**4)*sigma8rat**2*(I02+mu**2*I20+2.*k*k*Pl*(J02+mu**2*J20)) +1./4.*f*f*f*f*b1*b1*k**4*P00*((sigvL_squre+sigv1_squre/f**4)**2+sig4_squre/3. )
            P11=f*f*sigma8rat*(mu**2*(Pl+sigma8rat*(2.*I11+(b1+b1+b1+b1)*I22+b1*b1*I13+6.*k**2*Pl*(J11+(b1+b1)*J10)))  +b1*b1*sigma8rat*I31)
            P12=f*f*f*sigma8rat**2*(I12+mu**2*I21-b1*(I03+mu**2*I30)+2.*k**2*Pl*(J02+mu**2*J20))   -f*f*k*k*(sigvL_squre+sigv1_squre/f**4)* P01+2.*f*f*f*k*k*sigma8rat**2*(I01+I10+3.*k*k*Pl*(J01+J10))   *sigvf_squre 
            P13=-f*f*k*k*f*f*sigma8rat*((sigvL_squre+sigv2_squre/f**4)*mu**2*(Pl+sigma8rat*(2.*I11+(b1+b1+b1+b1)*I22 +6.*k*k*Pl*(J11+(b1+b1)*J10))) +(sigvL_squre+sigv1_squre/f**4)*b1*b1*sigma8rat*(mu**2*I13+I31) )
            P22=1./4.*f*f*f*f*sigma8rat**2*(I23+2.*mu**2*I32+mu**4*I33)+f*f*f*f*k**4*(sigvL_squre+sigv1_squre/f**4 )**2 *P00-f*f*k*k*(sigvL_squre+sigv1_squre/f**4 )*(2.*P02-f*f*sigma8rat**2*(b2*(K20+mu**2*K30)+bs*(Ks20+mu**2*Ks30)))          
        P_00[:,ss]=P00;   P_01[:,ss]=P01;   P_02[:,ss]=P02
        P_03[:,ss]=P03;   P_04[:,ss]=P04;   P_11[:,ss]=P11
        P_12[:,ss]=P12;   P_13[:,ss]=P13;  P_22[:,ss]=P22 
    return   P_00,P_01,P_02,P_03,P_04,P_11,P_12,P_13,P_22  
def Pk_MOD(params,sigma8_fid,PT,ps_type, Growthz ):
    
    fsigma8,  b1sigma8,b2sigma8,b3nlsigma8,   sigmav_squre  = params 
    
    sigmav_squre=np.abs(sigmav_squre)
    sigmavf_squre=np.abs(sigmav_squre)
    b1vsigma8,b2vsigma8 = b1sigma8, b2sigma8
    b3nlvsigma8 =  b3nlsigma8
    sigma_v1_squre ,sigma_v2_squre  = sigmav_squre , sigmav_squre
    sigmav_v1_squre,sigmav_v2_squre = sigmav_squre , sigmav_squre  
    sigmav_vf_squre=sigmavf_squre
    
    # normalize parms: 
    f = fsigma8/sigma8_fid
    sigma8rat = Growthz**2 
    # momentum field: 
    b1v       = b1vsigma8/sigma8_fid
    b2v       = b2vsigma8/sigma8_fid 
    bsv       = -4.0/7.0*(b1v-1.0)
    b3nlv     = b3nlvsigma8/sigma8_fid # 32.0/315.0*(b1v-1.0)    
    # density field:  
    b1        = b1sigma8/sigma8_fid 
    b2        = b2sigma8/sigma8_fid
    bs        = -4.0/7.0*(b1-1.0)
    b3nl      = b3nlsigma8/sigma8_fid #32.0/315.0*(b1-1.0)
    # model PS:
    P_model=[];PSTYP=[]
    mu = np.linspace(0.0, 1.0, 100)    
    if ('den' in ps_type):
        P_00,P_01,P_02,P_03,P_04,P_11,P_12,P_13,P_22 =PijFun( [f, b1,b2,bs,b3nl, sigma_v1_squre,sigma_v2_squre, sigmavf_squre,sigma8rat ],PT,mu, 'auto')
        P_den  = P_00 + mu**2*(2.0*P_01 + P_02 + P_11) + mu**4*(P_03 + P_04 + P_12 + P_13 + 1.0/4.0*P_22)
        P_0den = (2.0*0 + 1.0)  * sp.integrate.simps(P_den,mu,axis=1)
        P_2den = (2.0*2 + 1.0)  * sp.integrate.simps( ((3.*(mu**2)-1.)/2.)*P_den,mu,axis=1)
        P_4den = (2.0*4 + 1.0)  * sp.integrate.simps( ((35.*(mu**4)-30.*(mu**2)+3.)/8. )*P_den,mu,axis=1)  
        P_model.append([P_0den,P_2den,P_4den]) ; PSTYP.append( 'den' )
    if ('mom' in ps_type):
        P_00,P_01,P_02,P_03,P_04,P_11,P_12,P_13,P_22 =PijFun( [f, b1v,b2v,bsv,b3nlv, sigmav_v1_squre,sigmav_v2_squre, sigmav_vf_squre,sigma8rat ],PT,mu, 'auto')
        P_mom  = ((1.0e4)/(1.0**2)/(PT[0,0:]*PT[0,0:])*(P_11 + mu**2*(2.0*P_12 + 3.0*P_13 + P_22)).T).T
        P_0mom = (2.0*0 + 1.0)  * sp.integrate.simps(P_mom,mu,axis=1)
        P_2mom = (2.0*2 + 1.0)  * sp.integrate.simps( ((3.*(mu**2)-1.)/2.)*P_mom,mu,axis=1) 
        P_4mom = (2.0*4 + 1.0)  * sp.integrate.simps( ((35.*(mu**4)-30.*(mu**2)+3.)/8. )*P_mom,mu,axis=1)  
        P_model.append([P_0mom,P_2mom,P_4mom])  ;  PSTYP.append( 'mom' )
    if ('crs' in ps_type):
        P_00,P_01,P_02,P_03,P_04,P_11,P_12,P_13,P_22 =PijFun( [f, b1,b2,bs,b3nl, b1v,b2v,bsv,b3nlv, sigma_v1_squre,sigma_v2_squre,  sigmavf_squre,sigmav_v1_squre,sigmav_v2_squre, sigmav_vf_squre,sigma8rat],PT,mu, 'crs')
        P_crs  = (1.0e2/PT[0,0:]*((mu*(P_01 + P_02 + P_11 + mu**2*(3.0/2.0*P_03 + 2.0*P_04 + 2.0*P_12 + 3.0*P_13  + 1.0/2.0*P_22 ))).T)).T
        P_1crs = (2.0*1 + 1.0)  * sp.integrate.simps(mu*P_crs,mu,axis=1)
        P_3crs = (2.0*3 + 1.0)  * sp.integrate.simps( (1./2.*(5.*mu**3-3.*mu) )*P_crs,mu,axis=1)   
        P_model.append([P_1crs,P_3crs])  ; PSTYP.append( 'crs' )
    return P_model,PSTYP
def Pkmod_Fun(params,OnlyL0,Sig8_fid, k_obs, kmod_conv,WF_Kobs_Kmod,WF_rand,PT,ps_type, Growthz, WFCon_type):
    Pk_mods , psty= Pk_MOD(params,Sig8_fid,PT,ps_type, Growthz )
    Pk_mod   = []   ;  Pk_modunc   = []   ; kmod_uncon=[]
    for i in range(len(Pk_mods)):
        Pk_modss=np.concatenate(Pk_mods[i] )
        if(WFCon_type=='Ross'):
            Nic=3
            Nks=len(k_obs)
            if(len(Pk_mods)==2):Nks=Nks//2
            if(i==0):
                P_spline_0 = sp.interpolate.splrep(PT[0,:],Pk_modss[0:(len(Pk_modss)//Nic)],s=0)
                Pk_mod0    = sp.interpolate.splev(kmod_conv , P_spline_0)
                Pk_mod     = Pk_CONV_ross(k_obs[:Nks] ,WF_Kobs_Kmod,Pk_mod0,WF_rand[0],WF_rand[1]) 
            if(i==1): 
                P_spline_0 = sp.interpolate.splrep(PT[0,:],Pk_modss[0:(len(Pk_modss)//Nic)],s=0)
                Pk_modp0   = sp.interpolate.splev(kmod_conv , P_spline_0)
                Pk_modp    = Pk_CONV_ross(k_obs[:Nks] ,WF_Kobs_Kmod,Pk_modp0,WF_rand[0],WF_rand[1])
                Pk_mod     = np.concatenate((Pk_mod,Pk_modp))
                Pk_mod0    = np.concatenate((Pk_mod0,Pk_modp0))
                kmod_conv  = np.concatenate((kmod_conv,kmod_conv))
            if(len(Pk_mods)==1):  
                return Pk_mod, k_obs, Pk_mod0, kmod_conv
            else:
                if(i==1):
                    return Pk_mod, k_obs, Pk_mod0, kmod_conv
        if(WFCon_type!='Ross')and(OnlyL0):    
            Nic=3
            Nks=len(k_obs)
            if(i==0):
                P_spline_0 = sp.interpolate.splrep(PT[0,:],Pk_modss[0:(len(Pk_modss)//Nic)],s=0)
                Pk_mod0    = sp.interpolate.splev(kmod_conv[i] , P_spline_0) 
                Pk_mod     = Pk_CONV_multi( WF_Kobs_Kmod[i],[Pk_mod0,Pk_mod0,Pk_mod0],  psty[i])[0]
            if(i==1): 
                P_spline_0 = sp.interpolate.splrep(PT[0,:],Pk_modss[0:(len(Pk_modss)//Nic)],s=0)
                Pk_modp0   = sp.interpolate.splev(kmod_conv[i] , P_spline_0) 
                Pk_modp    = Pk_CONV_multi( WF_Kobs_Kmod[i],[Pk_modp0,Pk_modp0,Pk_modp0],  psty[i])[0]
                Pk_mod     = np.concatenate((Pk_mod,Pk_modp))
                Pk_mod0    = np.concatenate((Pk_mod0,Pk_modp0))
                kmod_con   = np.concatenate((kmod_conv[i],kmod_conv[i] ))
            if(len(Pk_mods)==1):  
                return Pk_mod, k_obs, Pk_mod0, kmod_conv[0]
            else:
                if(i==1):
                    return Pk_mod, k_obs, Pk_mod0, kmod_con 
        if(WFCon_type!='Ross')and(not OnlyL0): 
            if(psty[i]=='den')or(psty[i]=='mom'):
                Nic=3
            if(psty[i]=='crs'):
                Nic=2
            P_spline_0     = sp.interpolate.splrep(PT[0,:],Pk_modss[0                     :(len(Pk_modss)//Nic)  ],s=0)
            Pk_mod0        = sp.interpolate.splev(kmod_conv[i], P_spline_0)
            P_spline_2     = sp.interpolate.splrep(PT[0,:],Pk_modss[(len(Pk_modss)//Nic)  :(len(Pk_modss)//Nic*2)],s=0)
            Pk_mod2        = sp.interpolate.splev(kmod_conv[i], P_spline_2)            
            if(psty[i]=='den')or(psty[i]=='mom'):
                P_spline_4 = sp.interpolate.splrep(PT[0,:],Pk_modss[(len(Pk_modss)//Nic*2):(len(Pk_modss)//Nic*3)],s=0)
                Pk_mod4    = sp.interpolate.splev(kmod_conv[i], P_spline_4) 
            if(psty[i]=='den'):
                Pk_modm0 ,Pk_modm2 ,Pk_modm4    = Pk_CONV_multi( WF_Kobs_Kmod[0],[Pk_mod0,Pk_mod2,Pk_mod4],  psty[i])
                Pk_mod.append(Pk_modm0)  ; Pk_mod.append(Pk_modm2)        ; Pk_mod.append(Pk_modm4) 
                Pk_modunc.append(Pk_mod0); Pk_modunc.append(Pk_mod2); Pk_modunc.append(Pk_mod4)
                kmod_uncon.append(kmod_conv[0]);kmod_uncon.append(kmod_conv[0]);kmod_uncon.append(kmod_conv[0])
            if(psty[i]=='mom'):
                Pk_modm0 ,Pk_modm2 ,Pk_modm4    = Pk_CONV_multi( WF_Kobs_Kmod[1],[Pk_mod0,Pk_mod2,Pk_mod4],  psty[i])
                Pk_mod.append(Pk_modm0)  ; Pk_mod.append(Pk_modm2)        ; Pk_mod.append(Pk_modm4) 
                Pk_modunc.append(Pk_mod0); Pk_modunc.append(Pk_mod2); Pk_modunc.append(Pk_mod4)
                kmod_uncon.append(kmod_conv[1]);kmod_uncon.append(kmod_conv[1]);kmod_uncon.append(kmod_conv[1])
            if(psty[i]=='crs'):
                Pk_modm1 ,Pk_modm3      = Pk_CONV_multi( WF_Kobs_Kmod[2],[Pk_mod0,Pk_mod2 ],  psty[i])
                Pk_mod.append(Pk_modm1)  ; Pk_mod.append(Pk_modm3) 
                Pk_modunc.append(Pk_mod0); Pk_modunc.append(Pk_mod2)
                kmod_uncon.append(kmod_conv[2]);kmod_uncon.append(kmod_conv[2])
    if(WFCon_type!='Ross')and(not OnlyL0):
        return   np.concatenate(Pk_mod),   k_obs,     np.concatenate(Pk_modunc),  np.concatenate(kmod_uncon)    
                
                
                
                
   




#==============================================================================
#                               Fitting Parameters
# 6 the parameter estimation main code:----------------------------------------
def PSmulti_pickup(ind_fit,ps_type,Pk_obs,Pk_modc):
    Ps_type = ps_type.split(' ')
    mt=[]
    for i in range(len(Ps_type) ):
        ps_type0=Ps_type[i].split('-')
        tmp=[]
        for j in range(len(ps_type0[1])):
            tmp.append( int(ps_type0[1][j]) )
        mt.append(np.array(tmp))  
    if('den' in ps_type)and('mom' not in ps_type)and('crs' not in ps_type): nn=3               
    if('mom' in ps_type)and('den' not in ps_type)and('crs' not in ps_type): nn=3
    if('crs' in ps_type)and('den' not in ps_type)and('den' not in ps_type): nn=2    
    if('den' in ps_type)and('mom' in ps_type)and('crs' not in ps_type): nn=6
    if('den' in ps_type)and('mom' in ps_type)and('crs'     in ps_type): nn=8
    if('den' in ps_type)and('mom' not in ps_type)and('crs' in ps_type): nn=5
    if('den' not in ps_type)and('mom' in ps_type)and('crs' in ps_type): nn=5                    
    Nf=len(ind_fit)//nn
    Pk_obs=Pk_obs[ind_fit]
    Pk_modc=Pk_modc[ind_fit]
    Pko=[]; Pkm=[]  ; I=0;J=1
    for i in range(len(mt)):  
      if(i==0): 
        if(0 in mt[i])or(1 in mt[i]):
            Pko.append(Pk_obs[    I*Nf:J*Nf]) 
            Pkm.append(Pk_modc[   I*Nf:J*Nf])
        I=I+1;J=J+1
        if(2 in mt[i])or(3 in mt[i]):
            Pko.append(Pk_obs[    I*Nf:J*Nf])
            Pkm.append(Pk_modc[   I*Nf:J*Nf])
        I=I+1;J=J+1
        if(4 in mt[i]):     
            Pko.append(Pk_obs[    I*Nf:J*Nf])
            Pkm.append(Pk_modc[   I*Nf:J*Nf]) 
        I=I+1;J=J+1     
      if(i==1):    
        if(0 in mt[i])or(1 in mt[i]):
            Pko.append(Pk_obs[    I*Nf:J*Nf]) 
            Pkm.append(Pk_modc[   I*Nf:J*Nf])
        I=I+1;J=J+1
        if(2 in mt[i])or(3 in mt[i]):
            Pko.append(Pk_obs[    I*Nf:J*Nf])
            Pkm.append(Pk_modc[   I*Nf:J*Nf])
        I=I+1;J=J+1
        if(4 in mt[i]):     
            Pko.append(Pk_obs[    I*Nf:J*Nf])
            Pkm.append(Pk_modc[   I*Nf:J*Nf])
        I=I+1;J=J+1                 
      if(i==2):     
        if(1 in mt[i]):
            Pko.append(Pk_obs[    I*Nf:J*Nf]) 
            Pkm.append(Pk_modc[   I*Nf:J*Nf])
        I=I+1;J=J+1
        if(3 in mt[i]):
            Pko.append(Pk_obs[    I*Nf:J*Nf])
            Pkm.append(Pk_modc[   I*Nf:J*Nf])
        I=I+1;J=J+1
    Pk_obs=np.concatenate(Pko)
    Pk_modc=np.concatenate(Pkm)
    return Pk_obs,Pk_modc
  
def CHI2(params,    ind_fit,Cov_inv,OnlyL0,Sig8_fid, k_obs, Pk_obs,kmod_conv, WF_Kobs_Kmod,WF_rand,PT,ps_type, Growthz, WFCon_type ):
    Pk_modc, kmodc, Pk_mod , kmod =Pkmod_Fun(params,OnlyL0,Sig8_fid, k_obs, kmod_conv, WF_Kobs_Kmod,WF_rand,PT,ps_type, Growthz, WFCon_type )
    if(OnlyL0):
        VC=np.dot((Pk_obs[ind_fit]-Pk_modc[ind_fit]),Cov_inv)
        chi_squared=np.dot(VC,(Pk_obs[ind_fit]-Pk_modc[ind_fit]))
    if(not OnlyL0):
        Pk_obs,Pk_modc=PSmulti_pickup(ind_fit,ps_type,Pk_obs,Pk_modc)
        VC=np.dot((Pk_obs -Pk_modc ),Cov_inv)
        chi_squared=np.dot(VC,(Pk_obs -Pk_modc ))           
    return chi_squared
# 5. the posterior and prior of MCMC:--------------------------------------
def lnprior(params):
    fsigma8,  b1sigma8,b2sigma8,b3nlsigma8,   sigmav_squre    = params   
    if(0.  <fsigma8<1.5): pdf=1.
    else:return -np.inf
    if(0.  <b1sigma8<3.):pdf=1.
    else:return -np.inf
    if(-3.  <b2sigma8<3.):pdf=1.
    else:return -np.inf
    if(-3.  <b3nlsigma8<3.):pdf=1.
    else:return -np.inf
    if( 0.<sigmav_squre<250.):pdf=1.
    else:return -np.inf
    return 0.0
# 5.1: the post of MCMC:
def lnpost_chi2(params    ,   Cp,Nmock,ind_fit,Cov_inv,OnlyL0,Sig8_fid, k_obs, Pk_obs,kmod_conv, WF_Kobs_Kmod,WF_rand,PT,ps_type, Growthz, WFCon_type  ):
    prior = lnprior(params )
    if not np.isfinite(prior):
        return -np.inf
    like = -0.5* CHI2(params ,ind_fit,Cov_inv,OnlyL0,Sig8_fid, k_obs, Pk_obs,kmod_conv, WF_Kobs_Kmod,WF_rand,PT,ps_type, Growthz, WFCon_type )
    return prior + like
 
 
        
        
        
        
        
        
        
        
        
        
        
        

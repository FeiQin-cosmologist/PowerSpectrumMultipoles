import numpy as np
import math
from scipy.special import sph_harm
import scipy as sp
from scipy.interpolate import splev, splrep
from scipy import interpolate
from scipy import integrate
# Speed of light in km/s
LightSpeed = 299792.458
# Calculates H(z)/H0
def Ez(redshift, omega_m, omega_lambda, omega_rad, w0, wa, ap):
    fz = ((1.0+redshift)**(3*(1.0+w0+wa*ap)))*math.exp(-3*wa*(redshift/(1.0+redshift)))
    omega_k = 1.0-omega_m-omega_lambda-omega_rad
    return math.sqrt(omega_rad*(1.0+redshift)**4+omega_m*(1.0+redshift)**3+omega_k*(1.0+redshift)**2+omega_lambda*fz)
# The Comoving Distance Integrand
def DistDcIntegrand(redshift, omega_m, omega_lambda, omega_rad, w0, wa, ap):
    return 1.0/Ez(redshift, omega_m, omega_lambda, omega_rad, w0, wa, ap)
# The Comoving Distance in Mpc
def DistDc(redshift, omega_m, omega_lambda, omega_rad, Hubble_Constant, w0, wa, ap):
    return (LightSpeed/Hubble_Constant)*integrate.quad(DistDcIntegrand, 0.0, redshift, args=(omega_m, omega_lambda, omega_rad, w0, wa, ap))[0]
def spline_Dis2Rsf(OmegaM,OmegaA,Hub,nbin=10000,redmax=2.5):
    dist = np.empty(nbin);red = np.empty(nbin)
    for j in range(nbin):
        red[j] = j*redmax/nbin
        dist[j] = DistDc(red[j],OmegaM,OmegaA, 0.0,Hub,-1.0, 0.0, 0.0)
    rsf_spline = sp.interpolate.splrep(dist,red,  s=0)
    return rsf_spline
def spline_Rsf2Dis(OmegaM,OmegaA,Hub,nbin=10000,redmax=2.5):
    dist = np.empty(nbin);red = np.empty(nbin)
    for j in range(nbin):
        red[j] = j*redmax/nbin
        dist[j] = DistDc(red[j],OmegaM,OmegaA, 0.0,Hub,-1.0, 0.0, 0.0)
    dist_spline = sp.interpolate.splrep(red, dist, s=0)
    return dist_spline
def DisRsfConvert(xdt,types,OmegaM,OmegaA,Hub,nbin=10000,redmax=2.5):
    if(types=='z2d'):
        spl_fun=spline_Rsf2Dis(OmegaM,OmegaA,Hub,nbin,redmax)
        Distc        = splev(xdt, spl_fun)
        return Distc
    if(types=='d2z'):
        spl_fun=spline_Dis2Rsf(OmegaM,OmegaA,Hub,nbin,redmax)
        RSFs        = splev(xdt, spl_fun)
        return RSFs
def Sky2Cat(ra,dec,rsft,OmegaM , OmegaA ,Hub,nbin=1000,redmax=2.):
    disz= DisRsfConvert(rsft,'z2d',OmegaM,OmegaA,Hub,nbin,redmax)
    X = disz*np.cos(dec/180.*math.pi)*np.cos(ra/180.*math.pi)
    Y = disz*np.cos(dec/180.*math.pi)*np.sin(ra/180.*math.pi)
    Z = disz*np.sin(dec/180.*math.pi)    
    return X,Y,Z 
# The Linear Growth Factor Integrand assuming GR
def GrowthFactorGRIntegrand(scale_factor, omega_m, omega_lambda, omega_rad, w0, wa, ap):
    redshift = (1.0/scale_factor)-1.0
    return 1.0/(scale_factor*Ez(redshift, omega_m, omega_lambda, omega_rad, w0, wa, ap))**3
# The Linear Growth Factor assuming GR
def GrowthFactorGR(redshift, omega_m, omega_lambda, omega_rad, Hubble_Constant, w0, wa, ap):
    prefac = Ez(redshift, omega_m, omega_lambda, omega_rad, w0, wa, ap)
    scale_factor = 1.0/(1.0+redshift)
    return prefac*integrate.quad(GrowthFactorGRIntegrand, 0.0, scale_factor, args=(omega_m, omega_lambda, omega_rad, w0, wa, ap))[0]
# Omega_M at a given redshift
def Omega_m_z(redshift, omega_m, omega_lambda, omega_rad, w0, wa, ap):
    return (omega_m*(1.0+redshift)**3)/(Ez(redshift, omega_m, omega_lambda, omega_rad, w0, wa, ap)**2)
# The Linear Growth Factor Integrand for an arbitrary value of gamma
def GrowthFactorGammaIntegrand(scale_factor, gamma, omega_m, omega_lambda, omega_rad, w0, wa, ap):
    redshift = (1.0/scale_factor)-1.0
    return (Omega_m_z(redshift, omega_m, omega_lambda, omega_rad, w0, wa, ap)**gamma)/scale_factor
# The Linear Growth Factor for an arbitrary value of gamma
def GrowthFactorGamma(gamma, redshift_low, redshift_high, omega_m, omega_lambda, omega_rad, Hubble_Constant, w0, wa, ap):
    scale_factor_low = 1.0/(1.0+redshift_low)
    scale_factor_high = 1.0/(1.0+redshift_high)
    return math.exp(integrate.quad(GrowthFactorGammaIntegrand, scale_factor_low, scale_factor_high, args=(gamma, omega_m, omega_lambda, omega_rad, w0, wa, ap))[0])
def Fsigma8_Fun(gamma,RsFeff,Sig8_fid, OmegaM, OmegaA,  Hub):
    sig8_fac  = GrowthFactorGR(1000., OmegaM, OmegaA, 0.0, Hub, -1.0, 0.0, 0.0)/GrowthFactorGR(0.0, OmegaM, OmegaA, 0.0, Hub, -1.0, 0.0, 0.0)
    gamma_fac = GrowthFactorGamma(gamma, 1000.,  RsFeff ,  OmegaM, OmegaA, 0.0, Hub, -1.0, 0.0, 0.0)
    OmegaMz   = Omega_m_z(RsFeff, OmegaM, OmegaA, 0.0, -1.0, 0.0, 0.0)
    fsig8     = OmegaMz**gamma*Sig8_fid*sig8_fac*gamma_fac 
    return fsig8



 

def getkspec(nx,ny,nz,lx,ly,lz):
    kx = 2.*np.pi*np.fft.fftfreq(nx,d=lx/nx)
    ky = 2.*np.pi*np.fft.fftfreq(ny,d=ly/ny)
    kz = 2.*np.pi*np.fft.fftfreq(nz,d=lz/nz)
    indep = np.full((nx,ny,nz),True,dtype=bool)
    indep[0,0,0] = False
    kspec = np.sqrt(kx[:,np.newaxis,np.newaxis]**2 + ky[np.newaxis,:,np.newaxis]**2 + kz[np.newaxis,np.newaxis,:]**2)
    kspec[0,0,0] = 1.
    muspec = np.absolute(kx[:,np.newaxis,np.newaxis])/kspec
    kspec[0,0,0] = 0.
    return kspec,muspec,indep  

def binpk(pkspec,nx,ny,nz,lx,ly,lz,kmin,kmax,nkbin):
    kspec,muspec,indep = getkspec(nx,ny,nz,lx,ly,lz)
    pkspec = pkspec[indep == True]
    kspec = kspec[indep == True]
    ikbin = np.digitize(kspec,np.linspace(kmin,kmax,nkbin+1))
    nmodes,pk = np.zeros(nkbin,dtype=int),np.full(nkbin,-1.)
    for ik in range(nkbin):
      nmodes[ik] = int(np.sum(np.array([ikbin == ik+1])))
      if (nmodes[ik] > 0):
        pk[ik] = np.mean(pkspec[ikbin == ik+1])
      if(ik==0)and(pk[ik]==-1.):pk[ik]=0.# the first kbin is not used. 
    return pk,nmodes 

def ConvMat_Blake(kmin,kmax,nkbin,kmaxc,nkbinc,nx,ny,nz,lx,ly,lz,w,Win,w2,Win2,PStype):
    # Initializations
    if(PStype!='crs'):nlconv = 3 
    if(PStype=='crs'):nlconv = 2
    dkc = (kmaxc-kmin)/nkbinc
    k1 = np.linspace(kmin,kmaxc-dkc,nkbinc)
    k2 = np.linspace(kmin+dkc,kmaxc,nkbinc)
    convmat = np.zeros((nlconv*(nkbin-1),nlconv*nkbinc))
    # Convolve series of unit vectors
    iconv = -1
    for imult in range(nlconv):
      for ik in range(nkbinc):
        iconv += 1
        print( 'Obtaining convolution for bin',iconv+1,'of',nlconv*nkbinc,'...')
        pk0,pk1,pk2, pk3,pk4 = 0.,0.,0.,0.,0.
        if(PStype!='crs'):
            if (imult == 0):   pk0 = 1.
            elif (imult == 1): pk2 = 1.          
            elif (imult == 2): pk4 = 1.
            kmin1,kmax1 = k1[ik],k2[ik]
            pk=[pk0, pk2,  pk4]
            pk0con, pk2con, pk4con = Conv_Fun(nx,ny,nz,lx,ly,lz,w,Win,w2,Win2,kmin,kmax,nkbin,kmin1,kmax1,pk,PStype)
            convmat[:,iconv] = np.concatenate((pk0con, pk2con, pk4con))
        if(PStype=='crs'): 
            if (imult == 0):   pk1 = 1.
            elif (imult == 1): pk3 = 1.          
            kmin1,kmax1 = k1[ik],k2[ik]
            pk=[pk1,  pk3]
            pk1con, pk3con = Conv_Fun(nx,ny,nz,lx,ly,lz,w,Win,w2,Win2,kmin,kmax,nkbin,kmin1,kmax1,pk,PStype)
            convmat[:,iconv] = np.concatenate((pk1con, pk3con))            
    kd=np.linspace(kmin,kmax,nkbin+1)[:-1]
    kc=np.concatenate((k1,np.array([k2[-1]])))
    return convmat,kd[:-1]+0.5*np.diff(kd),kc[:-1]+0.5*np.diff(kc)
# https://arxiv.org/pdf/1801.04969.pdf 
def Conv_Fun(nx,ny,nz,lx,ly,lz,w,Win,w2,Win2,kmin,kmax,nkbin,kmin1,kmax1,pk,PStype ):
    if(PStype!='crs'):
        pk0, pk2, pk4=pk
        nlmod = 3 
        nl=3
        uselp = np.full(nlmod,True,dtype=bool)
        if (pk0 == 0.): uselp[0] = False
        if (pk2 == 0.): uselp[1] = False
        if (pk4 == 0.): uselp[2] = False 
    if(PStype=='crs'): 
        pk1, pk3=pk
        nlmod = 2 
        nl=2
        uselp = np.full(nlmod,True,dtype=bool)
        if (pk1 == 0.): uselp[0] = False
        if (pk3 == 0.): uselp[1] = False
    # grid cells' possition and wave numbers. 
    dx,dy,dz = lx/nx,ly/ny,lz/nz
    x  = dx*np.arange(nx) - lx/2. + 0.5*dx
    y  = dy*np.arange(ny) - ly/2. + 0.5*dy
    z  = dz*np.arange(nz) - lz/2. + 0.5*dz
    kx = 2.*np.pi*np.fft.fftfreq(nx,d=dx)
    ky = 2.*np.pi*np.fft.fftfreq(ny,d=dy)
    kz = 2.*np.pi*np.fft.fftfreq(nz,d=dz)
    # Obtain spherical polar angles over the grid
    rgrid  = np.sqrt(x[:,np.newaxis,np.newaxis]**2 + y[np.newaxis,:,np.newaxis]**2 + z[np.newaxis,np.newaxis,:]**2)
    rtheta = np.arctan2(z[np.newaxis,np.newaxis,:],y[np.newaxis,:,np.newaxis])
    rphi   = np.where(rgrid>0.,np.arccos(x[:,np.newaxis,np.newaxis]/rgrid),0.)
    kgrid  = np.sqrt(kx[:,np.newaxis,np.newaxis]**2 + ky[np.newaxis,:,np.newaxis]**2 + kz[np.newaxis,np.newaxis,:]**2)
    ktheta = np.arctan2(kz[np.newaxis,np.newaxis,:],ky[np.newaxis,:,np.newaxis])
    kphi   = np.where(kgrid>0.,np.arccos(kx[:,np.newaxis,np.newaxis]/kgrid),0.)     
    mask   = (kgrid >= kmin1) & (kgrid < kmax1)
    # normalization factor:
    nw  = w*Win/np.sum(Win) 
    nwb = np.fft.fftn( nw )
    Nc  = nx*ny*nz
    I   = Nc*np.sum(nw**2)  
    nw2 = w2*Win2/np.sum(Win2) 
    nwb2= np.fft.fftn( nw2 )
    I2  = Nc*np.sum(nw2**2) 
    # calculate convolution:
    pk0con,pk1con,pk2con,pk3con,pk4con = np.zeros(nkbin-1),np.zeros(nkbin-1),np.zeros(nkbin-1),np.zeros(nkbin-1),np.zeros(nkbin-1)
    for il in range(nl): 
        if(PStype!='crs'): l = 2*il
        if(PStype=='crs'): l = 2*(il+1)-1
        pkcon = np.zeros((nx,ny,nz))
        for m in range(-l,l+1):
            Ylm_r=sph_harm(m,l,rtheta ,rphi)
            Ylm_k=sph_harm(m,l,ktheta ,kphi) 
            for ilp in range(nlmod):
              if (uselp[ilp]):     
                if(PStype!='crs'):lp = 2*ilp
                if(PStype=='crs'):lp = 2*(ilp+1)-1               
                norm = (4.*np.pi)**2/(2.*lp+1.) /I   
                norm2= (4.*np.pi)**2/(2.*lp+1.) /I2 
                Pl = np.zeros((nx,ny,nz))
                if(PStype!='crs'):
                    if(ilp == 0):   Pl[mask] = pk0
                    elif(ilp == 1): Pl[mask] = pk2
                    elif(ilp == 2): Pl[mask] = pk4 
                if(PStype=='crs'):
                    if(ilp == 0):   Pl[mask] = pk1
                    elif(ilp == 1): Pl[mask] = pk3    
                for mpr in range(-lp,lp+1):
                    Ylm_kp = sph_harm(mpr,lp,ktheta ,kphi)
                    Ylm_rp = sph_harm(mpr,lp,rtheta ,rphi)
                    Slmlm  = np.fft.fftn( nw *Ylm_r*np.conj(Ylm_rp) )
                    Slmlm2 = np.fft.fftn( nw2*Ylm_r*np.conj(Ylm_rp) )
                    PY     = np.fft.fftn( np.conj(Ylm_kp) * Pl )
                    nS     = np.fft.fftn(nwb * np.conj(Slmlm))
                    nS2    = np.fft.fftn(nwb2 * np.conj(Slmlm2))
                    pkcon  = pkcon + np.sqrt(norm*norm2) * np.real( Ylm_k * np.fft.ifftn( PY* 0.5*(nS+nS2) ) )  
        pkc ,nmodes = binpk(pkcon,nx,ny,nz,lx,ly,lz,kmin,kmax,nkbin)
        if(PStype!='crs'):
            if(il == 0):   pk0con = pkc[:-1]  
            elif(il == 1): pk2con = pkc[:-1]   
            elif(il == 2): pk4con = pkc[:-1]             
        if(PStype=='crs'): 
            if(il == 0):   pk1con = pkc[:-1]    
            elif(il == 1): pk3con = pkc[:-1]             
    if(PStype!='crs'):
        return   pk0con, pk2con, pk4con 
    if(PStype=='crs'):
        return   pk1con, pk3con 

  
    
     
    
    
#==============================================================================  
 
    

 



 




def Pk_CONV_multi(convmat,pkm,psyp):
    if(psyp!='crs'):
        pk0unconv,pk2unconv,pk4unconv=pkm
        pkunconvlst = np.concatenate((pk0unconv,pk2unconv,pk4unconv))
        pkconvlst   = np.dot(convmat,pkunconvlst)    
        nlmod=3
        N=len(pkconvlst)//nlmod 
        pk0conv,pk2conv,pk4conv = pkconvlst[:N],pkconvlst[N:2*N],pkconvlst[2*N:3*N]
        return pk0conv,pk2conv,pk4conv 
    if(psyp=='crs'):    
        pk1unconv,pk3unconv=pkm
        pkunconvlst = np.concatenate((pk1unconv,pk3unconv))
        pkconvlst   = np.dot(convmat,pkunconvlst)    
        nlmod=2
        N=len(pkconvlst)//nlmod 
        pk1conv,pk3conv = pkconvlst[:N],pkconvlst[N:2*N] 
        return pk1conv,pk3conv  
    
    
    
#==============================================================================    
    
    
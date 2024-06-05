import numpy as np
import math
from scipy.special import sph_harm
import scipy as sp
from scipy.interpolate import splev, splrep
from scipy import interpolate
from scipy import integrate
# Speed of light in km/s
LightSpeed = 299792.458
# Calculate distance:
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




def discret_rand(xpos,ypos,zpos,nx,ny,nz,lx,ly,lz, typ, wei=1.):
    if(typ=='weit'):
        datgrid,edges = np.histogramdd(np.vstack([xpos,ypos,zpos]).transpose(),bins=(nx,ny,nz),range=((-lx/2.,lx/2.),(-ly/2.,ly/2.),(-lz/2.,lz/2.)),weights=wei)
    else:    
        datgrid,edges = np.histogramdd(np.vstack([xpos,ypos,zpos]).transpose(),bins=(nx,ny,nz),range=((-lx/2.,lx/2.),(-ly/2.,ly/2.),(-lz/2.,lz/2.))   )     
    return datgrid
'''def Prep_winFUN(nrand,  fkp,rand_x,rand_y,rand_z,pvR,sigv,nx,ny,nz,lx,ly,lz,types ):
    vol = lx*ly*lz
    randgrid = discret_rand(rand_x,rand_y,rand_z,nx,ny,nz,lx,ly,lz,'no-weit' )
    wingrid = randgrid
    weigrid=np.zeros((nx,ny,nz))
    if(types=='den'):  
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    weigrid[i,j,k] = 1./(1.0+(randgrid[i,j,k]/nrand*float(nx*ny*nz)/vol)*fkp)
                    if(randgrid[i,j,k]==0): weigrid[i,j,k]=0
    if(types=='mom'):
        pverrgrid = discret_rand(rand_x,rand_y,rand_z,nx,ny,nz,lx,ly,lz,'weit',pvR)
        pverrgrid[randgrid > 0] /= randgrid[randgrid > 0]
        pverrgrid += sigv**2
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    weigrid[i,j,k] = 1./(pverrgrid[i,j,k]+(randgrid[i,j,k]/nrand*float(nx*ny*nz)/vol)*fkp)
                    if(randgrid[i,j,k]==0.): weigrid[i,j,k]=0
    return wingrid,weigrid''' 
def Prep_winFUN(nrand,  fkp,rand_x,rand_y,rand_z,pvR,sigv,nx,ny,nz,lx,ly,lz,types,weitpy ):
    if(types!='crs'):
        vol = lx*ly*lz
        randgrid = discret_rand(rand_x,rand_y,rand_z,nx,ny,nz,lx,ly,lz,'no-weit' )
        wingrid = randgrid
    if(types=='den'):  
        weigrid=np.zeros((nx,ny,nz))
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    weigrid[i,j,k] = 1./(1.0+(randgrid[i,j,k]/nrand*float(nx*ny*nz)/vol)*fkp)
                    if(randgrid[i,j,k]==0): weigrid[i,j,k]=0
                    if( weitpy  =='unity')and(randgrid[i,j,k]!=0): weigrid[i,j,k]=1.         
        return wingrid,weigrid
    if(types=='mom'):
        pverrgrid = discret_rand(rand_x,rand_y,rand_z,nx,ny,nz,lx,ly,lz,'weit',pvR)
        pverrgrid[randgrid > 0] /= randgrid[randgrid > 0]
        pverrgrid += sigv**2
        weigrid=np.zeros((nx,ny,nz))
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    weigrid[i,j,k] = 1./(pverrgrid[i,j,k]+(randgrid[i,j,k]/nrand*float(nx*ny*nz)/vol)*fkp)
                    if(randgrid[i,j,k]==0.): weigrid[i,j,k]=0
                    if( weitpy  =='unity')and(randgrid[i,j,k]!=0): weigrid[i,j,k]=1.
        return wingrid,weigrid
    if(types=='crs'):
        fkpg  ,fkpv   =fkp
        nrandg,nrandv = nrand      
        randg_x,randv_x=rand_x
        randg_y,randv_y=rand_y
        randg_z,randv_z=rand_z
        lgx,lvx=lx 
        lgy,lvy=ly
        lgz,lvz=lz
        ngx,nvx=nx
        ngy,nvy=ny 
        ngz,nvz=nz    
        volg      = lgx*lgy*lgz
        randgridg = discret_rand(randg_x,randg_y,randg_z,ngx,ngy,ngz,lgx,lgy,lgz,'no-weit' )
        wingridg  = randgridg
        weigridg=np.zeros((ngx,ngy,ngz))
        for i in range(ngx):
            for j in range(ngy):
                for k in range(ngz):
                    weigridg[i,j,k] = 1./(1.0+(randgridg[i,j,k]/nrandg*float(ngx*ngy*ngz)/volg)*fkpg) 
                    if(randgridg[i,j,k]==0): weigridg[i,j,k]=0 
                    if( weitpy  =='unity')and(randgridg[i,j,k]!=0): weigridg[i,j,k]=1.
        volv      = lvx*lvy*lvz
        randgridv = discret_rand(randv_x,randv_y,randv_z,nvx,nvy,nvz,lvx,lvy,lvz,'no-weit' )
        wingridv  = randgridv       
        pverrgrid = discret_rand(randv_x,randv_y,randv_z,nvx,nvy,nvz,lvx,lvy,lvz,'weit',pvR)
        pverrgrid[randgridv > 0] /= randgridv[randgridv > 0]
        pverrgrid += sigv**2       
        weigridv=np.zeros((nvx,nvy,nvz))
        for i in range(nvx):
            for j in range(nvy):
                for k in range(nvz):
                    weigridv[i,j,k] =1./(pverrgrid[i,j,k]+(randgridv[i,j,k]/nrandv*float(nvx*nvy*nvz)/volv)*fkpv) #1./(pverrgrid[i,j,k]+(randgrid[i,j,k]/nrand*float(nx*ny*nz)/vol)*fkp)
                    if(randgridv[i,j,k]==0.): weigridv[i,j,k]=0 
                    if( weitpy  =='unity')and(randgridv[i,j,k]!=0): weigridv[i,j,k]=1.
        weigrid=np.zeros((nvx,nvy,nvz))            
        for i in range(nvx):
            for j in range(nvy):
                for k in range(nvz):                    
                    weigrid[i,j,k]  = 1.0/np.sqrt( 1.0 + weigridg[i,j,k]*weigridv[i,j,k]*((randgridg[i,j,k]/nrandg*randgridv[i,j,k]/nrandv*float(ngx*ngy*ngz)/volg*float(nvx*nvy*nvz)/volv)*fkpv*fkpg)**2)
                    if(randgridv[i,j,k]==0.)and(randgridg[i,j,k]==0): weigrid[i,j,k]=0
                    if(randgridv[i,j,k]==0.)and(randgridg[i,j,k]!=0): weigrid[i,j,k]=weigridg[i,j,k]
                    if(randgridv[i,j,k]!=0.)and(randgridg[i,j,k]==0): weigrid[i,j,k]=weigridv[i,j,k]
                    if( weitpy  =='unity'):
                        if(randgridv[i,j,k]!=0)or(randgridg[i,j,k]!=0): weigrid[i,j,k]=1.    
                    #weigrid  = 1.0/np.sqrt( 1.0 + weigridg*weigridv*((randgridg/nrandg *float(ngx*ngy*ngz)/volg )*fkpv*fkpg)**2)
        return wingridg,wingridv,weigrid

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
    #print 'Binning in angle-averaged bins...'
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
    if(PStype!='crs'):nlconv = 3 # Number of multipoles to include in convolution
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
      #if (uselp[il ]):    # ------ >>>>> ??????
        if(PStype!='crs'): l = 2*il
        if(PStype=='crs'): l = 2*(il+1)-1
        pkcon = np.zeros((nx,ny,nz))
        for m in range(-l,l+1):
            Ylm_r=sph_harm(m,l,rtheta ,rphi)
            Ylm_k=sph_harm(m,l,ktheta ,kphi) 
            for ilp in range(nlmod):
              if (uselp[ilp]):     # ------ >>>>> ??????
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
 
    
def PS_extender(k,Pk,k_start,k_end):
    #k_store=k
    XX=k[k_start:k_end]
    YY=Pk[k_start:k_end]
    x=np.log10(XX)
    y=np.log10(YY)
    KK,BB=np.polyfit(x,y,1)
    dks=0.001
    krs=np.zeros(int((10.0-max(k))/dks))
    for i in range(int((10.0-max(k))/dks)):
        krs[i]=max(k)+(i+1)*dks
    KRAND=k
    PKRAND=Pk
    k=np.zeros(len(KRAND)+int((10.0-max(k))/dks))
    Pk=np.zeros(len(k))
    k[0:len(KRAND)]=KRAND
    k[len(KRAND):len(KRAND)+int((10.0-max(k))/dks)]=krs
    Pk[0:len(KRAND)]=PKRAND
    Pk[len(KRAND):len(k)]=10.0**(KK*np.log10(krs)+BB)
    return k,Pk

def ConvMat_Ross(Ki,kjmin,kjmax,nkbinc,WF_rand,Ncosa,epsilon_par):  
    #kjmin=min(Ki);kjmax=max(Ki)
    Kj=np.linspace(kjmin,kjmax,nkbinc+1)
    
    dcosa=2.0/ Ncosa   
    Nki=len(Ki)  ;  #ki_min=min(Ki)  ;  ki_max=max(Ki)      
    Nkj=len(Kj)  ;  kj_min=min(Kj)  ;  kj_max=max(Kj)   
    WF_kikj=np.zeros((Nki,Nkj))  #WF_kikj[KI,KJ]
    
    WF_rand[0][0]=0.0  ;  WF_rand[1][0]=1.0
    Pwin_spline = sp.interpolate.splrep(WF_rand[0],WF_rand[1], s=0) 
    
    Neps=epsilon_par[2];  eps_min=epsilon_par[0];  eps_max=epsilon_par[1]
    deps=(eps_max-eps_min)/Neps
    # https://wenku.baidu.com/view/719b47fa0740be1e640e9a97.html
    cosa=np.zeros((Ncosa+1))
    for i_cosa in range(Ncosa+1):
        cosa[i_cosa]=-1.0+i_cosa*dcosa
        
    eps=np.zeros((Neps+1))
    for i_eps in range(Neps+1):
        eps[i_eps]=eps_min+i_eps*deps
    Pwin = sp.interpolate.splev(eps, Pwin_spline, der=0)
    
    for i in range(Nki):
       FUN=np.zeros(  (Nkj,Neps+1,Ncosa+1)  )# FUN[Nkj,Neps+1,Ncosa+1]
       for i_eps in range(Neps+1):  
          for i_cosa in range(Ncosa+1):
             if (np.abs(Ki[i]**2+eps[i_eps]**2-2.0*Ki[i]*eps[i_eps]*cosa[i_cosa])<10.e-15) and (cosa[i_cosa]==1.0):
                 R_eps=0.0
             else:
                 R_eps= np.sqrt(Ki[i]**2+eps[i_eps]**2-2.0*Ki[i]*eps[i_eps]*cosa[i_cosa])
             bins = int( Nkj*(R_eps-kj_min)/(kj_max-kj_min) )
             if((bins >= 0) and (bins <Nkj)): 
                 FUN[bins,i_eps,i_cosa]= eps[i_eps]* eps[i_eps] * Pwin[i_eps]  
       #--------------------           
       for j in range(Nkj):    
          F1=FUN[j,0,0]+FUN[j,Neps,0]+FUN[j,0,Ncosa]+FUN[j,Neps,Ncosa] 
          F2=sum(FUN[j,:,0])    
          F3=sum(FUN[j,:,Ncosa])    
          F4=sum(FUN[j,0,:])    
          F5=sum(FUN[j,Neps,:]) 
          F6=sum(sum(FUN[j,1:Neps,1:Ncosa]))
          #F6=0.
          #for i_eps in range(Neps-1): 
          #   for i_cosa in range(Ncosa-1):
          #      F6=F6+FUN[j,i_eps+1,i_cosa+1]
          WF_kikj[i,j]=deps*dcosa*(0.25*F1+0.5*(F2+F3+F4+F5)+F6)           
    
    # normalization WF_kikj[i,j]:
    for i in range(len(WF_kikj[:,1])):
        WF_kikj[i,:]=WF_kikj[i,:]/sum(WF_kikj[i,:])
        
    return WF_kikj,Kj

#============================================

def Pk_CONV_ross(Ki,WFkikj,Pkj,k_random,WF_random):    
    
    Psum=np.zeros(len(Ki))
    for i in range(len(Psum)):
        Psum[i]=sum(WFkikj[i,:]*Pkj)

    Pwin_spline = sp.interpolate.splrep(k_random,WF_random, s=0)
    PwinKi  = sp.interpolate.splev(Ki, Pwin_spline, der=0)
    Pwin0 = sp.interpolate.splev(0.0, Pwin_spline, der=0)

    P0 = sum(WFkikj[0,:]*Pkj) / Pwin0
    Pm = Psum - P0*PwinKi
    return Pm


#============================================    
    





def ConvMat_Beutler( Xr,Yr,Zr, randoms_weights,Ngrid,wnorm,NrandRat,kmin,kmax,nkbin,kminc,kmaxc,nkbinc,ells ,boxsize,boxsizes,frac_nyq,pstype,WFCon_type=False):
  if(WFCon_type=='Beutler'):  
    from pypower import CatalogSmoothWindow,PowerSpectrumSmoothWindow,PowerSpectrumSmoothWindowMatrix,Projection 
    kout    = np.linspace(kmin,kmax,nkbin+1)[:-1]
    kout= kout[:-1]+0.5*np.diff(kout)
    edges = {'step': 2. * np.pi / boxsize}
    if(pstype!='crs'):
        wnorm=wnorm*NrandRat**2
        randoms_positions=np.zeros((len(Xr),3))
        randoms_positions[:,0]=Xr ; randoms_positions[:,1]=Yr ; randoms_positions[:,2]=Zr
        projss=[Projection(ell=0, wa_order=0), Projection(ell=2, wa_order=0), Projection(ell=4, wa_order=0), Projection(ell=6, wa_order=0),  Projection(ell=8, wa_order=0), 
                Projection(ell=1, wa_order=1), Projection(ell=3, wa_order=1), Projection(ell=5, wa_order=1), Projection(ell=7, wa_order=1)  ]
        window_large = CatalogSmoothWindow(randoms_positions1=randoms_positions, randoms_weights1=randoms_weights,
                                   projs=projss, nmesh=Ngrid,wnorm=wnorm, edges=edges, boxsize=boxsize,  boxcenter=np.array([0,0,0]),position_type='pos', dtype='f8').poles
        window_small = CatalogSmoothWindow(randoms_positions1=randoms_positions, randoms_weights1=randoms_weights,
                                   projs=projss, nmesh=Ngrid,wnorm=wnorm, edges=edges, boxsize=boxsizes, boxcenter=np.array([0,0,0]), position_type='pos', dtype='f8').poles
        window  = PowerSpectrumSmoothWindow.concatenate_x(window_large, window_small, frac_nyq=frac_nyq) # Here we remove the range above 0.9 Nyquist (which may not be reliable) 
        projsin = [Projection(ell=0, wa_order=0), Projection(ell=2, wa_order=0), Projection(ell=4, wa_order=0),  
                   Projection(ell=1, wa_order=1), Projection(ell=3, wa_order=1) ] 
    if(pstype=='crs'):
        Xr1,Xr2=Xr
        Yr1,Yr2=Yr
        Zr1,Zr2=Zr
        rand1=np.zeros((len(Xr1),3));rand2=np.zeros((len(Xr2),3))
        rand1[:,0]=Xr1 ; rand1[:,1]=Yr1 ; rand1[:,2]=Zr1
        rand2[:,0]=Xr2 ; rand2[:,1]=Yr2 ; rand2[:,2]=Zr2
        rweit1,rweit2=randoms_weights  
        Nr1,Nr2=NrandRat       
        wnorm=wnorm*Nr1*Nr2  
        projss=[Projection(ell=0, wa_order=0), Projection(ell=1, wa_order=0), Projection(ell=2, wa_order=0),Projection(ell=3, wa_order=0),  Projection(ell=4, wa_order=0), Projection(ell=5, wa_order=0),  
                Projection(ell=0, wa_order=1), Projection(ell=1, wa_order=1), Projection(ell=2, wa_order=1),Projection(ell=3, wa_order=1),  Projection(ell=4, wa_order=1), Projection(ell=5, wa_order=1), Projection(ell=6, wa_order=1)  ]
        window_large = CatalogSmoothWindow(randoms_positions1=rand1,randoms_positions2=rand2, randoms_weights1=rweit1,randoms_weights2=rweit2,
                                   projs=projss, nmesh=Ngrid,wnorm=wnorm, edges=edges, boxsize=boxsize,  boxcenter=np.array([0,0,0]), position_type='pos', dtype='f8').poles
        window_small = CatalogSmoothWindow(randoms_positions1=rand1,randoms_positions2=rand2, randoms_weights1=rweit1,randoms_weights2=rweit2,
                                   projs=projss, nmesh=Ngrid,wnorm=wnorm, edges=edges, boxsize=boxsizes, boxcenter=np.array([0,0,0]), position_type='pos', dtype='f8').poles
        window  = PowerSpectrumSmoothWindow.concatenate_x(window_large, window_small, frac_nyq=frac_nyq) # Here we remove the range above 0.9 Nyquist (which may not be reliable)    
        projsin = [   Projection(ell=0, wa_order=0), Projection(ell=1, wa_order=1),Projection(ell=2, wa_order=0),  Projection(ell=3, wa_order=1)  ]  
    sep     = np.geomspace(1e-4, 4e3, nkbinc *4 )
    wawm    = PowerSpectrumSmoothWindowMatrix(kout, projsin=projsin, projsout=ells , window=window,kin_lim=(kminc,kmaxc),   sep=sep )
    kin     = wawm.xin[0]
    wawm.resum_input_odd_wide_angle() 
    if(pstype=='crs'):
        return  -((wawm.value).T)[:,:int(len(ells)*len(kin))], kout,kin
    else:
        return  (wawm.value).T, kout,kin



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
    
    
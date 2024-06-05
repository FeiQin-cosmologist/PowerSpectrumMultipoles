!==============================================================================
!                              INTRODUCTION
!------------------------------------------------------------------------------
! Please run the following: 
! gfortran -I/usr/local/include -L/usr/local/lib /Users/fei/WSP/Scie/Proj12/Prog/2_PS/PSest/2_PSrand.f90 /Users/fei/WSP/Scie/Proj12/Prog/2_PS/PSest/momPS.f90  -O3 -o Proj -lm -lgsl -lgslcblas -lfftw3
! ./Proj
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Warnning: before you running this code, you shoudl firstly open the 'momPS.f90' code
! to update the grid parameters.
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! This code will compute the density and momentum power spectrum of the random catalogue,
! it will be used to calculate the window function convolution. 
! Input_rand      : the same as the 'outdir_surveyrand' in the '1_PrepareData.py' .  
! output          : the output directory of the power spectrum of random data. 
! OmegaM          : the density parameter
! OmegaA          : the dark energy parameter. 
! obs_pos         : observer's positions. 
! The above 5 variables are the only things you need to modify in this code.
! Do not modify the rest of the code. 
!------------------------------------------------------------------------------
program WF_PS!
use randPS_interface
implicit none!
!##############################################################################













! Please set the cosmological parameters and the directories of the input and output
! Do not modify the rest of the code. 
integer(8),parameter ::    do_mom=1
character(1000)      ::    file_dir  = '/Users/fei/WSP/Scie/Proj12/Data/Prod/'
real(8),parameter    ::    OmegaM = 0.3121
real(8),parameter    ::    OmegaA = 1.- OmegaM
real(8),dimension(3),parameter::    obs_pos=(/ 0., 0., 0./)!mpc/h















!#################################  The  End   ################################






































































































! Do not change the following code:
!==============================================================================
!
!                                  CODE
!
!------------------------------------------------------------------------------
!--------------------------------------------------------------------------
real(8),parameter ::    mathPi=3.1415926535897932
real(8),parameter ::    Cl = 299792.458
real(8),parameter ::    Hub    = 100.!67.0
integer(8)::ndataR
real(8),allocatable,dimension(:):: Xr,Yr,Zr,nbr,wr
character(1000) ::output_dir,Input_rand,output
! the above is for sub.PS.do not change them.
!--------------------------------------------------------------------------
!Users can define variables in the
!following space, the above variables
!can not be change since they are
!related to the subroutines, if the
!above variables are changed, it will
!cause errors when compile the code file.
character(3)::varch
integer :: i
real(8),allocatable,dimension(:):: longitude_rand,latitude_rand,Vcmb_rand
real(8),allocatable,dimension(:):: nbar_rand,wden_rand,wmom_rand
real(8),allocatable,dimension(:):: disr
real(8):: var
!==========================    MAIN PROGRAM    ============================


if(do_mom==1)then
    Input_rand= trim(adjustl(file_dir))//trim(adjustl('Rand_pv'))  
endif
if(do_mom==0)then
    Input_rand= trim(adjustl(file_dir))//trim(adjustl('Rand_gal'))  
endif
output    = trim(adjustl(file_dir))//trim(adjustl('Rand_'))


! 1. Read random points:
open (unit=134,file=Input_rand,form='formatted')
read(134,*) varch ,ndataR
allocate(longitude_rand(ndataR));allocate(latitude_rand(ndataR));
allocate(Vcmb_rand(ndataR));allocate(nbar_rand(ndataR));
allocate(wden_rand(ndataR));allocate(wmom_rand(ndataR));
DO i=1,ndataR
if(do_mom==1)then
    read(134,*)longitude_rand(i),latitude_rand(i),Vcmb_rand(i),&
    nbar_rand(i),wmom_rand(i)
endif  
if(do_mom==0)then
    read(134,*)longitude_rand(i),latitude_rand(i),Vcmb_rand(i),&
    nbar_rand(i),wden_rand(i) 
endif      
ENDDO
close(134)



! 2. calculate distance:
! random:
allocate(disr(ndatar));
disr=0.
do i=1,ndatar
    call CosmDist(300,OmegaM,OmegaA,Hub,Cl,Vcmb_rand(i)/cl,Disr(i))
enddo
Xr  =disr*cos(latitude_rand/180.*mathPi)*cos(longitude_rand/180.*mathPi)
Yr  =disr*cos(latitude_rand/180.*mathPi)*sin(longitude_rand/180.*mathPi)
Zr  =disr*sin(latitude_rand/180.*mathPi)




! 3.calculate power spectrum:
if(do_mom==1)output_dir=trim(adjustl(output))//trim(adjustl('mom')) 
if(do_mom==0)output_dir=trim(adjustl(output))//trim(adjustl('den')) 

allocate(nbr(ndataR));allocate(wr(ndataR))

if(do_mom==1)then
    wR  =wmom_rand
endif
if(do_mom==0)then
    wR  =wden_rand
endif

print*,'Randoms'
nbr=nbar_rand
call rand_PS(do_mom,obs_pos,ndataR,nbr,Xr,Yr,Zr,wr,output_dir)





end program WF_PS
!==================





!--------------------------    subroutine.2    ----------------------------
subroutine CosmDist(Nprec,OmegaM,OmegaA,Hub,Cl,Rsf,Dis)
implicit none
integer,intent(in) ::  Nprec
real(8),intent(in) ::  OmegaM,OmegaA,Hub,Cl
real(8),intent(in) ::  Rsf ! input redshift data z rather than cz.
real(8),intent(out) :: Dis ! output comoving distance data.
real(8) :: z(Nprec+1),comdis(Nprec+1)
integer :: j,k

do j=1,Nprec+1
z(j)=0.+(j-1) * abs(Rsf - 0.) / Nprec
enddo
z(1)=0.
z(Nprec+1)=Rsf

Dis=0.
comdis=Cl/Hub * 1./sqrt(OmegaM*(1+z)**3+OmegaA);
Dis=sum((comdis(1:Nprec)+comdis(2:(Nprec+1)))*abs(z(2)-z(1))/2);

end subroutine CosmDist

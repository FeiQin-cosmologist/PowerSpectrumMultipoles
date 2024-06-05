!==============================================================================
!                              INTRODUCTION
!------------------------------------------------------------------------------
! Please run the following: 
! gfortran -I/usr/local/include -L/usr/local/lib /Users/fei/WSP/Scie/Proj12/Prog/2_PS/PSest/1_PS.f90 /Users/fei/WSP/Scie/Proj12/Prog/2_PS/PSest/momPS.f90  -O3 -o Proj -lm -lgsl -lgslcblas -lfftw3
! ./Proj
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Warnning: before you running this code, you shoudl firstly open the 'momPS.f90' code
! to update the grid parameters.
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Reference to GNU Lib: 
! https://people.sc.fsu.edu/~jburkardt/f_src/fftw3/fftw3.html
! This code will compute the density and momentum power spectrum.
! OmegaM          : the density parameter
! OmegaA          : the dark energy parameter. 
! obs_pos         : observer's positions. 
! PS_Multi        :   do you want to calculate non-zeros multipoles? if so, set it to be 'Yes'. If set it to be 'NO', only calculate l=0. 
!                   set it to be 'ALL' if you want to calculate all multipoles, including the zeros multipoles.
! Nmock        :   how many mocks you will use? 
! Nrand        :   how many random catalogue for mocks you will use?  Nrand=1 or Nmock 
! Input_dir        = outdir_galsurvey
! Input_dir_pv     = outdir_pvsurvey
! Input_rand_dir   = outdir_galsurveyrand
! Input_mock_dir   = outdir_galmock
! Input_mock_dir_pv= outdir_pvmock
! Input_rand_for_mock_dir = outdir_galmockrand
! Output_dir       = the output directory of the power spectrum of survey data.
! Output_mock_dir  = the output directory of the power spectrum of mock data.
! The above 14 variables are the only things you need to modify in this code.
! Do not modify the rest of the code. 
!------------------------------------------------------------------------------
! Define the Variables
program denmomPS!
use momPS_interface
implicit none!
real(8),parameter ::    mathPi = 3.1415926535897932
real(8),parameter ::    Cl     = 299792.458
real(8),allocatable,dimension(:)::longitude_rand  ,latitude_rand  ,Vcmb_rand  ,nbr  ,wr  ,disr  ,xr  ,yr  ,zR 
real(8),allocatable,dimension(:)::longitude_randmk,latitude_randmk,Vcmb_randmk,nbrmk,wrmk,disrmk,xrmk,yrmk,zRmk     
real(8),allocatable,dimension(:)::longitude,latitude,Vcmb,nbden,wden,disz,x,y,z,NAN,wrho
real(8),allocatable,dimension(:)::longitudepv,latitudepv,Vcmbpv,logd,nbmom,wmom,diszpv,xpv,ypv,zpv,Vmod,Vpec 
character(2000)::Input_dir,Input_dir_pv,Input_rand_dir,Output,Output_dir,Input_dirF,Input_rand_for_mock_dir
character(2000)::Input_mock_dir_pv,input_mock_dir,output_mock_dir,file_dir
character(10)::varch ,file_Ino,PS_Multi 
integer(8)  :: ndataR,ndata,ndatapv,do_mom,do_mock,icod,i_mockS,i_mockSs,NMOCK,i_mock,ndataRmk,nrand,i,tmp
real(8)::deccel,var,iRUN
!##############################################################################


 


 



! Please set the cosmological parameters :  
real(8),parameter ::                OmegaM = 0.3121
real(8),parameter ::                OmegaA = 1.- OmegaM
real(8),parameter ::                Hub    = 100.!67.0
real(8),dimension(3),parameter::    obs_pos   = (/ 0., 0., 0./)!mpc/h
! Please setting the directories of the input and output.
! These are the only thing you need to modify in this code.
! Do not modify the rest of the code. 
PS_Multi         = 'YES' ! 'ALL'  !    'NO'   !       do you want to calculate l=1,2,3,4 ? 
Nmock            =  1024!2048           ! how many mocks you will use? 
Nrand            =  1
file_dir         = '/Users/fei/WSP/Scie/Proj12/Data/Prod/'
! if Nrand=1, this measn you only have one random catalogue, in this case
! you do not need to set Input_rand_for_mock_dir.
! just leave it empty  Input_rand_for_mock_dir=''. 









!#################################  The  End   ################################












































 




! Do not change the following code:
!==============================================================================
!
!                                  CODE
!
!------------------------------------------------------------------------------
Input_dir        = trim(adjustl(file_dir))//trim(adjustl('Data_gal'))
Input_dir_pv     = trim(adjustl(file_dir))//trim(adjustl('Data_pv'))
Input_rand_dir   = trim(adjustl(file_dir))//trim(adjustl('Rand_gal'))
Input_mock_dir   = trim(adjustl(file_dir))//trim(adjustl('Mock_gal'))
Input_mock_dir_pv= trim(adjustl(file_dir))//trim(adjustl('Mock_pv'))
Input_rand_for_mock_dir = trim(adjustl(file_dir))//trim(adjustl('MockRand_gal'))
Output_dir       = trim(adjustl(file_dir))//trim(adjustl('Data_'))
Output_mock_dir  = trim(adjustl(file_dir))//trim(adjustl('Mock_'))
! 1: Read in the random catalogue : -------------------------------------------
open (unit=134,file=Input_rand_dir,form='formatted')
read(134,*) varch ,ndataR
allocate(longitude_rand(ndataR));allocate(latitude_rand(ndataR));
allocate(Vcmb_rand(ndataR));allocate(nbr(ndataR));allocate(wr(ndataR))
DO i=1,ndataR
    read(134,*)longitude_rand(i),latitude_rand(i),Vcmb_rand(i),nbr(i),wr(i)
ENDDO
close(134)
allocate(disr(ndatar));
disr=0.
do i=1,ndatar
    call CosmDist(100,OmegaM,OmegaA,Hub,Cl,Vcmb_rand(i)/cl,Disr(i))
enddo
allocate(xr(ndatar));allocate(yr(ndataR));allocate(zR(ndataR))
Xr  =disr*cos(latitude_rand/180.*mathpi)*cos(longitude_rand/180.*mathpi)
Yr  =disr*cos(latitude_rand/180.*mathpi)*sin(longitude_rand/180.*mathpi)
Zr  =disr*sin(latitude_rand/180.*mathpi)
         

! 2: Calculating the power spectrum : -----------------------------------------
tmp=-1
Do icod= 1,1!6

    do_mom = 2 ; do_mock= 1   ; iRUN = 1  
 
    !if(icod==1)then
    !    do_mom = 1 ; do_mock= 0  
    !endif
    !if(icod==2)then
    !    do_mom = 1 ; do_mock= 1  
    !endif
    !if(icod==3)then
    !    do_mom = 0 ; do_mock= 0  
    !endif
    !if(icod==4)then
    !    do_mom = 0 ; do_mock= 1  
    !endif
    !if(icod==5)then
    !    do_mom = 2 ; do_mock= 0  
    !endif
    !if(icod==6)then
    !    do_mom = 2 ; do_mock= 1  
    !endif
    


    !########################   survey ######################################## 
    IF(do_mock==0)then
        print*,'Real Survey'
        if((do_mom == 0) .or. (do_mom==2))then          
            ! 2.1. Read data:
            open (unit=11,file=Input_dir,form='formatted')
            read(11,*) varch ,ndata
            allocate(longitude(ndata));allocate(latitude(ndata));allocate(Vcmb(ndata));
            allocate(nbden(ndata)); allocate(wden(ndata)) 
            DO i=1,ndata
                read(11,*)longitude(i),latitude(i),Vcmb(i),nbden(i),wden(i) 
            ENDDO
            close(11)
            ! 2.2. Calculate distance:
            Output=trim(adjustl(Output_dir))//trim(adjustl('den')) 
            allocate(x(ndata)); allocate(y(ndata)) ; allocate(z(ndata));allocate(disz(ndata));
            disz=0. 
            do i=1,ndata
                call CosmDist(100,OmegaM,OmegaA,Hub,Cl,Vcmb(i)/cl,Disz(i))
            enddo         
            X   =disz*cos(latitude/180.*mathpi)*cos(longitude/180.*mathpi)
            Y   =disz*cos(latitude/180.*mathpi)*sin(longitude/180.*mathpi)
            Z   =disz*sin(latitude/180.*mathpi)
            ! 2.3. Calculate density power:
            if(do_mom==0)then 
                allocate(NAN(ndata)) 
                call cosmic_PS(do_mom ,obs_pos, PS_Multi,output,   &
                                       tmp,NAN,NAN,NAN,NAN,NAN,NAN, NAN, &                                                                    
                                       ndata ,  nbden,  X,  Y,  Z,  wden, &
                                       ndataR,  nbr,    Xr, Yr, Zr, wr )                          
                deallocate(NAN)                                                 
                deallocate(longitude);deallocate(latitude);deallocate(Vcmb); deallocate(nbden); 
                deallocate(x); deallocate(y) ; deallocate(z);deallocate(disz);deallocate(wden);    
            endif       
        endif  
        if((do_mom == 1) .or. (do_mom==2))then   
            ! 2.1. Read data:
            open (unit=125,file=Input_dir_pv,form='formatted')
            read(125,*) varch ,ndatapv
            allocate(longitudepv(ndatapv));allocate(latitudepv(ndatapv));allocate(Vcmbpv(ndatapv));
            allocate(logd(ndatapv));allocate(nbmom(ndatapv));allocate(wmom(ndatapv));allocate(wrho(ndatapv)) 
            DO i=1,ndatapv
                read(125,*)longitudepv(i),latitudepv(i),Vcmbpv(i),logd(i),nbmom(i),wmom(i),wrho(i)
            ENDDO
            close(125)                   
            ! 2.2. Calculate distance:
            Output=trim(adjustl(Output_dir))//trim(adjustl('mom')) 
            allocate(xpv(ndatapv)); allocate(ypv(ndatapv)) ; allocate(zpv(ndatapv));allocate(diszpv(ndatapv));
            diszpv=0. 
            do i=1,ndatapv
                call CosmDist(100,OmegaM,OmegaA,Hub,Cl,Vcmbpv(i)/cl,Diszpv(i))
            enddo         
            Xpv   =diszpv*cos(latitudepv/180.*mathpi)*cos(longitudepv/180.*mathpi)
            Ypv   =diszpv*cos(latitudepv/180.*mathpi)*sin(longitudepv/180.*mathpi)
            Zpv   =diszpv*sin(latitudepv/180.*mathpi)    
            allocate(Vmod(ndatapv));allocate(Vpec(ndatapv))
            deccel=3.0*OmegaM/2.0 - 1.0
            vmod  = Vcmbpv*(1.0+0.5*(1.0-deccel)*Vcmbpv/cl-(2.0-deccel-3.0*deccel*deccel)*Vcmbpv/cl*Vcmbpv/cl/6.0)
            vpec  =  log(10.)*Vmod/(1.0+Vmod/cl) * logd 
            ! 2.3. Calculate momentum power: 
            if(do_mom==1)then 
                allocate(NAN(ndata)) 
                call cosmic_PS(do_mom ,obs_pos, PS_Multi,output,  &
                                       ndatapv,nbmom,Xpv,Ypv,Zpv,Vpec,wmom, wrho,          &                                                                     
                                       tmp ,  NAN,  NAN,  NAN,  NAN,  NAN,   &
                                       tmp,  NAN,    NAN, NAN, NAN, NAN )                                                                                                               
                deallocate(NAN)                                                  
                deallocate(longitudepv);deallocate(latitudepv);deallocate(Vcmbpv);
                deallocate(logd); deallocate(nbmom);deallocate(Vpec);deallocate(diszpv);deallocate(wrho);
                deallocate(wmom); deallocate(xpv); deallocate(ypv) ; deallocate(zpv);deallocate(Vmod); 
            endif                         
        endif  
        if(do_mom==2)then  
            Output=trim(adjustl(Output_dir))//trim(adjustl('crs')) 
            call cosmic_PS(do_mom ,obs_pos, PS_Multi,output,  &
                                   ndatapv, nbmom,  Xpv,Ypv,Zpv,Vpec,wmom,wrho,  &                                                                 
                                   ndata ,  nbden,  X,  Y,  Z,  wden ,  &
                                   ndataR,  nbr,    Xr, Yr, Zr, wr )                                                            
            deallocate(longitudepv);deallocate(latitudepv);deallocate(Vcmbpv);
            deallocate(logd); deallocate(nbmom);deallocate(Vpec);deallocate(diszpv);
            deallocate(wmom); deallocate(xpv); deallocate(ypv) ; deallocate(zpv);deallocate(Vmod); 
            deallocate(longitude);deallocate(latitude);deallocate(Vcmb); deallocate(nbden); 
            deallocate(x); deallocate(y) ; deallocate(z);deallocate(disz);deallocate(wden);deallocate(wrho);                                                                                                                                          
        endif
    ENDIF 



    !########################   mocks #########################################
    IF(do_mock==1)then      
        DO i_mockSs=1,NMOCK
        
        
            i_mockS=i_mockSs+NMOCK*iRUN     
            
                         
            print*,' i_mock    =',i_mockS-1
            i_mock=(i_mockS-1) 
            write(file_Ino,'(i10)')i_mock
            
            ! if we have more than one random catalogue-----    
            if((Nrand .gt. 1).and.(do_mom .ne. 1))then  ! read randoms for mocks:
                Input_dirF =trim(adjustl(Input_rand_for_mock_dir))//trim(adjustl(file_Ino)) 
                open (unit=1348,file=Input_dirF,form='formatted')
                read(1348,*) varch ,ndataRmk
                allocate(longitude_randmk(ndataRmk));allocate(latitude_randmk(ndataRmk));
                allocate(Vcmb_randmk(ndataRmk));allocate(nbrmk(ndataRmk));allocate(wrmk(ndataRmk))
                DO i=1,ndataRmk
                    read(1348,*)longitude_randmk(i),latitude_randmk(i),Vcmb_randmk(i),nbrmk(i),wrmk(i)
                ENDDO
                close(1348)
                allocate(disrmk(ndataRmk));
                disrmk=0.
                do i=1,ndataRmk
                    call CosmDist(100,OmegaM,OmegaA,Hub,Cl,Vcmb_randmk(i)/cl,disrmk(i))
                enddo
                allocate(xrmk(ndataRmk));allocate(yrmk(ndataRmk));allocate(zRmk(ndataRmk))
                Xrmk  =disrmk*cos(latitude_randmk/180.*mathpi)*cos(longitude_randmk/180.*mathpi)
                Yrmk  =disrmk*cos(latitude_randmk/180.*mathpi)*sin(longitude_randmk/180.*mathpi)
                Zrmk  =disrmk*sin(latitude_randmk/180.*mathpi)
                deallocate(longitude_randmk);deallocate(latitude_randmk);deallocate(disrmk);deallocate(Vcmb_randmk) 
            endif    
            
            ! density power :   --------------------------
            if((do_mom == 0) .or. (do_mom==2))then          
                ! 2.1. Read data:
                Input_dirF =trim(adjustl(Input_mock_dir))//trim(adjustl(file_Ino)) 
                open (unit=14,file=Input_dirF,form='formatted')
                read(14,*) varch ,ndata
                allocate(longitude(ndata));allocate(latitude(ndata));allocate(Vcmb(ndata));
                allocate(nbden(ndata)); allocate(wden(ndata)) 
                DO i=1,ndata
                    read(14,*)longitude(i),latitude(i),Vcmb(i),nbden(i),wden(i) 
                ENDDO
                close(14)
                ! 2.2. Calculate distance:
                Output= trim(adjustl(output_mock_dir))//  &
                        trim(adjustl('den'))//trim(adjustl(file_Ino)) 
                allocate(x(ndata)); allocate(y(ndata)) ; allocate(z(ndata));allocate(disz(ndata));
                disz=0. 
                do i=1,ndata
                    call CosmDist(100,OmegaM,OmegaA,Hub,Cl,Vcmb(i)/cl,Disz(i))
                enddo         
                X   =disz*cos(latitude/180.*mathpi)*cos(longitude/180.*mathpi)
                Y   =disz*cos(latitude/180.*mathpi)*sin(longitude/180.*mathpi)
                Z   =disz*sin(latitude/180.*mathpi)
                ! 2.3. Calculate density power:  
                if(do_mom==0)then 
                    allocate(NAN(ndata)) 
                    if(Nrand == 1)then
                        call cosmic_PS(do_mom ,obs_pos, PS_Multi,output,  &
                                               tmp,NAN,NAN,NAN,NAN,NAN,NAN,NAN,  &                                                                    
                                               ndata ,  nbden,  X,  Y,  Z,  wden,   &
                                               ndataR,  nbr,    Xr, Yr, Zr, wr )  
                    else
                        call cosmic_PS(do_mom ,obs_pos, PS_Multi,output,   &
                                               tmp,NAN,NAN,NAN,NAN,NAN,NAN,NAN,  &                                                                    
                                               ndata ,  nbden,  X,  Y,  Z,  wden,   &
                                               ndataRmk,  nbrmk,    Xrmk, Yrmk, Zrmk, wrmk )    
                        deallocate(nbrmk);deallocate(Xrmk);deallocate(Yrmk);deallocate(Zrmk);deallocate(wrmk)                 
                    endif                                                
                    deallocate(NAN)                                                     
                    deallocate(longitude);deallocate(latitude);deallocate(Vcmb); deallocate(nbden); 
                    deallocate(x); deallocate(y) ; deallocate(z);deallocate(disz);deallocate(wden); 
                endif               
            endif        
    
            ! momentum power :   --------------------------    
            if((do_mom == 1) .or. (do_mom==2))then   
                ! 2.1. Read data:
                Input_dirF =trim(adjustl(Input_mock_dir_pv))//trim(adjustl(file_Ino)) 
                open (unit=1675,file=Input_dirF,form='formatted')
                read(1675,*) varch ,ndatapv
                allocate(longitudepv(ndatapv));allocate(latitudepv(ndatapv));allocate(Vcmbpv(ndatapv));
                allocate(logd(ndatapv));allocate(nbmom(ndatapv));allocate(wmom(ndatapv));allocate(wrho(ndatapv))
                DO i=1,ndatapv
                    read(1675,*)longitudepv(i),latitudepv(i),Vcmbpv(i),logd(i),nbmom(i),wmom(i),wrho(i)
                ENDDO
                close(1675)                   
                ! 2.2. Calculate distance:
                Output= trim(adjustl(output_mock_dir))//  &
                        trim(adjustl('mom'))//trim(adjustl(file_Ino)) 
                allocate(xpv(ndatapv)); allocate(ypv(ndatapv)) ; allocate(zpv(ndatapv));allocate(diszpv(ndatapv));
                diszpv=0. 
                do i=1,ndatapv
                    call CosmDist(100,OmegaM,OmegaA,Hub,Cl,Vcmbpv(i)/cl,Diszpv(i))
                enddo         
                Xpv   =diszpv*cos(latitudepv/180.*mathpi)*cos(longitudepv/180.*mathpi)
                Ypv   =diszpv*cos(latitudepv/180.*mathpi)*sin(longitudepv/180.*mathpi)
                Zpv   =diszpv*sin(latitudepv/180.*mathpi)    
                allocate(Vmod(ndatapv));allocate(Vpec(ndatapv))
                deccel=3.0*OmegaM/2.0 - 1.0
                vmod  = Vcmbpv*(1.0+0.5*(1.0-deccel)*Vcmbpv/cl-(2.0-deccel-3.0*deccel*deccel)*Vcmbpv/cl*Vcmbpv/cl/6.0)
                vpec  =  log(10.)*Vmod/(1.0+Vmod/cl) * logd 
                ! 2.3. Calculate momentum power: 
                if(do_mom==1)then 
                    allocate(NAN(ndata)) 
                    call cosmic_PS(do_mom ,obs_pos, PS_Multi,output,   &
                                           ndatapv,nbmom,Xpv,Ypv,Zpv,Vpec,wmom,wrho,&                                                                     
                                           tmp ,  NAN,  NAN,  NAN,  NAN,  NAN,   &
                                           tmp,  NAN,    NAN, NAN, NAN, NAN )                                                                                                               
                    deallocate(NAN)                                                  
                    deallocate(longitudepv);deallocate(latitudepv);deallocate(Vcmbpv);deallocate(wrho)
                    deallocate(logd); deallocate(nbmom);deallocate(Vpec);deallocate(diszpv);
                    deallocate(wmom); deallocate(xpv); deallocate(ypv) ; deallocate(zpv);deallocate(Vmod); 
                endif                                          
            endif      
    
            ! cross power :   -------------------------- 
            if(do_mom==2)then  
                Output=trim(adjustl(output_mock_dir))//  &
                        trim(adjustl('crs'))//trim(adjustl(file_Ino)) 
                if(Nrand == 1)then        
                    call cosmic_PS(do_mom ,obs_pos, PS_Multi,output,   &
                                             ndatapv,nbmom,Xpv,Ypv,Zpv,Vpec,wmom, wrho,&                                                                     
                                               ndata ,  nbden,  X,  Y,  Z,  wden ,  &
                                               ndataR,  nbr,    Xr, Yr, Zr, wr )      
                else
                    call cosmic_PS(do_mom ,obs_pos, PS_Multi,output,  &
                                             ndatapv,nbmom,Xpv,Ypv,Zpv,Vpec,wmom, wrho,&                                                                     
                                               ndata ,  nbden,  X,  Y,  Z,  wden ,  &
                                               ndataRmk,  nbrmk,    Xrmk, Yrmk, Zrmk, wrmk )   
                    deallocate(nbrmk);deallocate(Xrmk);deallocate(Yrmk);deallocate(Zrmk);deallocate(wrmk)                                          
                endif                                                                                  
                deallocate(longitudepv);deallocate(latitudepv);deallocate(Vcmbpv);deallocate(wrho)
                deallocate(logd); deallocate(nbmom);deallocate(Vpec);deallocate(diszpv);
                deallocate(wmom); deallocate(xpv); deallocate(ypv) ; deallocate(zpv);deallocate(Vmod); 
                deallocate(longitude);deallocate(latitude);deallocate(Vcmb); deallocate(nbden); 
                deallocate(x); deallocate(y) ; deallocate(z);deallocate(disz);deallocate(wden);                                                                                                                                                       
            endif    
         
        Enddo             
    ENDIF 



ENDDO 










      
end program denmomPS

!--------------------------    subroutine.3    ----------------------------
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
 
            
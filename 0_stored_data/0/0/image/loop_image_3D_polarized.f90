!***********************************************************************
! Global parameters; FDUMP must be taken from the simulation 
!***********************************************************************

MODULE GLOBALS

IMPLICIT NONE

INTEGER, PUBLIC				:: FDUMP = 1      ! Frequency of dumps (see mod_input)
INTEGER, PUBLIC				:: NALPHA = 10		! Number of points in viewing angle
INTEGER, PUBLIC				:: NPHASE = 24		! Number of points in phase
DOUBLE PRECISION, PUBLIC    :: pi=ACOS(-1.0d0)
INTEGER, PUBLIC				:: NTHR = 1		! Number of threads

! Boundaries of grid in viewing angle
DOUBLE PRECISION, PUBLIC    :: amin=0.0
DOUBLE PRECISION, PUBLIC    :: amax=0.5*ACOS(-1.0d0)

! Boundaries of grid in viewing angle
DOUBLE PRECISION, PUBLIC    :: pmin=0.0
DOUBLE PRECISION, PUBLIC    :: pmax=2.0*ACOS(-1.0d0)

! Resolution of observer's screen in X and Y
INTEGER, PUBLIC				:: NX = 300		
INTEGER, PUBLIC				:: NY = 300

! Resolution of observer's screen in X and Y for the polarization
INTEGER, PUBLIC				:: NXp = 300	
INTEGER, PUBLIC				:: NYp = 300


! Boundaries of observers's screen in X and Y
DOUBLE PRECISION, PUBLIC    :: xmax=10.0
DOUBLE PRECISION, PUBLIC    :: ymax=10.0

! Boundaries of frquency bands, in units of B0
DOUBLE PRECISION, PUBLIC    :: nu1=0.1
DOUBLE PRECISION, PUBLIC    :: nu2=10.0

END MODULE GLOBALS

!***********************************************************************
! Main program, with parallelized loop on all output geokerr files
!***********************************************************************

PROGRAM LOOP_IMAGE

USE GLOBALS
USE OMP_LIB

IMPLICIT NONE

INTEGER								:: NFILE,N0,i0,i2,ix,iy,ia,ip,thread_id,ixp,iyp
CHARACTER(LEN=80)					:: CNFILE,CN0,chain_in,chain_out,ci,CNX
DOUBLE PRECISION					:: start0,finish0,delta
DOUBLE PRECISION, ALLOCATABLE		:: im_E_cs_1(:,:,:,:),im_E_cs_2(:,:,:,:),im_E_cs_3(:,:,:,:)
DOUBLE PRECISION, ALLOCATABLE		:: im_E_pc_1(:,:,:,:),im_E_pc_2(:,:,:,:),im_E_pc_3(:,:,:,:)
DOUBLE PRECISION, ALLOCATABLE		:: im_Q_cs_1(:,:,:,:),im_Q_cs_2(:,:,:,:),im_Q_cs_3(:,:,:,:)
DOUBLE PRECISION, ALLOCATABLE		:: im_Q_pc_1(:,:,:,:),im_Q_pc_2(:,:,:,:),im_Q_pc_3(:,:,:,:)
DOUBLE PRECISION, ALLOCATABLE		:: im_U_cs_1(:,:,:,:),im_U_cs_2(:,:,:,:),im_U_cs_3(:,:,:,:)
DOUBLE PRECISION, ALLOCATABLE		:: im_U_pc_1(:,:,:,:),im_U_pc_2(:,:,:,:),im_U_pc_3(:,:,:,:)


DOUBLE PRECISION, ALLOCATABLE		:: local_cs_1(:,:,:,:),local_cs_2(:,:,:,:),local_cs_3(:,:,:,:)
DOUBLE PRECISION, ALLOCATABLE		:: local_pc_1(:,:,:,:),local_pc_2(:,:,:,:),local_pc_3(:,:,:,:)
DOUBLE PRECISION, ALLOCATABLE		:: local_Q_cs_1(:,:,:,:),local_Q_cs_2(:,:,:,:),local_Q_cs_3(:,:,:,:)
DOUBLE PRECISION, ALLOCATABLE		:: local_Q_pc_1(:,:,:,:),local_Q_pc_2(:,:,:,:),local_Q_pc_3(:,:,:,:)
DOUBLE PRECISION, ALLOCATABLE		:: local_U_cs_1(:,:,:,:),local_U_cs_2(:,:,:,:),local_U_cs_3(:,:,:,:)
DOUBLE PRECISION, ALLOCATABLE		:: local_U_pc_1(:,:,:,:),local_U_pc_2(:,:,:,:),local_U_pc_3(:,:,:,:)


DOUBLE PRECISION 					:: a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,dt



!***********************************************************************

ALLOCATE(im_E_cs_1(NX,NY,NALPHA,NPHASE),im_E_cs_2(NX,NY,NALPHA,NPHASE),im_E_cs_3(NX,NY,NALPHA,NPHASE))
ALLOCATE(im_E_pc_1(NX,NY,NALPHA,NPHASE),im_E_pc_2(NX,NY,NALPHA,NPHASE),im_E_pc_3(NX,NY,NALPHA,NPHASE))
ALLOCATE(local_cs_1(NX,NY,NALPHA,NPHASE),local_cs_2(NX,NY,NALPHA,NPHASE),local_cs_3(NX,NY,NALPHA,NPHASE))
ALLOCATE(local_pc_1(NX,NY,NALPHA,NPHASE),local_pc_2(NX,NY,NALPHA,NPHASE),local_pc_3(NX,NY,NALPHA,NPHASE))

ALLOCATE(im_Q_cs_1(NXp,NYp,NALPHA,NPHASE),im_Q_cs_2(NXp,NYp,NALPHA,NPHASE),im_Q_cs_3(NXp,NYp,NALPHA,NPHASE))
ALLOCATE(im_Q_pc_1(NXp,NYp,NALPHA,NPHASE),im_Q_pc_2(NXp,NYp,NALPHA,NPHASE),im_Q_pc_3(NXp,NYp,NALPHA,NPHASE))
ALLOCATE(local_Q_cs_1(NXp,NYp,NALPHA,NPHASE),local_Q_cs_2(NXp,NYp,NALPHA,NPHASE),local_Q_cs_3(NXp,NYp,NALPHA,NPHASE))
ALLOCATE(local_Q_pc_1(NXp,NYp,NALPHA,NPHASE),local_Q_pc_2(NXp,NYp,NALPHA,NPHASE),local_Q_pc_3(NXp,NYp,NALPHA,NPHASE))
ALLOCATE(im_U_cs_1(NXp,NYp,NALPHA,NPHASE),im_U_cs_2(NXp,NYp,NALPHA,NPHASE),im_U_cs_3(NXp,NYp,NALPHA,NPHASE))
ALLOCATE(im_U_pc_1(NXp,NYp,NALPHA,NPHASE),im_U_pc_2(NXp,NYp,NALPHA,NPHASE),im_U_pc_3(NXp,NYp,NALPHA,NPHASE))
ALLOCATE(local_U_cs_1(NXp,NYp,NALPHA,NPHASE),local_U_cs_2(NXp,NYp,NALPHA,NPHASE),local_U_cs_3(NXp,NYp,NALPHA,NPHASE))
ALLOCATE(local_U_pc_1(NXp,NYp,NALPHA,NPHASE),local_U_pc_2(NXp,NYp,NALPHA,NPHASE),local_U_pc_3(NXp,NYp,NALPHA,NPHASE))





! Reads dt from input_params
OPEN(9,FILE="../../../data/input_params.dat", ACTION="read")
READ(9,*)
READ(9,*) a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,dt 
CLOSE(9)

WRITE(CNX,'(i10)') NX
CNX=adjustl(CNX)

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(FDUMP,dt,im_E_cs_1,im_E_cs_2,im_E_cs_3,im_E_pc_1,im_E_pc_2,im_E_pc_3,&
!$OMP im_Q_cs_1,im_Q_cs_2,im_Q_cs_3,im_Q_pc_1,im_Q_pc_2,im_Q_pc_3,&
!$OMP im_U_cs_1,im_U_cs_2,im_U_cs_3,im_U_pc_1,im_U_pc_2,im_U_pc_3, delta)

thread_id = OMP_GET_THREAD_NUM()

CALL GET_COMMAND_ARGUMENT(1,CNFILE) ! Total number of files to process
CALL GET_COMMAND_ARGUMENT(2,CN0) ! Number of the first file to process

! This is needed to perform a loop from 1 to NFILE, and to retrieve the time step of the initial file

READ(CNFILE,'(i6)') NFILE
READ(CN0,'(i6)') N0

im_E_cs_1=0d0; im_E_cs_2=0d0; im_E_cs_3=0d0; im_E_pc_1=0d0; im_E_pc_2=0d0; im_E_pc_3=0d0
local_cs_1=0d0; local_cs_2=0d0; local_cs_3=0d0; local_pc_1=0d0; local_pc_2=0d0; local_pc_3=0d0

im_Q_cs_1=0d0; im_Q_cs_2=0d0; im_Q_cs_3=0d0; im_Q_pc_1=0d0; im_Q_pc_2=0d0; im_Q_pc_3=0d0
local_Q_cs_1=0d0; local_Q_cs_2=0d0; local_Q_cs_3=0d0; local_Q_pc_1=0d0; local_Q_pc_2=0d0; local_Q_pc_3=0d0

im_U_cs_1=0d0; im_U_cs_2=0d0; im_U_cs_3=0d0; im_U_pc_1=0d0; im_U_pc_2=0d0; im_U_pc_3=0d0
local_U_cs_1=0d0; local_U_cs_2=0d0; local_U_cs_3=0d0; local_U_pc_1=0d0; local_U_pc_2=0d0; local_U_pc_3=0d0


!***********************************************************************

start0=OMP_GET_WTIME()

!$OMP DO

DO i0=1,NFILE

	i2=N0+(i0-1)*FDUMP
	WRITE(ci,'(i10.6)') i2
	ci=adjustl(ci)

	chain_in="../../../data/geodesics_sync/geodesics_input_"//trim(ci)//".dat"
	chain_out="../../../data/geodesics_sync/geodesics_output_"//trim(ci)//".dat"

	CALL IMAGE_GEODESICS_PARALLEL(chain_in,chain_out,NFILE, &
		local_pc_1,local_pc_2,local_pc_3,local_cs_1,local_cs_2,local_cs_3, &
		local_Q_pc_1,local_Q_pc_2,local_Q_pc_3,local_Q_cs_1,local_Q_cs_2,local_Q_cs_3, &
		local_U_pc_1,local_U_pc_2,local_U_pc_3,local_U_cs_1,local_U_cs_2,local_U_cs_3, dt,thread_id)

!	write(*,*) 'Step',i0,NFILE,thread_id

ENDDO

!$OMP END DO

!$OMP CRITICAL

im_E_cs_1=im_E_cs_1+local_cs_1
im_E_cs_2=im_E_cs_2+local_cs_2
im_E_cs_3=im_E_cs_3+local_cs_3
im_E_pc_1=im_E_pc_1+local_pc_1
im_E_pc_2=im_E_pc_2+local_pc_2
im_E_pc_3=im_E_pc_3+local_pc_3

im_Q_cs_1=im_Q_cs_1+local_Q_cs_1
im_Q_cs_2=im_Q_cs_2+local_Q_cs_2
im_Q_cs_3=im_Q_cs_3+local_Q_cs_3
im_Q_pc_1=im_Q_pc_1+local_Q_pc_1
im_Q_pc_2=im_Q_pc_2+local_Q_pc_2
im_Q_pc_3=im_Q_pc_3+local_Q_pc_3

im_U_cs_1=im_U_cs_1+local_U_cs_1
im_U_cs_2=im_U_cs_2+local_U_cs_2
im_U_cs_3=im_U_cs_3+local_U_cs_3
im_U_pc_1=im_U_pc_1+local_U_pc_1
im_U_pc_2=im_U_pc_2+local_U_pc_2
im_U_pc_3=im_U_pc_3+local_U_pc_3

!$OMP END CRITICAL

finish0=OMP_GET_WTIME()
delta=finish0-start0

!$OMP END PARALLEL

!write(*,*) "Total time:",delta

!***********************************************************************
! Dump the total data

OPEN(9,FILE="../../../data/image/im_E_cs_1_phase.dat")
    DO ip=1,NPHASE    
		DO ia=1,NALPHA    
			DO iy=1,NY
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_E_cs_1(ix,iy,ia,ip),ix=1,NX)
			ENDDO		
		ENDDO
    ENDDO
CLOSE(9)

OPEN(9,FILE="../../../data/image/im_E_cs_2_phase.dat")
    DO ip=1,NPHASE  
		DO ia=1,NALPHA    
			DO iy=1,NY
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_E_cs_2(ix,iy,ia,ip),ix=1,NX)
			ENDDO		
		ENDDO
    ENDDO
CLOSE(9)

OPEN(9,FILE="../../../data/image/im_E_cs_3_phase.dat")
    DO ip=1,NPHASE  
		DO ia=1,NALPHA    
			DO iy=1,NY
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_E_cs_3(ix,iy,ia,ip),ix=1,NX)
			ENDDO		
		ENDDO
    ENDDO
CLOSE(9)

OPEN(9,FILE="../../../data/image/im_E_pc_1_phase.dat")
    DO ip=1,NPHASE  
		DO ia=1,NALPHA    
			DO iy=1,NY
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_E_pc_1(ix,iy,ia,ip),ix=1,NX)
			ENDDO		
		ENDDO
    ENDDO
CLOSE(9)

OPEN(9,FILE="../../../data/image/im_E_pc_2_phase.dat")
    DO ip=1,NPHASE  
		DO ia=1,NALPHA    
			DO iy=1,NY
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_E_pc_2(ix,iy,ia,ip),ix=1,NX)
			ENDDO		
		ENDDO
    ENDDO
CLOSE(9)

OPEN(9,FILE="../../../data/image/im_E_pc_3_phase.dat")
    DO ip=1,NPHASE  
		DO ia=1,NALPHA    
			DO iy=1,NY
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_E_pc_3(ix,iy,ia,ip),ix=1,NX)
			ENDDO		
		ENDDO
    ENDDO
CLOSE(9)

!**********************************8

OPEN(9,FILE="../../../data/image/im_Q_cs_1_phase.dat")
    DO ip=1,NPHASE    
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_Q_cs_1(ix,iy,ia,ip),ix=1,NXp)
			ENDDO		
		ENDDO
    ENDDO
CLOSE(9)

OPEN(9,FILE="../../../data/image/im_Q_cs_2_phase.dat")
    DO ip=1,NPHASE  
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_Q_cs_2(ix,iy,ia,ip),ix=1,NXp)
			ENDDO		
		ENDDO
    ENDDO
CLOSE(9)

OPEN(9,FILE="../../../data/image/im_Q_cs_3_phase.dat")
    DO ip=1,NPHASE  
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_Q_cs_3(ix,iy,ia,ip),ix=1,NXp)
			ENDDO		
		ENDDO
    ENDDO
CLOSE(9)

OPEN(9,FILE="../../../data/image/im_Q_pc_1_phase.dat")
    DO ip=1,NPHASE  
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_Q_pc_1(ix,iy,ia,ip),ix=1,NXp)
			ENDDO		
		ENDDO
    ENDDO
CLOSE(9)

OPEN(9,FILE="../../../data/image/im_Q_pc_2_phase.dat")
    DO ip=1,NPHASE  
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_Q_pc_2(ix,iy,ia,ip),ix=1,NXp)
			ENDDO		
		ENDDO
    ENDDO
CLOSE(9)

OPEN(9,FILE="../../../data/image/im_Q_pc_3_phase.dat")
    DO ip=1,NPHASE  
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_Q_pc_3(ix,iy,ia,ip),ix=1,NXp)
			ENDDO		
		ENDDO
    ENDDO
CLOSE(9)
!**************************************88
OPEN(9,FILE="../../../data/image/im_U_cs_1_phase.dat")
    DO ip=1,NPHASE    
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_U_cs_1(ix,iy,ia,ip),ix=1,NXp)
			ENDDO		
		ENDDO
    ENDDO
CLOSE(9)

OPEN(9,FILE="../../../data/image/im_U_cs_2_phase.dat")
    DO ip=1,NPHASE  
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_U_cs_2(ix,iy,ia,ip),ix=1,NXp)
			ENDDO		
		ENDDO
    ENDDO
CLOSE(9)

OPEN(9,FILE="../../../data/image/im_U_cs_3_phase.dat")
    DO ip=1,NPHASE  
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_U_cs_3(ix,iy,ia,ip),ix=1,NXp)
			ENDDO		
		ENDDO
    ENDDO
CLOSE(9)

OPEN(9,FILE="../../../data/image/im_U_pc_1_phase.dat")
    DO ip=1,NPHASE  
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_U_pc_1(ix,iy,ia,ip),ix=1,NXp)
			ENDDO		
		ENDDO
    ENDDO
CLOSE(9)

OPEN(9,FILE="../../../data/image/im_U_pc_2_phase.dat")
    DO ip=1,NPHASE  
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_U_pc_2(ix,iy,ia,ip),ix=1,NXp)
			ENDDO		
		ENDDO
    ENDDO
CLOSE(9)

OPEN(9,FILE="../../../data/image/im_U_pc_3_phase.dat")
    DO ip=1,NPHASE  
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_U_pc_3(ix,iy,ia,ip),ix=1,NXp)
			ENDDO		
		ENDDO
    ENDDO
CLOSE(9)





CONTAINS

!***********************************************************************
! Subroutine which produces the image data of a single output geokerr file
!***********************************************************************

SUBROUTINE IMAGE_GEODESICS_PARALLEL(chain_in,chain_out,NFILE, &
	im_E_pc_1,im_E_pc_2,im_E_pc_3,im_E_cs_1,im_E_cs_2,im_E_cs_3, &
	im_Q_pc_1,im_Q_pc_2,im_Q_pc_3,im_Q_cs_1,im_Q_cs_2,im_Q_cs_3, &
	im_U_pc_1,im_U_pc_2,im_U_pc_3,im_U_cs_1,im_U_cs_2,im_U_cs_3, dt,thread_id)

USE GLOBALS

IMPLICIT NONE

DOUBLE PRECISION 				   :: dalpha,dphase,end_bt,vx,vy,vz,alphap,t_out,t_delay,r0,th0,rh,rin
DOUBLE PRECISION 				   :: suf,rf,phif,spec,vph,vth,car,vr,cosa,sina,cosp,sinp
DOUBLE PRECISION 				   :: a,q2,l,u0,mu0,uf,muf,su,TPR,SM,wt,tbl,phibl,E,dtf,dphf,dt,dphi,phif_mod
DOUBLE PRECISION 				   :: start,finish,dx,dy,x,y,uif,kappa1,kappa2,dxp,dyp,nu,scale,fa,fb,stokesU,stokesQ
CHARACTER(LEN=60)				   :: chain_in,chain_out
INTEGER 						   :: i,j,N_UNDER,NT,NGEO,NUP,STANDARD,QL,NTIME,NFILE,OK,OK1,OK2,OK3,OK4,OK5
INTEGER							   :: ialpha,iphase,ix,iy,thread_id,unit_in,unit_out
DOUBLE PRECISION, ALLOCATABLE	   :: MUFI(:)
DOUBLE PRECISION, DIMENSION(1:NX,1:NY,1:NALPHA,1:NPHASE) :: im_E_pc_1,im_E_pc_2,im_E_pc_3,im_E_cs_1,im_E_cs_2,im_E_cs_3
DOUBLE PRECISION, DIMENSION(1:NXp,1:NYp,1:NALPHA,1:NPHASE) :: im_Q_pc_1,im_Q_pc_2,im_Q_pc_3,im_Q_cs_1,im_Q_cs_2,im_Q_cs_3
DOUBLE PRECISION, DIMENSION(1:NXp,1:NYp,1:NALPHA,1:NPHASE) :: im_U_pc_1,im_U_pc_2,im_U_pc_3,im_U_cs_1,im_U_cs_2,im_U_cs_3
dalpha=(amax-amin)/NALPHA

dx=2*xmax/NX
dy=2*ymax/NY
dxp=2*xmax/NXp
dyp=2*ymax/NYp


dphase=(pmax-pmin)/NPHASE

! These units must differ for every thread
unit_in=thread_id
unit_out=NTHR*(thread_id+1)

!!!!!!!!!!!!!!!!!!!!!
OPEN(unit_in,FILE=chain_in, ACTION="read") ! Reads input file for geokerr

READ(unit_in,*) STANDARD
READ(unit_in,*) QL
READ(unit_in,*) NGEO ! Number of geodesics
READ(unit_in,*) a ! BH spin
READ(unit_in,*) NUP ! Number of points per geodesic

NUP=3

ALLOCATE(mufi(1:NUP+1)) ! cos(theta)

rh=1+SQRT(1-a*a)
rin=6.0

OPEN(unit_out,FILE=chain_out, ACTION="read") ! Reads output file from geokerr

READ(unit_out,*) NGEO, a

OK=0
i=1

!DO i=1,NGEO-1
DO WHILE ((i.LE.(NGEO-1)).AND.(OK.EQ.0))

	READ(unit_in,*,IOSTAT=OK1) q2, l, u0, mu0, uf, su, TPR, SM, wt, tbl, phibl, E, kappa1, kappa2 ! One line per geodesic in input
	READ(unit_out,*,IOSTAT=OK2) ! Unnecessary recap line in output

	mufi(1)=mu0

	DO j=1,NUP-2
		READ(unit_out,*,IOSTAT=OK3) uf,mufi(1+j)
	ENDDO 

	READ(unit_out,*,IOSTAT=OK4) uif, mufi(NUP)

	READ(unit_out,*,IOSTAT=OK5) uf, muf, dtf, dphf ! Reads last line of output for the ith geodesic 

	OK=OK1+OK2+OK3+OK4+OK5

	mufi(NUP+1)=muf

	suf=SQRT(1-muf*muf)
	rf=1/uf	
	phif=dphf+phibl

	phif_mod=phif-FLOOR(phif/(2.*pi))*2.*pi
	
	! Checks pth^2 remains positive along geodesic, as it should
	! spec=MINVAL(q2+mufi*mufi*(a*a-q2-l*l)-(a*mufi*mufi)**2)
	! vph=l/(rf*suf)

	! r0=1/u0
	! th0=ACOS(mu0)

	! Checks pth^2 remains positive along geodesic, as it should
	spec=MINVAL(q2+mufi*mufi*(a*a-q2-l*l)-(a*mufi*mufi)**2) !TODO: SLOW

	r0=1/u0
	th0=ACOS(mu0) !TODO: SLOW

	!  Computes vth, if it can, and car needs to be <1 to compute vr
	vth=SQRT(q2+(muf*a)**2 - (muf*l/suf)**2) * SIGN(1d0,mufi(NUP)-muf) / rf 
	vph=l/(rf*suf)
	car=vth*vth+vph*vph
	vr=SQRT(1-car)
	vz=muf*vr-suf*vth
	cosa=vz
	sina=SQRT(1-cosa*cosa)
	alphap=ACOS(cosa)

	x=-l/suf
	y=vth*rf

!	vx=vr*suf*COS(phif)+vth*muf*COS(phif)-vph*SIN(phif)	
!	vy=vr*suf*SIN(phif)+vth*muf*SIN(phif)+vph*COS(phif)
!	cosp=vx/SQRT(vx*vx+vy*vy)
!	sinp=vy/SQRT(vx*vx+vy*vy)

!	alphap=ACOS(muf)
!	sina=SQRT(1-cosa*cosa)

	! Checks that geodesics is correct, and that photon actually got out of the box
	IF ((rf>5.0).AND.(spec>0.0).AND.(car<1.0).AND.(ABS(x).LT.xmax).AND.(ABS(y).LT.ymax).AND.(dtf.GT.(0d0)) &
		.AND.(uif/SQRT(1-mufi(NUP)*mufi(NUP)).GT.uf/suf).AND.(alphap.LT.amax).AND.&
		(alphap.GT.amin).AND.(phif_mod.GT.pmin).AND.(phif_mod.LT.pmax)) THEN !.AND.(E.LT.(1d4))

		ialpha=FLOOR((alphap-amin)/dalpha)+1
		iphase=FLOOR((phif_mod-pmin)/dphase)+1

		ix=FLOOR((x+xmax)/dx)+1
		iy=FLOOR((y+ymax)/dy)+1		

		ixp=FLOOR((x+xmax)/dxp)+1
		iyp=FLOOR((y+ymax)/dyp)+1	

		! compute the observed polarization
		nu=-(x + a * suf)
		scale=SQRT((kappa1*kappa1 + kappa2*kappa2)*(y*y+nu*nu))
		! scale=SQRT((y*y+nu*nu))
		fa=(y*kappa2-nu*kappa1)/scale
		fb=(y*kappa1+nu*kappa2)/scale
		
		stokesQ=-fa*fa+fb*fb 
		stokesU=-2*fa*fb 
                
		IF ((r0.LT.(0.8*rin)).AND.(r0.GT.rh).AND.(ABS(th0-0.5*pi).LT.(0.5))) THEN
!		IF ((r0.GT.rh).AND.(ABS(th0-0.5*pi).LT.(0.5))) THEN

			IF (E.LT.nu1) THEN
				im_E_cs_1(ix,iy,ialpha,iphase)=im_E_cs_1(ix,iy,ialpha,iphase)+wt/((sina*dalpha)*dx*dy*dphase)
				im_U_cs_1(ixp,iyp,ialpha,iphase)=im_U_cs_1(ixp,iyp,ialpha,iphase)+wt/((sina*dalpha)*dxp*dyp*dphase)*stokesU
				im_Q_cs_1(ixp,iyp,ialpha,iphase)=im_Q_cs_1(ixp,iyp,ialpha,iphase)+wt/((sina*dalpha)*dxp*dyp*dphase)*stokesQ

			ELSE IF (E.GT.nu2) THEN
				im_E_cs_3(ix,iy,ialpha,iphase)=im_E_cs_3(ix,iy,ialpha,iphase)+wt/((sina*dalpha)*dx*dy*dphase)
				im_U_cs_3(ixp,iyp,ialpha,iphase)=im_U_cs_3(ixp,iyp,ialpha,iphase)+wt/((sina*dalpha)*dxp*dyp*dphase)*stokesU
				im_Q_cs_3(ixp,iyp,ialpha,iphase)=im_Q_cs_3(ixp,iyp,ialpha,iphase)+wt/((sina*dalpha)*dxp*dyp*dphase)*stokesQ
			ELSE
				im_E_cs_2(ix,iy,ialpha,iphase)=im_E_cs_2(ix,iy,ialpha,iphase)+wt/((sina*dalpha)*dx*dy*dphase)
				im_U_cs_2(ixp,iyp,ialpha,iphase)=im_U_cs_2(ixp,iyp,ialpha,iphase)+wt/((sina*dalpha)*dxp*dyp*dphase)*stokesU
				im_Q_cs_2(ixp,iyp,ialpha,iphase)=im_Q_cs_2(ixp,iyp,ialpha,iphase)+wt/((sina*dalpha)*dxp*dyp*dphase)*stokesQ
			END IF

		ELSE IF ((r0.GT.rh).AND.(ABS(th0-0.5*pi).GT.(0.5))) THEN

			IF (E.LT.nu1) THEN
				im_E_pc_1(ix,iy,ialpha,iphase)=im_E_pc_1(ix,iy,ialpha,iphase)+wt/((sina*dalpha)*dx*dy*dphase)
				im_U_pc_1(ixp,iyp,ialpha,iphase)=im_U_pc_1(ixp,iyp,ialpha,iphase)+wt/((sina*dalpha)*dxp*dyp*dphase)*stokesU 
				im_Q_pc_1(ixp,iyp,ialpha,iphase)=im_Q_pc_1(ixp,iyp,ialpha,iphase)+wt/((sina*dalpha)*dxp*dyp*dphase)*stokesQ
			ELSE IF (E.GT.nu2) THEN
				im_E_pc_3(ix,iy,ialpha,iphase)=im_E_pc_3(ix,iy,ialpha,iphase)+wt/((sina*dalpha)*dx*dy*dphase)
				im_U_pc_3(ixp,iyp,ialpha,iphase)=im_U_pc_3(ixp,iyp,ialpha,iphase)+wt/((sina*dalpha)*dxp*dyp*dphase)*stokesU
				im_Q_pc_3(ixp,iyp,ialpha,iphase)=im_Q_pc_3(ixp,iyp,ialpha,iphase)+wt/((sina*dalpha)*dxp*dyp*dphase)*stokesQ
			ELSE
				im_E_pc_2(ix,iy,ialpha,iphase)=im_E_pc_2(ix,iy,ialpha,iphase)+wt/((sina*dalpha)*dx*dy*dphase)
				im_U_pc_2(ixp,iyp,ialpha,iphase)=im_U_pc_2(ixp,iyp,ialpha,iphase)+wt/((sina*dalpha)*dxp*dyp*dphase)*stokesU
				im_Q_pc_2(ixp,iyp,ialpha,iphase)=im_Q_pc_2(ixp,iyp,ialpha,iphase)+wt/((sina*dalpha)*dxp*dyp*dphase)*stokesQ
			END IF
		END IF

	END IF

i=i+1

ENDDO

! write(*,*) 'Final rank',i,NGEO,OK1,OK2,OK3,OK4,OK5

CLOSE(unit_in)
CLOSE(unit_out)

END SUBROUTINE IMAGE_GEODESICS_PARALLEL

END PROGRAM LOOP_IMAGE

!***********************************************************************
! Global parameters; FDUMP must be taken from the simulation 
!***********************************************************************


MODULE GLOBALS

#ifndef nxy
#define nxy=200
#endif

#ifndef FRAC
#define FRAC=1
#endif


#ifndef DIRECTORY
#define DIRECTORY="data"
#endif

INTEGER, PUBLIC				:: FDUMP = 1      ! Frequency of dumps (see mod_input)
INTEGER, PUBLIC				:: NALPHA = 40		! Number of points in viewing angle
DOUBLE PRECISION, PUBLIC    :: pi=ACOS(-1.0d0)
INTEGER, PUBLIC				:: NTHR = 2		! Number of threads

! Boundaries of grid in viewing angle
DOUBLE PRECISION, PUBLIC    :: amin=0.0
DOUBLE PRECISION, PUBLIC    :: amax=0.5*ACOS(-1.0d0)
DOUBLE PRECISION, PUBLIC    :: fraction=FRAC


! Resolution of observer's screen in X and Y
INTEGER, PUBLIC				:: NX = nxy		
INTEGER, PUBLIC				:: NY = nxy

! Resolution of observer's screen in X and Y for the polarization
INTEGER, PUBLIC				:: NXp = nxy	
INTEGER, PUBLIC				:: NYp = nxy

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

INTEGER								:: NFILE,N0,i0,i2,ix,iy,ia,thread_id,log,ixp,iyp
CHARACTER(LEN=80)					:: CNFILE,CN0,chain_in,chain_out,ci,CNX,CNXp
DOUBLE PRECISION					:: start0,finish0,delta
DOUBLE PRECISION, ALLOCATABLE		:: im_E_cs_1(:,:,:),im_E_cs_2(:,:,:),im_E_cs_3(:,:,:)
DOUBLE PRECISION, ALLOCATABLE		:: im_E_pc_1(:,:,:),im_E_pc_2(:,:,:),im_E_pc_3(:,:,:)
DOUBLE PRECISION, ALLOCATABLE		:: im_Q_cs_1(:,:,:),im_Q_cs_2(:,:,:),im_Q_cs_3(:,:,:)
DOUBLE PRECISION, ALLOCATABLE		:: im_Q_pc_1(:,:,:),im_Q_pc_2(:,:,:),im_Q_pc_3(:,:,:)
DOUBLE PRECISION, ALLOCATABLE		:: im_U_cs_1(:,:,:),im_U_cs_2(:,:,:),im_U_cs_3(:,:,:)
DOUBLE PRECISION, ALLOCATABLE		:: im_U_pc_1(:,:,:),im_U_pc_2(:,:,:),im_U_pc_3(:,:,:)


DOUBLE PRECISION 					:: a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,dt

!***********************************************************************

ALLOCATE(im_E_cs_1(NX,NY,NALPHA),im_E_cs_2(NX,NY,NALPHA),im_E_cs_3(NX,NY,NALPHA))
ALLOCATE(im_E_pc_1(NX,NY,NALPHA),im_E_pc_2(NX,NY,NALPHA),im_E_pc_3(NX,NY,NALPHA))


ALLOCATE(im_Q_cs_1(NXp,NYp,NALPHA),im_Q_cs_2(NXp,NYp,NALPHA),im_Q_cs_3(NXp,NYp,NALPHA))
ALLOCATE(im_Q_pc_1(NXp,NYp,NALPHA),im_Q_pc_2(NXp,NYp,NALPHA),im_Q_pc_3(NXp,NYp,NALPHA))

ALLOCATE(im_U_cs_1(NXp,NYp,NALPHA),im_U_cs_2(NXp,NYp,NALPHA),im_U_cs_3(NXp,NYp,NALPHA))
ALLOCATE(im_U_pc_1(NXp,NYp,NALPHA),im_U_pc_2(NXp,NYp,NALPHA),im_U_pc_3(NXp,NYp,NALPHA))


! Reads dt from input_params
OPEN(9,FILE="../../../"//DIRECTORY//"/input_params.dat", ACTION="read")
READ(9,*)
READ(9,*) a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,dt 
CLOSE(9)

WRITE(CNX,'(i10)') NX
CNX=adjustl(CNX)

WRITE(CNXp,'(i10)') NXp
CNXp=adjustl(CNXp)


thread_id = OMP_GET_THREAD_NUM()



CALL GET_COMMAND_ARGUMENT(1,CNFILE) ! Total number of files to process
CALL GET_COMMAND_ARGUMENT(2,CN0) ! Number of the first file to process

! This is needed to perform a loop from 1 to NFILE, and to retrieve the time step of the initial file

READ(CNFILE,'(i6)') NFILE
READ(CN0,'(i6)') N0

im_E_cs_1=0d0; im_E_cs_2=0d0; im_E_cs_3=0d0; im_E_pc_1=0d0; im_E_pc_2=0d0; im_E_pc_3=0d0

im_Q_cs_1=0d0; im_Q_cs_2=0d0; im_Q_cs_3=0d0; im_Q_pc_1=0d0; im_Q_pc_2=0d0; im_Q_pc_3=0d0

im_U_cs_1=0d0; im_U_cs_2=0d0; im_U_cs_3=0d0; im_U_pc_1=0d0; im_U_pc_2=0d0; im_U_pc_3=0d0
!***********************************************************************



DO i0=1,NFILE



	i2=N0+(i0-1)*FDUMP
	WRITE(ci,'(i10.6)') i2
	ci=adjustl(ci)

	chain_in="../../../"//DIRECTORY//"/geodesics_sync/geodesics_input_"//trim(ci)//".dat"
	chain_out="../../../"//DIRECTORY//"/geodesics_sync/geodesics_output_"//trim(ci)//".dat"


	! PRINT *, "../../../data/geodesics_sync/geodesics_input_"//trim(ci)//".dat"

	CALL IMAGE_GEODESICS_PARALLEL(chain_in,chain_out,NFILE, &
	im_E_pc_1,im_E_pc_2,im_E_pc_3,im_E_cs_1,im_E_cs_2,im_E_cs_3, &
	im_Q_pc_1,im_Q_pc_2,im_Q_pc_3,im_Q_cs_1,im_Q_cs_2,im_Q_cs_3, &
	im_U_pc_1,im_U_pc_2,im_U_pc_3,im_U_cs_1,im_U_cs_2,im_U_cs_3, dt)

!	write(*,*) 'Step',i0,NFILE,thread_id

ENDDO






! im_E_cs_1=im_E_cs_1+local_cs_1
! im_E_cs_2=im_E_cs_2+local_cs_2
! im_E_cs_3=im_E_cs_3+local_cs_3
! im_E_pc_1=im_E_pc_1+local_pc_1
! im_E_pc_2=im_E_pc_2+local_pc_2
! im_E_pc_3=im_E_pc_3+local_pc_3

! im_Q_cs_1=im_Q_cs_1+local_Q_cs_1
! im_Q_cs_2=im_Q_cs_2+local_Q_cs_2
! im_Q_cs_3=im_Q_cs_3+local_Q_cs_3
! im_Q_pc_1=im_Q_pc_1+local_Q_pc_1
! im_Q_pc_2=im_Q_pc_2+local_Q_pc_2
! im_Q_pc_3=im_Q_pc_3+local_Q_pc_3

! im_U_cs_1=im_U_cs_1+local_U_cs_1
! im_U_cs_2=im_U_cs_2+local_U_cs_2
! im_U_cs_3=im_U_cs_3+local_U_cs_3
! im_U_pc_1=im_U_pc_1+local_U_pc_1
! im_U_pc_2=im_U_pc_2+local_U_pc_2
! im_U_pc_3=im_U_pc_3+local_U_pc_3



! finish0=OMP_GET_WTIME()

! delta=finish0-start0


! OPEN(log, FILE="log2.log", ACTION="write")
! WRITE(log, *) im_Q_cs_2
!write(*,*) "Total time:",delta

!***********************************************************************
! Dump the total data


OPEN(9,FILE="../../../"//DIRECTORY//"/image/im_E_cs_1.dat")
		DO ia=1,NALPHA    
			DO iy=1,NY
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_E_cs_1(ix,iy,ia),ix=1,NX)
			ENDDO		
		ENDDO
CLOSE(9)

OPEN(9,FILE="../../../"//DIRECTORY//"/image/im_E_cs_2.dat")
		DO ia=1,NALPHA    
			DO iy=1,NY
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_E_cs_2(ix,iy,ia),ix=1,NX)
			ENDDO		
		ENDDO
CLOSE(9)

OPEN(9,FILE="../../../"//DIRECTORY//"/image/im_E_cs_3.dat")
		DO ia=1,NALPHA    
			DO iy=1,NY
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_E_cs_3(ix,iy,ia),ix=1,NX)
			ENDDO		
		ENDDO
CLOSE(9)

OPEN(9,FILE="../../../"//DIRECTORY//"/image/im_E_pc_1.dat")
		DO ia=1,NALPHA    
			DO iy=1,NY
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_E_pc_1(ix,iy,ia),ix=1,NX)
			ENDDO		
		ENDDO
CLOSE(9)

OPEN(9,FILE="../../../"//DIRECTORY//"/image/im_E_pc_2.dat")
		DO ia=1,NALPHA    
			DO iy=1,NY
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_E_pc_2(ix,iy,ia),ix=1,NX)
			ENDDO		
		ENDDO
CLOSE(9)

OPEN(9,FILE="../../../"//DIRECTORY//"/image/im_E_pc_3.dat")
		DO ia=1,NALPHA    
			DO iy=1,NY
	   			WRITE(9,'(' // trim(CNX) // 'E26.16E3)') (im_E_pc_3(ix,iy,ia),ix=1,NX)
			ENDDO		
		ENDDO
CLOSE(9)

!***********
OPEN(9,FILE="../../../"//DIRECTORY//"/image/im_Q_cs_1.dat")
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNXp) // 'E26.16E3)') (im_Q_cs_1(ix,iy,ia),ix=1,NXp)
			ENDDO		
		ENDDO
CLOSE(9)

OPEN(9,FILE="../../../"//DIRECTORY//"/image/im_Q_cs_2.dat")
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNXp) // 'E26.16E3)') (im_Q_cs_2(ix,iy,ia),ix=1,NXp)
			ENDDO		
		ENDDO
CLOSE(9)

OPEN(9,FILE="../../../"//DIRECTORY//"/image/im_Q_cs_3.dat")
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNXp) // 'E26.16E3)') (im_Q_cs_3(ix,iy,ia),ix=1,NXp)
			ENDDO		
		ENDDO
CLOSE(9)

OPEN(9,FILE="../../../"//DIRECTORY//"/image/im_Q_pc_1.dat")
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNXp) // 'E26.16E3)') (im_Q_pc_1(ix,iy,ia),ix=1,NXp)
			ENDDO		
		ENDDO
CLOSE(9)

OPEN(9,FILE="../../../"//DIRECTORY//"/image/im_Q_pc_2.dat")
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNXp) // 'E26.16E3)') (im_Q_pc_2(ix,iy,ia),ix=1,NXp)
			ENDDO		
		ENDDO
CLOSE(9)

OPEN(9,FILE="../../../"//DIRECTORY//"/image/im_Q_pc_3.dat")
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNXp) // 'E26.16E3)') (im_Q_pc_3(ix,iy,ia),ix=1,NXp)
			ENDDO		
		ENDDO
CLOSE(9)
!***********
OPEN(9,FILE="../../../"//DIRECTORY//"/image/im_U_cs_1.dat")
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNXp) // 'E26.16E3)') (im_U_cs_1(ix,iy,ia),ix=1,NXp)
			ENDDO		
		ENDDO
CLOSE(9)

OPEN(9,FILE="../../../"//DIRECTORY//"/image/im_U_cs_2.dat")
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNXp) // 'E26.16E3)') (im_U_cs_2(ix,iy,ia),ix=1,NXp)
			ENDDO		
		ENDDO
CLOSE(9)

OPEN(9,FILE="../../../"//DIRECTORY//"/image/im_U_cs_3.dat")
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNXp) // 'E26.16E3)') (im_U_cs_3(ix,iy,ia),ix=1,NXp)
			ENDDO		
		ENDDO
CLOSE(9)

OPEN(9,FILE="../../../"//DIRECTORY//"/image/im_U_pc_1.dat")
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNXp) // 'E26.16E3)') (im_U_pc_1(ix,iy,ia),ix=1,NXp)
			ENDDO		
		ENDDO
CLOSE(9)

OPEN(9,FILE="../../../"//DIRECTORY//"/image/im_U_pc_2.dat")
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNXp) // 'E26.16E3)') (im_U_pc_2(ix,iy,ia),ix=1,NXp)
			ENDDO		
		ENDDO
CLOSE(9)

OPEN(9,FILE="../../../"//DIRECTORY//"/image/im_U_pc_3.dat")
		DO ia=1,NALPHA    
			DO iy=1,NYp
	   			WRITE(9,'(' // trim(CNXp) // 'E26.16E3)') (im_U_pc_3(ix,iy,ia),ix=1,NXp)
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
im_U_pc_1,im_U_pc_2,im_U_pc_3,im_U_cs_1,im_U_cs_2,im_U_cs_3, dt)

USE GLOBALS

IMPLICIT NONE

DOUBLE PRECISION 				   :: dalpha,end_bt,vx,vy,vz,alphap,t_out,t_delay,r0,th0,rh,rand
DOUBLE PRECISION 				   :: suf,rf,phif,spec,vph,vth,car,vr,cosa,sina,cosp,sinp
DOUBLE PRECISION 				   :: a,q2,l,u0,mu0,uf,muf,su,TPR,SM,wt,tbl,phibl,E,dtf,dphf,dt,kappa1,kappa2
DOUBLE PRECISION 				   :: start,finish,dx,dy,x,y,uif,fa,fb,stokesU,stokesQ,nu,scale,chi,stokesI,dxp,dyp
CHARACTER(LEN=70)				   :: chain_in,chain_out
INTEGER 						   :: i,j,N_UNDER,NT,NGEO,NUP,STANDARD,QL,NTIME,NFILE,OK,OK1,OK2,OK3,OK4,OK5
INTEGER							   :: ialpha,ix,iy,thread_id,unit_in,unit_out,log
DOUBLE PRECISION, ALLOCATABLE	   :: MUFI(:)
DOUBLE PRECISION, DIMENSION(1:NX,1:NY,1:NALPHA) :: im_E_pc_1,im_E_pc_2,im_E_pc_3,im_E_cs_1,im_E_cs_2,im_E_cs_3
DOUBLE PRECISION, DIMENSION(1:NXp,1:NYp,1:NALPHA) :: im_Q_pc_1,im_Q_pc_2,im_Q_pc_3,im_Q_cs_1,im_Q_cs_2,im_Q_cs_3
DOUBLE PRECISION, DIMENSION(1:NXp,1:NYp,1:NALPHA) :: im_U_pc_1,im_U_pc_2,im_U_pc_3,im_U_cs_1,im_U_cs_2,im_U_cs_3
DOUBLE PRECISION, ALLOCATABLE		:: local_cs_1(:,:,:),local_cs_2(:,:,:),local_cs_3(:,:,:)
DOUBLE PRECISION, ALLOCATABLE		:: local_pc_1(:,:,:),local_pc_2(:,:,:),local_pc_3(:,:,:)

DOUBLE PRECISION, ALLOCATABLE		:: local_Q_cs_1(:,:,:),local_Q_cs_2(:,:,:),local_Q_cs_3(:,:,:)
DOUBLE PRECISION, ALLOCATABLE		:: local_Q_pc_1(:,:,:),local_Q_pc_2(:,:,:),local_Q_pc_3(:,:,:)
DOUBLE PRECISION, ALLOCATABLE		:: local_U_cs_1(:,:,:),local_U_cs_2(:,:,:),local_U_cs_3(:,:,:)
DOUBLE PRECISION, ALLOCATABLE		:: local_U_pc_1(:,:,:),local_U_pc_2(:,:,:),local_U_pc_3(:,:,:)


dalpha=(amax-amin)/NALPHA

dx=2*xmax/NX
dy=2*ymax/NY

dxp=2*xmax/NXp
dyp=2*ymax/NYp
thread_id = OMP_GET_THREAD_NUM()

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

OPEN(unit_out,FILE=chain_out, ACTION="read") ! Reads output file from geokerr
READ(unit_out,*) NGEO, a

OK=0

ALLOCATE(local_cs_1(NX,NY,NALPHA),local_cs_2(NX,NY,NALPHA),local_cs_3(NX,NY,NALPHA))
ALLOCATE(local_pc_1(NX,NY,NALPHA),local_pc_2(NX,NY,NALPHA),local_pc_3(NX,NY,NALPHA))
ALLOCATE(local_Q_cs_1(NXp,NYp,NALPHA),local_Q_cs_2(NXp,NYp,NALPHA),local_Q_cs_3(NXp,NYp,NALPHA))
ALLOCATE(local_Q_pc_1(NXp,NYp,NALPHA),local_Q_pc_2(NXp,NYp,NALPHA),local_Q_pc_3(NXp,NYp,NALPHA))
ALLOCATE(local_U_cs_1(NXp,NYp,NALPHA),local_U_cs_2(NXp,NYp,NALPHA),local_U_cs_3(NXp,NYp,NALPHA))
ALLOCATE(local_U_pc_1(NXp,NYp,NALPHA),local_U_pc_2(NXp,NYp,NALPHA),local_U_pc_3(NXp,NYp,NALPHA))

!!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(dalpha, dx, dy, dxp, dyp, unit_in, &
!!$OMP unit_out, STANDARD, QL, NGEO, a, NUP, rh, &
!!$OMP im_E_cs_1,im_E_cs_2,im_E_cs_3,im_E_pc_1,im_E_pc_2,im_E_pc_3,&
!!$OMP im_Q_cs_1,im_Q_cs_2,im_Q_cs_3,im_Q_pc_1,im_Q_pc_2,im_Q_pc_3,&
!!$OMP im_U_cs_1,im_U_cs_2,im_U_cs_3,im_U_pc_1,im_U_pc_2,im_U_pc_3)
local_cs_1=0d0; local_cs_2=0d0; local_cs_3=0d0; local_pc_1=0d0; local_pc_2=0d0; local_pc_3=0d0
local_Q_cs_1=0d0; local_Q_cs_2=0d0; local_Q_cs_3=0d0; local_Q_pc_1=0d0; local_Q_pc_2=0d0; local_Q_pc_3=0d0
local_U_cs_1=0d0; local_U_cs_2=0d0; local_U_cs_3=0d0; local_U_pc_1=0d0; local_U_pc_2=0d0; local_U_pc_3=0d0

CALL RANDOM_SEED()

!!$OMP DO
DO i=1, (NGEO-1)
! DO WHILE ((i.LE.(NGEO-1)).AND.(OK.EQ.0))

!!$OMP CRITICAL
	READ(unit_in,*,IOSTAT=OK1) q2, l, u0, mu0, uf, su, TPR, SM, wt, tbl, phibl, E, kappa1, kappa2 ! One line per geodesic in input
	READ(unit_out,*,IOSTAT=OK2) ! Unnecessary recap line in output

	mufi(1)=mu0

	DO j=1,NUP-2
		READ(unit_out,*,IOSTAT=OK3) uf,mufi(1+j)
	ENDDO 

	READ(unit_out,*,IOSTAT=OK4) uif, mufi(NUP)

	READ(unit_out,*,IOSTAT=OK5) uf, muf, dtf, dphf ! Reads last line of output for the ith geodesic 
!!$OMP END CRITICAL

	OK=OK1+OK2+OK3+OK4+OK5

	mufi(NUP+1)=muf

	suf=SQRT(1-muf*muf)
	rf=1/uf	
	phif=dphf+phibl

	r0=1/u0
	th0=ACOS(mu0)

	! Checks pth^2 remains positive along geodesic, as it should
	spec=MINVAL(q2+mufi*mufi*(a*a-q2-l*l)-(a*mufi*mufi)**2)

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

	rand = 0.0
	call random_number(rand)
	! Checks that geodesics is correct, and that photon actually got out of the box
	IF ((rf>5.0).AND.(spec>0.0).AND.(car<1.0).AND.(ABS(x).LT.xmax).AND.(ABS(y).LT.ymax) &
		.AND.(dtf.GT.(0d0)).AND.(alphap.LT.amax).AND.(alphap.GT.amin).AND.(rand.LT.fraction)) THEN
!		.AND.(uif/SQRT(1-mufi(NUP)*mufi(NUP)).GT.uf/suf) .AND.(E.LT.(1d4)) THEN

		ialpha=FLOOR((alphap-amin)/dalpha)+1

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
		
		! These both should be negative what they are 
		! thsy should be ie Q = fb^2 - fa^2, U = -2fa*fb
		stokesQ=fb*fb-fa*fa
		stokesU=-2*fa*fb 
		
		

		IF ((r0.GT.rh).AND.(ABS(th0-0.5*pi).LT.(0.5))) THEN

			IF (E.LT.(nu1)) THEN
				local_cs_1(ix,iy,ialpha)=local_cs_1(ix,iy,ialpha)+wt/((2*pi*sina*dalpha)*dx*dy)
				local_Q_cs_1(ixp,iyp,ialpha)=local_Q_cs_1(ixp,iyp,ialpha)+wt/((2*pi*sina*dalpha)*dxp*dyp)*stokesQ
				local_U_cs_1(ixp,iyp,ialpha)=local_U_cs_1(ixp,iyp,ialpha)+wt/((2*pi*sina*dalpha)*dxp*dyp)*stokesU

			ELSE IF (E.GT.(nu2)) THEN
				local_cs_3(ix,iy,ialpha)=local_cs_3(ix,iy,ialpha)+wt/((2*pi*sina*dalpha)*dx*dy)
				local_Q_cs_3(ixp,iyp,ialpha)=local_Q_cs_3(ixp,iyp,ialpha)+wt/((2*pi*sina*dalpha)*dxp*dyp)*stokesQ
				local_U_cs_3(ixp,iyp,ialpha)=local_U_cs_3(ixp,iyp,ialpha)+wt/((2*pi*sina*dalpha)*dxp*dyp)*stokesU
			ELSE
				local_cs_2(ix,iy,ialpha)=local_cs_2(ix,iy,ialpha)+wt/((2*pi*sina*dalpha)*dx*dy)
				local_Q_cs_2(ixp,iyp,ialpha)=local_Q_cs_2(ixp,iyp,ialpha)+wt/((2*pi*sina*dalpha)*dxp*dyp)*stokesQ
				local_U_cs_2(ixp,iyp,ialpha)=local_U_cs_2(ixp,iyp,ialpha)+wt/((2*pi*sina*dalpha)*dxp*dyp)*stokesU

			END IF

		ELSE IF ((r0.GT.rh).AND.(ABS(th0-0.5*pi).GT.(0.5))) THEN

			IF (E.LT.(nu1)) THEN
				local_pc_1(ix,iy,ialpha)=local_pc_1(ix,iy,ialpha)+wt/((2*pi*sina*dalpha)*dx*dy)
				local_Q_pc_1(ixp,iyp,ialpha)=local_Q_pc_1(ixp,iyp,ialpha)+wt/((2*pi*sina*dalpha)*dxp*dyp)*stokesQ
				local_U_pc_1(ixp,iyp,ialpha)=local_U_pc_1(ixp,iyp,ialpha)+wt/((2*pi*sina*dalpha)*dxp*dyp)*stokesU

			ELSE IF (E.GT.(nu2)) THEN
				local_pc_3(ix,iy,ialpha)=local_pc_3(ix,iy,ialpha)+wt/((2*pi*sina*dalpha)*dx*dy)
				local_Q_pc_3(ixp,iyp,ialpha)=local_Q_pc_3(ixp,iyp,ialpha)+wt/((2*pi*sina*dalpha)*dxp*dyp)*stokesQ
				local_U_pc_3(ixp,iyp,ialpha)=local_U_pc_3(ixp,iyp,ialpha)+wt/((2*pi*sina*dalpha)*dxp*dyp)*stokesU

			ELSE
				local_pc_2(ix,iy,ialpha)=local_pc_2(ix,iy,ialpha)+wt/((2*pi*sina*dalpha)*dx*dy)
				local_Q_pc_2(ixp,iyp,ialpha)=local_Q_pc_2(ixp,iyp,ialpha)+wt/((2*pi*sina*dalpha)*dxp*dyp)*stokesQ
				local_U_pc_2(ixp,iyp,ialpha)=local_U_pc_2(ixp,iyp,ialpha)+wt/((2*pi*sina*dalpha)*dxp*dyp)*stokesU

			END IF
		END IF

	END IF

ENDDO
!!$OMP END DO
!!$OMP CRITICAL


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


!!$OMP END CRITICAL
!!$OMP END PARALLEL


! write(*,*) 'Final rank',i,NGEO,OK1,OK2,OK3,OK4,OK5

CLOSE(unit_in)
CLOSE(unit_out)

END SUBROUTINE IMAGE_GEODESICS_PARALLEL

END PROGRAM LOOP_IMAGE


subroutine INIT_RANDOM_SEED(it,id)

	integer, intent(in), optional :: it,id
	integer :: i, n, clock
	integer, dimension(:), allocatable :: seed
	
	call RANDOM_SEED(size = n)
	
	allocate(seed(n))
	
	call SYSTEM_CLOCK(COUNT=clock)
	
	if ((.not. present(it)) .and. (.not. present(id))) then
	  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
	elseif (present(it) .and. present(id)) then
	  seed = clock + 37 * (/ (i - 1 + it + id*2*FDUMP*NDUMP , i = 1, n) /)
	else
	  print*, 'When initialize random seed, use either both or none arguments'
	  stop
	endif
	
	call RANDOM_SEED(PUT = seed)
	
	deallocate(seed)
	
end subroutine INIT_RANDOM_SEED
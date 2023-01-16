
MODULE MOD_INPUT

IMPLICIT NONE

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++ CONSTANTS +++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
INTEGER, PARAMETER, PUBLIC :: VECWIDTH = 4

! Speed of light 
DOUBLE PRECISION, PARAMETER, PUBLIC           :: c  = 1.0d0
! Fundamental charge [esu]
DOUBLE PRECISION, PARAMETER, PUBLIC           :: e  = 1.0d0
! Mass of the electron [g]
DOUBLE PRECISION, PARAMETER, PUBLIC           :: me = 1.0d0
! pi
DOUBLE PRECISION, PARAMETER, PUBLIC           :: pi = acos(-1.0d0)               
DOUBLE PRECISION, PARAMETER, PUBLIC           :: h  = 1.0d0
DOUBLE PRECISION, PARAMETER, PUBLIC           :: evtoerg = 1.0d0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++ SAVE/RESTORE ++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Save particle and field data for future restoration of the simulation
LOGICAL, PARAMETER, PUBLIC :: CHECKPOINT = .FALSE.

! Restore a simulation where it stopped
LOGICAL, PARAMETER, PUBLIC :: RESTORE = .TRUE.

! Give the wall time in seconds at which the simulation should save to disk 
DOUBLE PRECISION, PARAMETER, PUBLIC :: tstop = 86000

! Give the time step from which the simulation should restart
INTEGER, PARAMETER, PUBLIC :: time_ref = 41469

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++ BOUNDARY CONDITIONS ++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_PART_THMIN="REFLECT"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_PART_THMAX="REFLECT"

! Switch whether EM fields are absorbed or anchored in the PML
! ABSORBED is more convenient for a decreasing field, ANCHORED for Wald (uniform)
CHARACTER(LEN=10), PARAMETER, PUBLIC :: FIELDS_PML="ABSORBED"

! Angular aperture of the absorbing zone in the pml
DOUBLE PRECISION, PARAMETER, PUBLIC :: lambda=acos(-1.0d0)/(12d0)

! Is there a conductor in the simulation?
LOGICAL, PARAMETER, PUBLIC :: CONDUCTOR = .TRUE.

! Initial magnetic configuration: MONOPOLE or WALD or SPLIT or PARABOLOID
CHARACTER(LEN=10), PARAMETER, PUBLIC :: MAG_CONF="PARABOLOID"
DOUBLE PRECISION, PARAMETER, PUBLIC  :: r0=10.0
DOUBLE PRECISION, PARAMETER, PUBLIC  :: nu=3.0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++ INPUT PARAMETERS ++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
! Number of cells in R 
INTEGER*8, PARAMETER, PUBLIC :: NCR=1024

! Number of cells in THETA
INTEGER*8, PARAMETER, PUBLIC :: NCTH=256

! Number of cells in PHI
INTEGER*8, PARAMETER, PUBLIC :: NCPH=512

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Number of particles per cell per species
INTEGER*8, PARAMETER, PUBLIC :: PPC=1
 
! Number of process (domain decomposition in the R- and THETA- and PHI-directions)
INTEGER, PARAMETER, PUBLIC :: NPR=64
INTEGER, PARAMETER, PUBLIC :: NPTH=2
INTEGER, PARAMETER, PUBLIC :: NPPH=32

! Mass ratio IONS/ELECTRONS
DOUBLE PRECISION, PARAMETER, PUBLIC :: mass_ratio=1.0

! Dimensionless black hole spin parameter
DOUBLE PRECISION, PARAMETER, PUBLIC :: a=0.99

! Radius of the horizon in units of rg
DOUBLE PRECISION, PARAMETER, PUBLIC :: rh=1.0+dsqrt(1.0-a*a)

! Radius of the neutron star in units of rg of the black hole
DOUBLE PRECISION, PARAMETER, PUBLIC :: rns=1.0

! Light cylinder radius of the neutron star
DOUBLE PRECISION, PARAMETER, PUBLIC :: rlc=3.0

! Initial orbital separation in units of rg
DOUBLE PRECISION, PARAMETER, PUBLIC :: d0=4.0

! Pulsar magnetic obliquity
DOUBLE PRECISION, PARAMETER, PUBLIC :: chi=0.0*pi

! Spatial boundaries in the R-direction
DOUBLE PRECISION, PARAMETER, PUBLIC :: rmin=0.90*rh !0.98*rh
DOUBLE PRECISION, PUBLIC            :: rmax=10.0

! Radius of the absorbing layer in units of rmax
DOUBLE PRECISION, PUBLIC            :: rpml=0.9

! Maximum value of sigma * dt in the absorbing layer
! Set to 0 when using the characteristic variable BC (CVBC)
DOUBLE PRECISION, PARAMETER, PUBLIC :: sigma_pml = 40.0

! Spatial boundaries in the THETA-direction
DOUBLE PRECISION, PARAMETER, PUBLIC :: thmin=0.03*pi
DOUBLE PRECISION, PARAMETER, PUBLIC :: thmax=pi-thmin

! Spatial boundaries in the PHI-direction
DOUBLE PRECISION, PARAMETER, PUBLIC :: phmin=0.0
DOUBLE PRECISION, PARAMETER, PUBLIC :: phmax=2.*pi

! Dump data frequency in terms of timesteps
INTEGER, PARAMETER, PUBLIC :: FDUMP=1!500

! Number of data dumps
INTEGER, PARAMETER, PUBLIC :: NDUMP=100

! Number of time steps
INTEGER, PARAMETER, PUBLIC :: NT=FDUMP*NDUMP

! Poisson solver calling frequency in terms of timesteps
INTEGER, PARAMETER, PUBLIC :: FREQ_POISSON=25

! Number of iterations to solve Poisson's equation
INTEGER, PARAMETER, PUBLIC :: NIT=500

! Number of filter passes
INTEGER, PARAMETER, PUBLIC :: NFILTER=0

! Time step (as a fraction of the Courant-Friedrichs-Lewy timestep)
DOUBLE PRECISION, PUBLIC :: dt=0.5

! Radiation reaction limit for inverse compton losses
DOUBLE PRECISION, PARAMETER, PUBLIC :: gammaradIC=1d9

!+++++++++++++++++++++++++++++++++++++++++++++ Synchrotron photons parameters
!
! ! Radiation reaction limit for synchrotron losses
 DOUBLE PRECISION, PARAMETER, PUBLIC :: gammaradSC=70000.0
!
! ! Radiation reaction limit for inverse compton losses
! DOUBLE PRECISION, PARAMETER, PUBLIC :: gammaradIC=1d9

! Switch ON/OFF production of synchrotron photons
LOGICAL, PARAMETER, PUBLIC :: SYNC_PHOTONS=.TRUE.

! Rescaling parameter kappa= \tilde{B}_0^2 r_0 / r_g
!DOUBLE PRECISION, PARAMETER, PUBLIC :: kappa = 1e-2

! Proportion of dumped synchrotron photons
DOUBLE PRECISION, PARAMETER, PUBLIC :: frac_sync= 0.02

! Frequency boundaries for the synchrotron spectrum in unit of eB_0 / mc
DOUBLE PRECISION, PARAMETER, PUBLIC :: numin=1d-3,numax=1d5

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Temperature in unit of me*c^2/k of the pairs in the lab frame
DOUBLE PRECISION, PARAMETER, PUBLIC :: temp_inj = 1.

! Fiducial Larmor radius of the electrons
DOUBLE PRECISION, PARAMETER, PUBLIC :: rho0 = 0.0001

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Maximum radius of particle injection
DOUBLE PRECISION, PARAMETER, PUBLIC :: rinj2=4.

! Minimum radius of particle injection
DOUBLE PRECISION, PARAMETER, PUBLIC :: rinj1=rh

! Maximum th of particle injection
DOUBLE PRECISION, PARAMETER, PUBLIC :: thinj2=thmax

! Minimum th of particle injection
DOUBLE PRECISION, PARAMETER, PUBLIC :: thinj1=thmin

! Maximum ph of particle injection
DOUBLE PRECISION, PARAMETER, PUBLIC :: phinj2=phmax

! Minimum ph of particle injection
DOUBLE PRECISION, PARAMETER, PUBLIC :: phinj1=phmin

! Minimum plasma magnetization to have particle injection
DOUBLE PRECISION, PARAMETER, PUBLIC :: sigma_min_inject = 50.

! Minimum DdotB/B0^2 to have particle injection
DOUBLE PRECISION, PARAMETER, PUBLIC :: DdotB_min_inject =1000.0 ! 0.01 ! Too high: no injection

! Normalization of charge released every timestep
! dn = (rate/4 pi q r_g) * D.B/B
DOUBLE PRECISION, PARAMETER, PUBLIC :: rate = 0.5

! Injection rate for initial photons
DOUBLE PRECISION, PARAMETER, PUBLIC :: rate_ph = 0.0 !0.001

! Fiducial optical depth
DOUBLE PRECISION, PARAMETER, PUBLIC :: tau0=30.0

! Switch ON/OFF density test
LOGICAL, PARAMETER, PUBLIC :: TEST_DENSITY=.TRUE.

! Inner radius of conducting disk
DOUBLE PRECISION, PARAMETER, PUBLIC :: rin=6.0

! Innermost Stable Circular Orbit
DOUBLE PRECISION, PARAMETER, PUBLIC :: risco=99999.!6.0

! Accretion flow velocity
DOUBLE PRECISION, PARAMETER, PUBLIC :: v0=0.05!0.05

! Angular size of conducting disk
DOUBLE PRECISION, PARAMETER, PUBLIC :: thin=0.04

! Downsampling factor of photon production
DOUBLE PRECISION, PARAMETER, PUBLIC :: nic=1.0

! Efficiency of the IC cooling of a particle
!DOUBLE PRECISION, PARAMETER, PUBLIC :: alpha=1.0

! Energy boundaries for the particles
DOUBLE PRECISION, PARAMETER, PUBLIC :: gmin=1d0,gmax=1d5

! Lorentz factor boundaries for the particles in the computation of opacities
DOUBLE PRECISION, PARAMETER, PUBLIC :: gdmin=1.01,gdmax=1d4 ! 10.0 ?????

! Grid point for the Lorentz factor array
INTEGER, PARAMETER, PUBLIC :: NGAMD=500

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Characteristic scales which are set by spin a and rho0

! Angular velocity of the horizon
DOUBLE PRECISION, PARAMETER, PUBLIC :: OmegaH = a / (rh**2 + a**2)

! Keplerian angular velocity at the inner radius
DOUBLE PRECISION, PARAMETER, PUBLIC :: OmegaKin=1./((rin)**(3./2.)+a)

! Fiducial field strength
DOUBLE PRECISION, PARAMETER, PUBLIC :: B0 = me*c**2/(e*rho0)

! Fiducial plasma number density - set by approximate GJ density
DOUBLE PRECISION, PARAMETER, PUBLIC :: n0 = (OmegaH+OmegaKin)*B0/(4.*pi*c*e)

! ... or directly by the magnetic energy density i.e. the density giving sigma = 1
!DOUBLE PRECISION, PARAMETER, PUBLIC :: n0 = B0**2 / (4 * pi * me * c**2)

! The fiducial-state magnetization implied by B0 and n0 - note assuming charge-separated
DOUBLE PRECISION, PARAMETER, PUBLIC :: sigma0 = B0**2 / (4 * pi * n0 * me * c**2)

! The skin depth implied by n0; identical to de = c/omega_pe = c*SQRT(me/(4*pi*n0*e**2)) 
DOUBLE PRECISION, PARAMETER, PUBLIC :: de = dsqrt(sigma0) * rho0 

DOUBLE PRECISION, PARAMETER, PUBLIC :: density_threshold=100.0*n0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++ INITIAL PARTICLES ++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Shape of the background distribution: monoenergetic ("MONO") or Maxwellian ("MAXWELL")
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BACKGROUND_PARTICLES="MONO"

! If MONO
DOUBLE PRECISION, PARAMETER, PUBLIC :: gam_init=1.1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++ BACKGROUND RADIATION ++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Switch ON/OFF emission of photons
LOGICAL, PARAMETER, PUBLIC :: EMIT_PHOTONS=.TRUE.

! Switch ON/OFF pair production
LOGICAL, PARAMETER, PUBLIC :: CREATE_PAIRS=.TRUE.

! Shape of the background distribution: monoenergetic ("MONO") or power law ("POWER")
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BACKGROUND_PHOTONS="MONO"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! If MONO

! It must be much smaller than 1, otherwise the simpyfing approximation in the pair
! creation process breaks down and the code fails (max \approx 0.5)
DOUBLE PRECISION, PARAMETER, PUBLIC :: esoft=0.01

! Initial energy of the photons in the simulation, in units of m_e c^2
DOUBLE PRECISION, PARAMETER, PUBLIC :: eph = 1.1

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! If POWER

! Spectral index: dNph/de~e^-pp
DOUBLE PRECISION, PARAMETER, PUBLIC :: pp=2.0

! Energy boundaries for the background radiation in unit of mc^2
DOUBLE PRECISION, PARAMETER, PUBLIC :: e0min=1d-8,e0max=5d-1

! Grid point for the background photon energy array
INTEGER, PARAMETER, PUBLIC :: NEPS0=1000

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++ GAMMA-RAY RADIATION ++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Energy boundaries for the gamma-ray radiation in unit of mc^2
DOUBLE PRECISION, PARAMETER, PUBLIC :: e1min=1d-3,e1max=gmax

! Grid point for the gamma-ray photon energy array
INTEGER, PARAMETER, PUBLIC :: NEPS1=1000

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Shouldn't usually have to modify below this line,
!  unless switching analysis routines on or off.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Number of grid cells in R
INTEGER, PARAMETER, PUBLIC :: NR=NCR+2

! Number of grid cells in THETA
INTEGER, PARAMETER, PUBLIC :: NTH=NCTH+2

! Number of grid cells in PHI
INTEGER, PARAMETER, PUBLIC :: NPHI=NCPH+2

! Number of cell per process in the R-direction
INTEGER, PARAMETER, PUBLIC :: NCRP=NCR/NPR

! Number of cell per process in the THETA-direction
INTEGER, PARAMETER, PUBLIC :: NCTHP=NCTH/NPTH

! Number of cell per process in the PHI-direction
INTEGER, PARAMETER, PUBLIC :: NCPHP=NCPH/NPPH

! Number of grid cells per process in R
INTEGER, PARAMETER, PUBLIC :: NRP=NCRP+2
 
! Number of grid cells per process in THETA
INTEGER, PARAMETER, PUBLIC :: NTHP=NCTHP+2

! Number of grid cells per process in PHI
INTEGER, PARAMETER, PUBLIC :: NPHP=NCPHP+2

! Total number of particles per species NP=NCR*NCTH*NCPH*PPC
INTEGER*8, PARAMETER, PUBLIC :: NP=NCR*NCTH*NCPH*PPC

! Spatial step
DOUBLE PRECISION, PUBLIC            :: dr
DOUBLE PRECISION, PUBLIC            :: dth
DOUBLE PRECISION, PUBLIC            :: dph

#ifdef RLOG
DOUBLE PRECISION, PUBLIC            :: dlogr
#endif
#ifdef THCOS
DOUBLE PRECISION, PUBLIC            :: dcosth
#endif

! Initial field components
DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP,1:NPHP), PUBLIC :: Bruin, Bthuin, Bphuin
DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP,1:NPHP), PUBLIC :: Druin, Dthuin, Dphuin
DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP,1:NPHP), PUBLIC :: Hrdin, Hthdin, Hphdin
DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP,1:NPHP), PUBLIC :: Erdin, Ethdin, Ephdin

! Frozen field components
DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP,1:NPHP), PUBLIC :: Brufz, Bthufz, Bphufz
DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP,1:NPHP), PUBLIC :: Drufz, Dthufz, Dphufz
DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP,1:NPHP), PUBLIC :: Hrdfz, Hthdfz, Hphdfz
DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP,1:NPHP), PUBLIC :: Erdfz, Ethdfz, Ephdfz

! Nodal grid in each domain
DOUBLE PRECISION, DIMENSION(1:NRP), PUBLIC  :: rgp
DOUBLE PRECISION, DIMENSION(1:NTHP), PUBLIC :: thgp,dthetap
DOUBLE PRECISION, DIMENSION(1:NPHP), PUBLIC :: phgp

! Yee grid in each domain
DOUBLE PRECISION, DIMENSION(1:NRP), PUBLIC  :: ryeep
DOUBLE PRECISION, DIMENSION(1:NTHP), PUBLIC :: thyeep
DOUBLE PRECISION, DIMENSION(1:NPHP), PUBLIC :: phyeep

!***********************************************************************
! 1- Metric at (i,j)
! 2- Metric at (i+1/2,j)
! 3- Metric at (i+1/2,j+1/2)
! 4- Metric at (i,j+1/2)

! COVARIANT components
DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP), PUBLIC :: grd1,grphd1,gthd1,gphd1
DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP), PUBLIC :: gru1,grphu1,gthu1,gphu1

DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP), PUBLIC :: grd2,grphd2,gthd2,gphd2
DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP), PUBLIC :: gru2,grphu2,gthu2,gphu2

DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP), PUBLIC :: grd3,grphd3,gthd3,gphd3
DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP), PUBLIC :: gru3,grphu3,gthu3,gphu3

DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP), PUBLIC :: grd4,grphd4,gthd4,gphd4
DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP), PUBLIC :: gru4,grphu4,gthu4,gphu4

! Square-root of the determinant
DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP), PUBLIC :: sgam1,sgamtild1
DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP), PUBLIC :: sgam2,sgamtild2
DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP), PUBLIC :: sgam3,sgamtild3
DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP), PUBLIC :: sgam4,sgamtild4

! Lapse function alpha and CONTRAVARIANT shift vector (beta^r component only)
DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP), PUBLIC :: alpha1,alpha2,alpha3,alpha4
DOUBLE PRECISION, DIMENSION(1:NRP,1:NTHP), PUBLIC :: beta_ru1,beta_ru2,beta_ru3,beta_ru4

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++ DEFINE PARTICLE ARRAYS +++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Array of Structures of Arrays data structure for particles.
! p_species(:,1,:)=r      : r position of the particle
! p_species(:,2,:)=th     : theta-angle of the particle
! p_species(:,3,:)=ph     : phi-angle of the particle
! p_species(:,4,:)=ur     : r-component of the reduced particle 4-velocity
! p_species(:,5,:)=uth    : theta-component of the reduced particle 4-velocity
! p_species(:,6,:)=uph    : phi-component of the reduced particle 4-velocity
! p_species(:,7,:)=weight : weight of the particle
!
! Defition of the particle distribution data function: it gives extra
! information regarding the particles.
!
! pd_species(:,1,:)=vr : Coordinate velocities used for the deposit
! pd_species(:,2,:)=vth
! pd_species(:,3,:)=vph 
! 
!***********************************************************************

! ELECTRONS distribution function components
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: p_e(:,:,:)
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pd_e(:,:,:)
INTEGER*8, ALLOCATABLE, PUBLIC        :: tage(:,:)

! IONS distribution function components
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: p_i(:,:,:)
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pd_i(:,:,:)
INTEGER*8, ALLOCATABLE, PUBLIC        :: tagi(:,:)

! PHOTONS distribution function components
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: p_ph(:,:,:)
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pd_ph(:,:,:)
INTEGER*8, ALLOCATABLE, PUBLIC        :: tagph(:,:)

! PHOTONS distribution function components
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: p_bt(:,:,:)
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pd_bt(:,:,:)
INTEGER*8, ALLOCATABLE, PUBLIC        :: tagbt(:,:)

! Temporary particle distribution function
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: p_f(:,:,:)
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pd_f(:,:,:)
INTEGER*8, ALLOCATABLE, PUBLIC        :: tagf(:,:)

! Subsampling factor for low-energy photons push
!INTEGER, PARAMETER, PUBLIC :: sub_bt=1

! Subsampling factor for low-energy photons push
DOUBLE PRECISION, PARAMETER, PUBLIC :: frac_bt=0.2

! Starting time for geokerr dumping
DOUBLE PRECISION, PARAMETER, PUBLIC :: start_bt=0.0

! Size of the box for below threshold photons, must be greater than rmax
DOUBLE PRECISION, PUBLIC            :: rmax_bt=10.0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++ ANALYSIS ++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Switch PARALLEL dumping mode ON/OFF for the density, pressure, current,
! and field arrays.
! *Recommended for large simulations*
LOGICAL, PARAMETER, PUBLIC :: DUMP_FIELDS_PARALLEL=.FALSE.

! Switch ON/OFF write fields to disk
LOGICAL, PARAMETER, PUBLIC :: WRITE_FIELDS=.TRUE.

! Switch ON/OFF write divergence of D and B
LOGICAL, PARAMETER, PUBLIC :: CHECK_DIV=.FALSE.

! Switch ON/OFF write the charge current densities to disk
LOGICAL, PARAMETER, PUBLIC :: WRITE_RHOJ=.FALSE.

! Switch PARALLEL dumping mode ON/OFF for particle arrays.
! *Recommended for large simulations*
LOGICAL, PARAMETER, PUBLIC :: DUMP_PARTICLES_PARALLEL=.FALSE.

! Switch ON/OFF write a sub-sample of the raw particle data to disk
LOGICAL, PARAMETER, PUBLIC :: WRITE_PARTICLES=.FALSE.

! Switch ON/OFF particle tracker
LOGICAL, PARAMETER, PUBLIC :: TRACK_PARTICLES=.FALSE.

! Switch ON/OFF analysis of particle stress-energy tensor
LOGICAL, PARAMETER, PUBLIC :: ANALYZE_TMUNU=.FALSE.

! Switch ON/OFF analysis the plasma density
LOGICAL, PARAMETER, PUBLIC :: ANALYZE_DENSITIES=.FALSE.

! Switch ON/OFF analysis the plasma density on a regular Cartesian grid (useful for visu)
LOGICAL, PARAMETER, PUBLIC :: ANALYZE_DENSITIES_CARTESIAN=.TRUE.

! Switch ON/OFF analysis of particle energy spectra
LOGICAL, PARAMETER, PUBLIC :: ANALYZE_SPECTRA=.FALSE.!.TRUE.

! Switch ON/OFF analysis of phase space
LOGICAL, PARAMETER, PUBLIC :: ANALYZE_PHASE_SPACE=.FALSE.

! Switch ON/OFF analysis of the light curve
LOGICAL, PARAMETER, PUBLIC :: WRITE_LIGHTCURVE=.FALSE.

! Max extension of the Cartesian mesh
DOUBLE PRECISION, PUBLIC            :: xmaxi=6.0

! Fraction of sampled particles
DOUBLE PRECISION, PARAMETER, PUBLIC :: frac=1.0

!Angle at which phase space is sampled
DOUBLE PRECISION, PARAMETER, PUBLIC :: theta0=pi/6

!Angular aperture of the sampled phase space
DOUBLE PRECISION, PARAMETER, PUBLIC :: deltath=0.1

! Total number of tracked particles
INTEGER, PARAMETER, PUBLIC :: NSAMPLE=1000

! Give the time step from which the tracker should restart
INTEGER, PARAMETER, PUBLIC :: time_track=1500

! Particle's Lorentz factor grid for the spectrum
INTEGER, PARAMETER, PUBLIC :: NGAM=300

! Particle's photon energy grid for the spectrum
INTEGER, PARAMETER, PUBLIC :: NEPS=100

! Sub-sampling factor of the dumped particles
INTEGER, PARAMETER, PUBLIC :: SUB_PART=10

! Sub-sampling factor of the dumped grid quantities
INTEGER, PARAMETER, PUBLIC :: SUB_NR=4,SUB_NTH=1,SUB_NPH=2

! Reduced number of cell in the R-direction
INTEGER, PARAMETER, PUBLIC :: NRR=NCR/SUB_NR+2

! Reduced number of cell in the THETA-direction
INTEGER, PARAMETER, PUBLIC :: NTHR=NCTH/SUB_NTH+2

! Reduced number of cell in the PHI-direction
INTEGER, PARAMETER, PUBLIC :: NPHR=NCPH/SUB_NPH+2

! Reduced number of cell per process in the R-direction
INTEGER, PARAMETER, PUBLIC :: NRRP=NCRP/SUB_NR+2

! Reduced number of cell per process in the THETA-direction
INTEGER, PARAMETER, PUBLIC :: NTHRP=NCTHP/SUB_NTH+2

! Reduced number of cell per process in the PHI-direction
INTEGER, PARAMETER, PUBLIC :: NPHRP=NCPHP/SUB_NPH+2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Reduced global grid for data analysis
DOUBLE PRECISION, DIMENSION(1:NRR), PUBLIC  :: rr
DOUBLE PRECISION, DIMENSION(1:NTHR), PUBLIC :: thr
DOUBLE PRECISION, DIMENSION(1:NPHR), PUBLIC :: phr

! Reduced local grid for data analysis
DOUBLE PRECISION, DIMENSION(1:NRRP), PUBLIC  :: rrp
DOUBLE PRECISION, DIMENSION(1:NTHRP), PUBLIC :: thrp
DOUBLE PRECISION, DIMENSION(1:NPHRP), PUBLIC :: phrp

! Reduced global Cartesian grid for data analysis
DOUBLE PRECISION, DIMENSION(1:NRR), PUBLIC :: xr
DOUBLE PRECISION, DIMENSION(1:NRR), PUBLIC :: yr
DOUBLE PRECISION, DIMENSION(1:NRR), PUBLIC :: zr

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!********************** Lightcurve diagnostic
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Number of latitude angles in the lightcurve diagnostic
INTEGER, PARAMETER, PUBLIC :: NALPHA=20

! Number of light crossing times before the lightcurve diagnostic starts
DOUBLE PRECISION, PARAMETER, PUBLIC :: start_lc=0.0

! Number of time steps between every calculation of an instantaneous light curve
INTEGER, PARAMETER, PUBLIC :: N_LIGHT=1 ! N_LIGHT must divide FDUMP !

! Under-sampling of the time grid on which the lightcurve is computed
DOUBLE PRECISION, PARAMETER, PUBLIC :: N_UNDER=20.0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE MOD_INPUT

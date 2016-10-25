!=============================================
! Module to use SPT pol B mode data.
! Based on SPT publicly available likelihood
! at www. ......
! and cliklike.f90
! 
! Version: 0.1
!
! Written by Alex Zucca, azucca@sfu.ca
!=============================================

module SPT_pol_BB

 use CosmologyTypes
 use CosmoTheory
 use settings
 use Likelihood_Cosmology


use MatrixUtils
use Calculator_Cosmology
use IniObjects


 implicit none

 logical :: use_SPTpol = .false.

 integer, parameter :: dp = kind(1.d0)

 ! Bandpowers are read in this order (all BB):
 !            (90x90), (90x150), (150x150)
 integer :: nfreq, nband, nbin, nall, bands_per_freq
 double precision, dimension(:,:), allocatable :: spec, cov
 double precision, dimension(:,:), allocatable :: windows
 integer :: spt_windows_lmin, spt_windows_lmax
 double precision, dimension(:), allocatable :: cl_to_dl_conversion, dl_th, poisson_constant
 double precision, dimension(:), allocatable :: dust_template, tensor_spectrum, scalar_spectrum

 logical :: SuccessfulSPTInitialization
 logical :: use_dust_gaussian_prior


!Now introducing the SPT likelihood class
 type, extends(TCMBLikelihood) :: SPTpolLikelihood
    integer(kind=4) :: spt_lmax
 contains
   procedure :: LogLike => SPTpolLnLike
   procedure :: SPTpol_likeinit
 end type SPTpolLikelihood

 private

 public :: use_SPTpol, SPT_ReadParams

 contains

 subroutine SPT_ReadParams(LikeList,Ini)
   use IniFile
   use settings
   class(TLikelihoodList) :: LikeList
   class(TSettingIni) Ini
   class(SPTpolLikelihood), pointer :: Like
   character (LEN=Ini_max_string_len) :: desc_file
   character (LEN=Ini_max_string_len) :: bp_file
   character (LEN=Ini_max_string_len) :: cov_file
   character (LEN=Ini_max_string_len) :: dust_file
   character (LEN=Ini_max_string_len) :: tensor_file
   character (LEN=Ini_max_string_len) :: scalar_file
   character (LEN=Ini_max_string_len) :: win_folder

   write(*,*) "Reading SPT params.."

   if (.not. Ini%Read_Logical('use_SPTpol', .false.)) then
      return
   end if


   use_dust_gaussian_prior = Ini%Read_Logical('use_dust_gaussian_prior', .false.)

   allocate(SPTpolLikelihood::like)
   Like%speed = 10
   Like%LikelihoodType = 'CMB'
   Like%name = 'SPTpol_bb100d'
   Like%needs_background_functions = .true.

   print *, 'Likelihood Name: ', Like%name

   call Like%loadParamNames(trim(DataDir)//'SPTpol_mag.paramnames')
   call LikeList%Add(like)

   desc_file = Ini%Read_String_Default('sptpol_desc_file', &
               trim(DataDir)//'bb100d.desc_file')

   !Read in the bandpower and cov files, and determine where the windows are.
   bp_file     = Ini%Read_String_Default('sptpol_bp_file','')
   cov_file    = Ini%Read_String_Default('sptpol_cov_file','')
   dust_file   = Ini%Read_String_Default('sptpol_dust_file','')
   tensor_file = Ini%Read_String_Default('sptpol_tensor_file','')
   scalar_file = Ini%Read_String_Default('sptpol_scalar_file','')
   win_folder  = Ini%Read_String_Default('sptpol_window_dir','')

   if (bp_file=='' .or. desc_file=='' .or. cov_file=='' .or. dust_file=='' .or. tensor_file=='' .or. scalar_file=='' .or. win_folder=='' ) then
     print*,'Missing required sptpol key: received: ',bp_file,desc_file,cov_file,dust_file,tensor_file,scalar_file,win_folder
     stop
   endif

   call Like%SPTpol_likeinit(desc_file, bp_file, cov_file, dust_file, tensor_file, scalar_file, win_folder)

 end subroutine SPT_ReadParams



 subroutine SPTpol_likeinit(this,desc_file, bp_file, cov_file, dust_file, tensor_file, scalar_file, win_folder)
   use IniFile
   class (SPTpolLikelihood) :: this

   character(LEN=Ini_max_string_len) :: desc_file, bp_file, cov_file, dust_file, tensor_file, scalar_file
   character(LEN=Ini_max_string_len) :: win_folder
   integer i,j,k,dum,lwin
   real rtmp
  
   type(TTextFile) :: F

   !Obtain necessary info from the desc_file pertaining
   !to which freqs we're using, ell range, and number of bins per spectrum.
   call F%Open(desc_file)
   print *, 'Reading desc file...',trim(desc_file)

   read (F%Unit,*) nbin,nfreq !Number of bandpowers per spectrum, Number of freqs
   read (F%Unit,*) spt_windows_lmin, spt_windows_lmax !Min and Max ell in window files.

   call F%Close()
 
   print *, 'nbin: ', nbin
   print *, 'nfreq: ', nfreq
   print *, 'spt_windows_lmin: ', spt_windows_lmin
   print *, 'spt_windows_lmax: ', spt_windows_lmax
   print *, 'win_folder: ', trim(win_folder)

   !We have only BB for three cross-frequencies
   bands_per_freq = 1

   nband = bands_per_freq*(nfreq)
   nall = nbin*(nfreq*bands_per_freq)

   print *, 'nband: ', nband
   print *, 'nall: ', nall
   print *, 'bands_per_freq: ', bands_per_freq

 !if (lmax .lt. spt_windows_lmax) then
 ! B-modes should be in Cls(3,3)
 !if (CosmoSettings%cl_lmax(3,3) .lt. spt_windows_lmax) then
 !spt_windows_lmax=CosmoSettings%cl_lmax(3,3)
 !end if

   if (spt_windows_lmin < 2 .or. spt_windows_lmin >= spt_windows_lmax) then
     call mpistop('Invalid lranges for sptpol')
   end if

   if (spt_windows_lmin .lt. 2) then
     print *, "Unexpected lmin for windows!"
     stop
   end if

   allocate(windows(spt_windows_lmin:spt_windows_lmax,1:nall), &
   spec(1:nbin,nband))

   allocate(cov(1:nall,1:nall))

   allocate(poisson_constant(spt_windows_lmin:spt_windows_lmax))
   do j=spt_windows_lmin,spt_windows_lmax
      poisson_constant(j) = 1.0d0
   end do

   allocate(dust_template(spt_windows_lmin:spt_windows_lmax), &
   tensor_spectrum(spt_windows_lmin:spt_windows_lmax), &
   scalar_spectrum(spt_windows_lmin:spt_windows_lmax), &
   dl_th(spt_windows_lmin:spt_windows_lmax))

   !Read in bandpowers
   !Should be 90x90, 90x150, 150x150 BB in that order from SPTpol analysis.
   print *, 'Reading in bandpowers...'
   print *, 'Reading bp_file...',trim(bp_file)

   call F%Open(bp_file)
   do i=1,nband
      do j=1,nbin
         read (F%Unit,*) dum,spec(j,i)
      end do
   end do
  
   call F%Close()
   print *, 'Bandpowers are:', spec

   !Read in dust template (150 GHz).
   print *, 'Reading in dust template...'
   print *, 'Reading dust file...',trim(dust_file)

   call F%Open(dust_file)
   do j=spt_windows_lmin,spt_windows_lmax
      read (F%unit,*) dum,dust_template(j)
   end do
   call F%Close()

   !Read in pre-calculated tensor BB.
   ! This is useless for PMFs. I keep this part, but it will not be used.
   print *, 'Reading in tensor spectrum template...'
   call F%Open(tensor_file)
   do j=spt_windows_lmin,spt_windows_lmax
      read (F%unit,*) dum, tensor_spectrum(j)
   end do
   call F%Close()

   !Read in pre-calculated scalar BB.
    print *, 'Reading in scalar spectrum template...'
    call F%Open(scalar_file)
   do j=spt_windows_lmin,spt_windows_lmax
      read (F%unit,*) dum, scalar_spectrum(j)
   end do
   call F%Close()


   !Read in covariance
   !Should be 90x90, 90x150, 150x150 BB
   print *, 'Reading in cov_data matrix...'
   call F%Open(cov_file)
   do i=1,nall
      do j=1,nall
         read (F%unit,*) cov(i,j)
      end do
   end do

   call F%Close()

   call Matrix_CholeskyDouble(cov)

   ! Read in windows
   !Should be 90x90, 90x150, 150x150
   print *, 'Reading in windows...'
   do i=1,nall
      call F%Open(trim(win_folder)//trim(numcat('/window_',i)))
      do j=spt_windows_lmin,spt_windows_lmax
         read (F%unit,*) dum, windows(j,i)
      end do
      call F%Close()
   end do

   !Setting l_max
   allocate(this%cl_lmax(CL_B,CL_B), source = 0)
   this%cl_lmax = 0
   this%cl_lmax(CL_B,CL_B) = spt_windows_lmax
   where (this%cl_lmax<0)
         this%cl_lmax=0
   end where


   SuccessfulSPTInitialization = .true.

 end subroutine SPTpol_likeinit

 !Likelihood Function
 function SPTpolLnLike(this, CMB, Theory, DataParams)
    class(SPTpolLikelihood) :: this
    class(CMBParams) :: CMB

    class(TCosmoTheoryPredictions), target :: Theory !cl(lmax, num_tot_cls) TT, TE, EE, BB,(+phi) in that order from Camb.

    real(mcp) :: SPTpolLnLike
    real(mcp) :: A_dustLnLike
    real(mcp) :: dum
    real(mcp) :: DataParams(:) !Foreground parameters (plus wimp annihilation vector amplitude)
    double precision, parameter :: sigmaA_dust = 0.30
    double precision, parameter :: meanA_dust = 1.0 ! Scales all cross-frequencies.
    double precision, dimension(1:nall) :: deltacb
    double precision, dimension(1:nbin) :: tmpcb
    double precision, dimension(1) :: junk
    integer :: j,k,l, i !ALEX
    real(mcp) :: norm


    do j=1,(nfreq)
       tmpcb(:)=0

      !DO EVERYTHING IN CL.
      !DataParams: [A_lens, A_tens, A_dust, A_PS90, A_PS90150, A_PS150]
      !First get model CMB spectra (in Cl).
      !AZ
      !         dl_th(:) = DataParams(1)*scalar_spectrum(:) + DataParams(2)*tensor_spectrum(:)

      !AZ: dl_th must be substituted with my own prediction...

      do i = spt_windows_lmin, spt_windows_lmax
         dl_th(i) = Theory%Cls(3,3)%CL(i)*twopi/(i*(i+1))  ! it should be ok...
      end do

      !AZ: changing numbers for DataParams..
      !Now add foreground terms.
      !Case 1: 90x90
      if (j==1) then
          dl_th(:) = dl_th(:) + DataParams(1)*0.00169*dust_template(:) + DataParams(2)*4.65*10**(-8.0d0)*poisson_constant(:)
      end if
 
      !Case 2: 90x150
      if (j==2) then
         dl_th(:) = dl_th(:) + DataParams(1)*0.00447*dust_template(:) + DataParams(3)*2.27*10**(-8.0d0)*poisson_constant(:)
      end if

      !Case 2: 150x150
      if (j==3) then
         dl_th(:) = dl_th(:) + DataParams(1)*0.0118*dust_template(:) + DataParams(4)*1.11*10**(-8.0d0)*poisson_constant(:)
      end if
 
      !Now bin into bandpowers with the window functions - output is tmpcb (binned theory+foreground spectrum).
      call dgemv('T',spt_windows_lmax-spt_windows_lmin+1,nbin,1.0d0,&
      windows(:,1+nbin*(j-1):nbin+nbin*(j-1)),spt_windows_lmax-spt_windows_lmin+1,&
      dl_th,1,0d0,tmpcb,1)

 !Finally, get \Delta Cl's between our model Cls and our input spectra
 !for covariance calculation.
 !spec(:,1) = 90x90
 !spec(:,2) = 90x150
 !spec(:,3) = 150x150
      deltacb(1+(j-1)*nbin:nbin+(j-1)*nbin) = tmpcb(:) - spec(:,j)
    end do



    !Now calculate the SPTpol likelihood.

    !Here, the covariance doesn't change with parameter values, so we can just calculate chi^2 for
    !the likelihood.
    SPTpolLnLike = Matrix_GaussianHalfChisqCholDouble(cov,deltacb)

    !Calculate the Gaussian prior term for the dust, ignoring the extra log|cov| term.
    if (use_dust_gaussian_prior== .true.) then
       A_dustLnLike = 0.5*((DataParams(1) - meanA_dust)**2.0d0 / sigmaA_dust**2.0d0)
    else
       A_dustLnLike = 0.0d0
    end if

    print *, 'SPTpolLnLike = ', SPTpolLnLike
    print *, 'A_dustLnLike = ', A_dustLnLike

    SPTpolLnLike = SPTpolLnLike + A_dustLnLike

    if (feedback > 1)      print *, 'SPTpolLnLike = ', SPTpolLnLike
end function SPTpolLnLike



end module SPT_pol_BB

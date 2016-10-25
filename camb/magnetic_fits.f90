module magnetic
  
  use constants
  use AMLutils
  use ModelParams
  use CAMB
  use richsub
  use Transfer

!ALEX:
use InitialPower
!ALEX.

  ! Module for calculating the amplitudes of the magnetic perturbations
  ! to the energy momentum tensor from the magnetic field properties.
  
  implicit none

  !private
  
  character(len= 100) :: mi_filename = 'magnetic_integral.dat'
!  integer, parameter :: mi_nr = 26  
  integer, parameter :: mi_nr = 50
  integer, parameter :: mi_nc = 6

  ! Array of the spectral indices
  real(dl) :: specind(mi_nr)

  ! Array of the calculated angular integrals
  ! In order: delta, delta-pi, pi, pi vector, pi-tensor
  real(dl) :: intval(mi_nr,mi_nc-1)
  
  ! Array of second derivatives for splining
  real(dl) :: intval2(mi_nr,mi_nc-1)

  ! Array of constants for each PS
  real(dl), parameter :: psconst(mi_nc-1) = (/0.25, 0.25, 0.25, 2., 2./3. /)
  !real(dl), parameter :: psconst(mi_nc-1) = (/0.25, 0.5, 1., 1., 1./2. /)!, yun

  ! Constant for calculating the amplitudes = 1/(2*rho_gamma_0)^2
  real(dl), parameter :: amp_constant = 1.432e-12_dl

  !Define constants for specifying type of perturbation
  !integer, parameter :: mag_scal_passive = 1, mag_scal_comp = 2, &
  !     mag_vec_comp = 3, mag_tens_passive = 4, mag_tens_comp = 5
  integer, parameter :: mag_compensated = 1, mag_passive = 2, mag_and_prim = 3

  logical :: notloaded = .true.

contains


  ! Initialise the module. Read files etc.
  subroutine magamp_init

    implicit none

    integer :: i

    open(217,file=mi_filename)

    do i=1,mi_nr
       read(217,*) specind(i), intval(i,:)
    end do

    do i=1,mi_nc-1
       call spline_double(specind,intval(:,i),mi_nr,intval2(:,i))
    end do

    close(217)

    notloaded = .false.

  end subroutine magamp_init

  
  ! Function to calculate the common ampltiude
  ! B_lambda is in units of nano Gauss.
  real(dl) function mag_amplitude(spec_ind, b_lambda)
    
    !use AMLutils

    real(dl), intent(in) :: spec_ind
    real(dl), intent(in) :: b_lambda
    
    
    ! Scale by the spectral index dependent bits, and the field amplitude
    mag_amplitude = amp_constant * b_lambda**4 * (2*pi)**(2*spec_ind + 4) &
         / GAMMA((spec_ind + 3)/2._dl)**2
    
  end function mag_amplitude
    

    



  ! Get integral values by splining loaded files.
  real(dl) function mag_spline_int(mag_ind, corr_num)
    
    real(dl), intent(in) :: mag_ind
    integer, intent(in) :: corr_num

    if(notloaded) call magamp_init

    mag_spline_int = spline_val(mag_ind, specind, &
         intval(:,corr_num), intval2(:,corr_num), mi_nr)

  end function mag_spline_int

    
  ! Function to calculate the amplitude of a particular spectrum
  real(dl) function mag_psamp(spec_ind, b_lambda, corr_num)

    real(dl), intent(in) :: spec_ind
    real(dl), intent(in) :: b_lambda
    integer, intent(in) :: corr_num

    mag_psamp = mag_amplitude(spec_ind, b_lambda) * &
         mag_spline_int(spec_ind, corr_num) * psconst(corr_num)
  end function mag_psamp



  ! Set the parameters needed for adding in the Alfven velocity to get
  ! a magnetic Jeans mass. Must be done separately so we can use for
  ! the non magnetic modes.
  subroutine mag_set_alfven(mag_amp, mag_ind)
    real(dl), intent(in) :: mag_amp, mag_ind

    klambda = 2._dl * pi
    blambda = mag_amp
    nb_ind = mag_ind

  end subroutine mag_set_alfven

  ! Reset back to the default parameters.
  subroutine mag_reset_alfven
    call mag_set_alfven(0._dl, -3._dl)
  end subroutine mag_reset_alfven



 subroutine cls_from_params(CP, Cls, usedata)

   ! Barely needed, only for call to ClsFromTheoryData to get
   ! normalisation

   ! Contains all the parameters for CAMB
   type(CAMBParams), intent(in) :: CP
   
   ! Array to return Cls into
   real, intent(out) :: Cls(lmin:,:)

   logical, intent(in), optional :: usedata

   logical :: ud
   
   type(CAMBdata) :: data
   type(CAMBParams) :: tcp   

   integer :: error, l

   integer, parameter :: ScalClOrder(3) = (/1,3,2/), TensClOrder(4) = (/1,4,2,3/)
   
   ud = .false.
   if(present(usedata)) ud = usedata

   tcp = CP

   if(ud) then 
      ! Initialise data type for holding transfers etc.
      call CAMB_InitCAMBdata(data)
      data%Params = tcp
   
      ! Do hard work in calculating
      call CAMB_GetTransfers(tcp, data, error)
       !Print *, "ERRORS :", error!
      ! Calculate power spectra
      call CAMB_TransfersToPowers(data)

   else
      ! Just call getresults
      call CAMB_GetResults(tcp, error)
    !Print *, "In cls_from_params :", error!, tcp
   end if
   !check
   if(error /= 0) then
      write (*,*) "Eeek error:", error
      return
   end if


   ! Copy out Cls into the return array.
   if(CP%WantScalars) then
      ! do l = lmin, CP%Max_l
      if(CP%DoLensing) then
         Cls(lmin:CP%Max_l,C_Temp:C_last) =  Cl_lensed(lmin:CP%Max_l,1, C_Temp:C_last)
      else 
         Cls(lmin:CP%Max_l,C_Temp:C_last) =  Cl_scalar(lmin:CP%Max_l,1, C_Temp:C_last)
       !abs(Cl_scalar(lmin:CP%Max_l,1, C_Temp:C_last))
        ! print*, Cls(lmin:CP%Max_l,C_Temp:C_last)
         !print*, Cl_scalar(lmin:CP%Max_l,1, C_Temp:C_last) !yun
      end if
      !yun ,check
        !print*, "scalar Cls :"
        !pause
        !print*, Cl_scalar(lmin:CP%Max_l,1, C_Temp:C_last)
      !end do
   else if(CP%WantTensors) then
      do l = 2, CP%Max_l_tensor
         Cls(l,1:4) = Cl_tensor(l,1, TensClOrder(1:4))
      end do
   else if(CP%WantVectors) then
      do l = 2, CP%Max_l
         Cls(l,1:4) =  Cl_vector(l,1, TensClOrder(1:4))
      end do
   end if
   
   ! Free up everything created.
   if(ud) call CAMB_FreeCAMBdata(data)
   
 end subroutine cls_from_params



 ! Calculate the magnetic Cls 
 ! Largely this sets up all the relevant parameters and then calls
 ! cls_from_params, returning the results.
  subroutine mag_cls(CP, magmode, magamp, magind, maglrat, usedata)

    type(CAMBParams) ,intent(inout) :: CP

    integer, intent(in) :: magmode
    double precision, intent(in) :: magamp, magind, maglrat
    logical, intent(in), optional :: usedata
    
    double precision :: amp, camb_ind

    double precision :: Rg, Rv

    double precision :: lrat
    real :: Cltemp(lmin:CP%Max_l,1:4)
    real :: Cltemp2(lmin:CP%Max_l,1:4)
    integer :: i
    logical :: ud

    Rg = 0.0_dl; Rv = 0._dl

    ud = .false.
    if(present(usedata)) ud = usedata
    Cltemp(:,:) = 0.0_dl!yun

    ! Set photon and neutrino fractions
    Rg = 1._dl / (1._dl+7._dl/8. * (4._dl/11.)**(4._dl/3.) * (CP%Num_Nu_massless + CP%Num_Nu_massive))
    Rv = 1._dl - Rg

    if(magmode == 0 .or. .not. CP%WantCls) then
       write (*,*) "Can't calculate Mag Cls if mode = 0, or WantCls is false."
       return
    end if

    !Alex: COMMON STUFF for FITTING FUNCTIONS
    CP%InitPower%kDissip = (5.5d0 * 10.d0**4)**(1.d0/(magind + 5.d0)) * (magamp)**(-2.d0/(magind+5.d0)) * &
    (2*pi)**((magind+3)/(magind+5))*(CP%H0/100)**(1/(magind+5))&
    *((CP%omegab*(CP%H0/100)**2)/0.022)**(1/(magind+5))

    write(*,*) 'kDissip = ', CP%InitPower%kDissip
!this I can neglect from now...
    !Magnetic_Index = magind
    !write(*,*) "Magnetic_index = ", Magnetic_Index


    ! Set common power spectrum stuff
    CP%InitPower%nn = 1
    !CP%InitPower%rat(1) = 1._dl  !yun
    CP%InitPower%n_run(1) = 0._dl
    CP%InitPower%nt_run(1) = 0._dl!yun
    CP%InitPower%k_0_scalar = 2*pi
    CP%InitPower%k_0_tensor = 2*pi

    ! Setup common options
    CP%WantCls = .true.
    CP%WantTransfer = .false.
    !CP%WantScalars = .false.
    !CP%WantVectors = .false.
    !CP%WantTensors = .false.
    CP%OnlyTransfers = .false.

    ! Convert logarithm into natural log
    lrat = maglrat * log(10._dl)

    if(CP%WantTensors) then !yun
! Ensure tensor neutrinos are on (as mag field modes don't make
! any sense otherwise)
    DoTensorNeutrinos = .true.

     if(magmode == mag_passive) then
       CP%WantTensors = .true.!YUN
       delb = 0._dl ! Set up perturbations
       pib = 0._dl
       tens0 = 1._dl !yun
       !Choose between n>0 and n<0
       if(magind<-1.5) then
            if(Feedbacklevel>1) write(*,*) "USING INTERPOLATION TABLE"
            CP%InitPower%ant(1) = 2._dl*(3 + magind)
            CP%InitPower%TensorPowerAmp(1) = mag_psamp(magind, magamp, 5) * (6._dl*Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2!yun
            !amp = mag_psamp(magind, magamp, 5) * (3._dl*Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2 !yun, check
            !print*, CP%InitPower%TensorPowerAmp(1)
            !CP%InitPower%ant(1) = camb_ind
            CP%InitPower%CorrType = 0
        else if (magind .ge. -1.5) then
            !Use first fitting formula
            if(Feedbacklevel>1) write(*,*) "USING FITTING FORMULAS"
            CP%InitPower%ant(1) = magind
            CP%InitPower%TensorPowerAmp(1)=mag_amplitude(magind, magamp)*psconst(5)*(6._dl*Rg*(lrat+(5._dl/(8._dl*Rv)-1)))**2
            CP%InitPower%CorrType = 5
        end if

!ALEX: I must use the fitting functions even for this..
     else if(magmode == mag_compensated) then
       CP%WantTensors = .true.!YUN
       delb = 0._dl ! Set up perturbations
       pib = 1._dl
       tens0 = 0._dl
       camb_ind = 2._dl*(3 + magind)
       CP%InitPower%ant(1) = camb_ind
        !CP%InitPower%TensorPowerAmp(1)= amp!yun
       CP%InitPower%TensorPowerAmp(1) = mag_psamp(magind, magamp, 5)!yun
     end if

    end if  !yun

    if(CP%WantScalars) then
       ! Setup the parameters for calculating the Alfven velocity
       call mag_set_alfven(magamp, magind)

       if(magmode == mag_compensated) then
          ! This section is quite trickey as we need to combine the two
          ! correlated perturbation types. It essentially fetches each
          ! set of Cls with the correct amplitudes (taking into account
          ! the cross correlation), and then sums them up.
          
          ! Set up common parameters
          !camb_ind = 1._dl + 2._dl*(3 + magind)
          CP%Scalar_initial_condition = 6
          CP%WantScalars = .true.
          !CP%InitPower%an(1) = camb_ind
          !yun
!Delta-Delta
          !! The density perturbed mode (delb = 1, pib = 0)
          !CP%InitPower%ScalarPowerAmp(1)= mag_psamp(magind, magamp, 1) - mag_psamp(magind, magamp, 2)
          if(magind .ge. -1.5) then
                CP%InitPower%an(1) =  magind
                CP%InitPower%CorrType = 1
                write(*,*) "Using FITTING Functions"
                CP%InitPower%ScalarPowerAmp(1)=mag_amplitude(magind, magamp)!* psconst(1) - mag_amplitude(magind, magamp)* psconst(2)
          else
                CP%InitPower%an(1) = 1._dl + 2._dl*(3 + magind)
                CP%InitPower%CorrType = 0
                CP%InitPower%ScalarPowerAmp(1)= mag_psamp(magind, magamp, 1) - mag_psamp(magind, magamp, 2)
                write(*,*) "Using interpolation TABLE"
          endif
          delb = 1._dl ! Set up perturbations
          pib = 0._dl
           ! Cltemp = 0.
          call cls_from_params(CP,Cltemp2,ud)
          Cltemp(:,:) = Cltemp(:,:) + Cltemp2(:,:)
          !print*, Cltemp(lmin:CP%Max_l, C_Temp:C_last)!, CP, ud!yun

!Pi-Pi
          !! The stress perturbed mode (delb = 0, pib = 1)
          !CP%InitPower%ScalarPowerAmp(1)= mag_psamp(magind, magamp, 3) - mag_psamp(magind, magamp, 2)
          if(magind .ge. -1.5) then
                CP%InitPower%CorrType = 2
                write(*,*) "Using FITTING Functions"!mag_amplitude(magind, magamp)* psconst(
                CP%InitPower%ScalarPowerAmp(1)=mag_amplitude(magind, magamp)!* psconst(3) - mag_amplitude(magind, magamp)* psconst(2)
          else
                CP%InitPower%CorrType = 0
                write(*,*) "Using interpolation TABLE"
                CP%InitPower%ScalarPowerAmp(1)= mag_psamp(magind, magamp, 3) - mag_psamp(magind, magamp, 2)
          endif
          delb = 0._dl ! Set up perturbations
          pib = 1._dl
          call cls_from_params(CP,Cltemp2,ud)
          Cltemp(:,:) = Cltemp(:,:) + Cltemp2(:,:)
!Delta-Pi
          !! The combined perturbation (delb = 1, pib = 1) !yun, 2*mag_psamp(magind, magamp, 2)
          !CP%InitPower%ScalarPowerAmp(1)= mag_psamp(magind, magamp, 2)
          if(magind .ge. -1.5) then
                CP%InitPower%CorrType = 3
                write(*,*) "Using FITTING Functions"
                CP%InitPower%ScalarPowerAmp(1)= mag_amplitude(magind, magamp)* psconst(2)
          else
                CP%InitPower%CorrType = 0
                write(*,*) "Using interpolation TABLE"
                CP%InitPower%ScalarPowerAmp(1)= mag_psamp(magind, magamp, 2)
          endif
          delb = 1._dl ! Set up perturbations
          pib = 1._dl
          call cls_from_params(CP,Cltemp2,ud)
          Cltemp(:,:) = Cltemp(:,:) + Cltemp2(:,:)
!pause
          write (*,*) "End"
          ! Copy back into Cl_scalar array
          Cl_scalar(lmin:CP%Max_l, 1, C_Temp:C_last) = Cltemp(lmin:CP%Max_l, C_Temp:C_last)
          !print*, Cl_scalar(lmin:CP%Max_l, 1, C_Temp:C_last)!Cltemp(lmin:CP%Max_l, C_Temp:C_last)
          ! Reset the Alfven velocity related stuff
          call mag_reset_alfven

          return

       else if(magmode == mag_passive) then
          !amp = mag_psamp(magind, magamp, 3) * (Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2
          !camb_ind = 1._dl + 2._dl*(3 + magind)
          CP%Scalar_initial_condition = 1
          
          
          delb = 0._dl ! Set up perturbations
          pib = 0._dl
          
          CP%WantScalars = .true.
        if(magind.ge.-1.5) then
            if (Feedbacklevel>1) write(*,*) "Using FITTING FUNCTIONS"
                CP%InitPower%an(1) = magind
                CP%InitPower%CorrType = 5
                CP%InitPower%ScalarPowerAmp(1)=mag_amplitude(magind, magamp)*psconst(3)*&
                                            (Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2

        else
          if (Feedbacklevel>1) write(*,*) "Using INTERPOLATION TABLE"
          CP%InitPower%an(1) = 1._dl + 2._dl*(3 + magind)
          CP%InitPower%ScalarPowerAmp(1)= mag_psamp(magind, magamp, 3) * (Rg*(lrat + (5._dl/(8._dl*Rv) - 1)))**2
          CP%InitPower%CorrType = 0
        end if

       end if

    else if(CP%WantVectors) then 

       if(magmode == mag_compensated) then
            !Choose wheter to use the fitting functions or not
            if(magind .ge. -1.5d0) then
                if (Feedbacklevel>1) write(*,*) "USING FITTING FUNCTION"
                CP%InitPower%CorrType = 4 !VECTOR
                !CP%InitPower%an(1) = magind
                !amp = mag_amplitude(magind, magamp)* psconst(4)
                !write(*,*)
                CP%InitPower%ScalarPowerAmp(1)= mag_amplitude(magind, magamp)* psconst(4)
                CP%InitPower%an(1) = magind
                !write(*,*) "P.S. Amplitude: ", amp
            else
                if (Feedbacklevel>1) write(*,*) "Using Table for integrals"
                CP%InitPower%CorrType = 0 ! use the approximated results.
                !amp = mag_psamp(magind, magamp, 4)
                !camb_ind = 1._dl + 2._dl*(3 + magind)
                CP%InitPower%an(1) = 1._dl + 2._dl*(3 + magind)
                CP%InitPower%ScalarPowerAmp(1)= mag_psamp(magind, magamp, 4)
            end if
            vec_sig0 = 0._dl
            delb = 0._dl ! Set up perturbations
            pib = 1._dl
            CP%WantVectors = .true.
       else
          write (*,*) "There are no passive vector modes."
          return
       end if
    end if       !yun



    ! A little debug output
    if(DebugMsgs .and. Feedbacklevel > 0) then
       write(*,*) "\t=== Mag debug ==="
       write(*,*) "\tMag amp: ", magamp, "Pert amp: ", CP%InitPower%ScalarPowerAmp(1)
       if(CP%WantTensors) then
          write(*,*) "\tMag ind: ", magind, "CAMB ind: ", CP%InitPower%ant(1)
       else
          write(*,*) "\tMag ind: ", magind, "CAMB ind: ", CP%InitPower%an(1)
       end if
       write(*,*) "\tMag prd: ", maglrat
    end if

    ! Fetch results into temporary array. Generally we'll pull the
    ! results out of the CAMB Cl_scalar etc. variables
    call cls_from_params(CP,Cltemp,ud)

    ! Set mag perturbations back to zero.
    delb = 0._dl
    pib = 0._dl

    ! Reset the Alfven velocity related stuff
    call mag_reset_alfven

  end subroutine mag_cls

end module magnetic


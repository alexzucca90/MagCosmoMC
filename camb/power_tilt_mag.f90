!
!   MODULE FOR MAGNETIC POWER SPECTRA
!
!   Alex Zucca, azucca@sfu.ca
! 
!   version 1.0, introducing fits.
!


    !This module provides the initial power spectra, parameterized as an expansion in ln k
    !
    ! ln P_s = ln A_s + (n_s -1)*ln(k/k_0_scalar) + n_{run}/2 * ln(k/k_0_scalar)^2 + n_{runrun}/6 * ln(k/k_0_scalar)^3
    !
    ! so if n_{run} = 0, n_{runrun}=0
    !
    ! P_s = A_s (k/k_0_scalar)^(n_s-1)
    !
    !for the scalar spectrum, when n_s=an(in) is the in'th spectral index. k_0_scalar
    !is a pivot scale, fixed here to 0.05/Mpc (change it below as desired or via .ini file).
    !
    !The tensor spectrum has three different supported parameterizations giving
    !
    ! ln P_t = ln A_t + n_t*ln(k/k_0_tensor) + n_{t,run}/2 * ln(k/k_0_tensor)^2
    !
    ! tensor_parameterization==tensor_param_indeptilt (=1) (default, same as CAMB pre-April 2014)
    !
    ! A_t = r A_s
    !
    ! tensor_parameterization==tensor_param_rpivot (=2)
    !
    ! A_t = r P_s(k_0_tensor)
    !
    ! tensor_parameterization==tensor_param_AT (=3)
    !
    ! A_t =  tensor_amp
    !
    !The absolute normalization of the Cls is unimportant here, but the relative ratio
    !of the tensor and scalar Cls generated with this module will be correct for general models
    !
    !December 2003 - changed default tensor pivot to 0.05 (consistent with CMBFAST 4.5)
    !April 2014 added different tensor parameterizations, running of running and running of tensors

    module InitialPower
!------
!ALEX:
!Adding using magnetic, I'll need some parameters.
!use magnetic
!------
    use Precision
    use AMLutils
    implicit none

    private

    character(LEN=*), parameter :: Power_Name = 'power_tilt'

    integer, parameter :: nnmax= 5
    !Maximum possible number of different power spectra to use

    integer, parameter, public :: tensor_param_indeptilt=1,  tensor_param_rpivot = 2, tensor_param_AT = 3



    Type InitialPowerParams
        integer :: tensor_parameterization = tensor_param_indeptilt
        integer nn  !Must have this variable
        !The actual number of power spectra to use

        !For the default implementation return power spectra based on spectral indices
        real(dl) an(nnmax) !scalar spectral indices
        real(dl) n_run(nnmax) !running of spectral index
        real(dl) n_runrun(nnmax) !running of spectral index
        real(dl) ant(nnmax) !tensor spectral indices
        real(dl) nt_run(nnmax) !tensor spectral index running
        real(dl) rat(nnmax) !ratio of scalar to tensor initial power spectrum amplitudes
        real(dl) k_0_scalar, k_0_tensor !pivot scales
        real(dl) ScalarPowerAmp(nnmax)
        real(dl) TensorPowerAmp(nnmax) !A_T at k_0_tensor if tensor_parameterization==tensor_param_AT
        !MagCAMB: adding two new parameters..
        real(dl) :: kDissip ! k_D for magnetic fields...
        integer :: CorrType = 0 ! Default 0 to run standard code....
    end Type InitialPowerParams

    real(dl) curv  !Curvature contant, set in InitializePowers

    Type(InitialPowerParams), save :: P

!==================================================
!AZ:  COEFFICIENTS FOR FITTING FUNCTIONS
!-------------- SCALAR ---------------
!----Delta-Delta----
!n>0
!General model:
!f(x,y) = 4/(2*x+3)*y^3-y^4+(sdd1+sdd2*x+sdd3*x^2+sdd4*x^3)*y^5
!Coefficients (with 95% confidence bounds):
real(dl) :: sdd1 =      0.4709  !(-14.64, 15.58)
real(dl) :: sdd2 =      0.2305  !(-32.63, 33.09)
real(dl) :: sdd3 =      0.8443  !(-22.03, 23.72)
real(dl) :: sdd4 =      0.1948  !(-4.701, 5.09)
!Goodness of fit:
!SSE: 3.106e-18
!R-square: 1
!Adjusted R-square: 1
!RMSE: 1.643e-10

!-1.5<n<0
!Linear model Poly35:
!f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y
!+ p12*x*y^2 + p03*y^3 + p31*x^3*y + p22*x^2*y^2 + p13*x*y^3
!+ p04*y^4 + p32*x^3*y^2 + p23*x^2*y^3 + p14*x*y^4 +
!p05*y^5
!Coefficients (with 95% confidence bounds):
real(dl) :: ssdd00 =    -0.01177  !(-58.51, 58.48)
real(dl) :: ssdd10 =     0.01639  !(-17.79, 17.82)
real(dl) :: ssdd01 =       2.774  !(-46.82, 52.37)
real(dl) :: ssdd20 =       1.503  !(-5.841, 8.848)
real(dl) :: ssdd11 =      0.3287  !(-10.91, 11.56)
real(dl) :: ssdd02 =    -0.07678  !(-16.82, 16.67)
real(dl) :: ssdd30 =       0.675  !(-1.329, 2.679)
real(dl) :: ssdd21 =      0.6162  !(-2.619, 3.851)
real(dl) :: ssdd12 =    0.006669  !(-2.685, 2.699)
real(dl) :: ssdd03 =    -0.01118  !(-2.823, 2.801)
real(dl) :: ssdd31 =      0.4391  !(-0.2481, 1.126)
real(dl) :: ssdd22 =   -0.001858  !(-0.4843, 0.4806)
real(dl) :: ssdd13 =    0.001018  !(-0.2894, 0.2914)
real(dl) :: ssdd04 =  -0.0008377  !(-0.2358, 0.2342)
real(dl) :: ssdd32 =     0.01505  !(-0.04276, 0.07286)
real(dl) :: ssdd23 =   -0.001302  !(-0.02672, 0.02411)
real(dl) :: ssdd14 =   9.271e-05  !(-0.0118, 0.01199)
real(dl) :: ssdd05 =  -2.614e-05  !(-0.007839, 0.007787)
!Goodness of fit:
!SSE: 0.05841
!R-square: 1
!Adjusted R-square: 1
!RMSE: 0.01554

!---- Pi-Pi ----
!n>=0
!General model:
!f(x,y) = 28/(5*(2*x+3))*y^3-y^4+(spp1+spp2*x+spp3*x^2+spp4*x^3)*y^5
!Coefficients (with 95% confidence bounds):
real(dl) :: spp1 =      0.2625  !(-4.327, 4.852)
real(dl) :: spp2 =      0.5079  !(-16.67, 17.68)
real(dl) :: spp3 =      0.2967  !(-15.58, 16.17)
real(dl) :: spp4 =      0.7112  !(-3.294, 4.717)
!Goodness of fit:
!SSE: 8.409e-18
!R-square: 1
!Adjusted R-square: 1
!RMSE: 2.728e-10

!-1.5<=n<0
!Linear model Poly35:
!f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y
!+ p12*x*y^2 + p03*y^3 + p31*x^3*y + p22*x^2*y^2 + p13*x*y^3
!+ p04*y^4 + p32*x^3*y^2 + p23*x^2*y^3 + p14*x*y^4 +
!p05*y^5
!Coefficients (with 95% confidence bounds):
real(dl) :: sspp00 =       5.613  !(-428.3, 439.5)
real(dl) :: sspp10 =       3.014  !(-85.5, 91.53)
real(dl) :: sspp01 =       7.655  !(-387.8, 403.1)
real(dl) :: sspp20 =       7.482  !(-17.64, 32.6)
real(dl) :: sspp11 =      0.8361  !(-61.13, 62.8)
real(dl) :: sspp02 =       1.827  !(-142.1, 145.7)
real(dl) :: sspp30 =       2.371  !(-3.049, 7.791)
real(dl) :: sspp21 =       3.199  !(-9.33, 15.73)
real(dl) :: sspp12 =     -0.2467  !(-16.65, 16.16)
real(dl) :: sspp03 =      0.3756  !(-25.75, 26.5)
real(dl) :: sspp31 =       1.096  !(-0.8987, 3.09)
real(dl) :: sspp22 =      0.3132  !(-1.815, 2.441)
real(dl) :: sspp13 =    -0.06395  !(-2.011, 1.883)
real(dl) :: sspp04 =     0.03923  !(-2.327, 2.405)
real(dl) :: sspp32 =     0.07151  !(-0.1103, 0.2534)
real(dl) :: sspp23 =    0.009197  !(-0.116, 0.1344)
real(dl) :: sspp14 =   -0.003868  !(-0.09127, 0.08353)
real(dl) :: sspp05 =    0.001642  !(-0.08389, 0.08718)
!Goodness of fit:
!SSE: 0.1325
!R-square: 0.9999
!Adjusted R-square: 0.9999
!RMSE: 0.02293

!---- Delta-Pi -----
!n>=0
!General model:
!f(x,y) = (1/4)*y^4+(-8/15+sdp1*x)*y^5
!Coefficients (with 95% confidence bounds):
real(dl) :: sdp1 =      0.3922  !(0.3169, 0.4676)
!Goodness of fit:
!SSE: 5.597e-20
!R-square: 0.9991
!Adjusted R-square: 0.9991
!RMSE: 2.169e-11

!-1.5<n<0
!General model:
!f(x,y) = (sdp00+sdp01*n+sdp02*n2+sdp03*n3)*t**(-2*n-2 )+(sdp10+sdp11*n+
!sdp12*n2+sdp13*n3)*t**(-2*n-1)
!Coefficients (with 95% confidence bounds):
real(dl) :: ssdp00 =      0.2294  !(0.06578, 0.393)
real(dl) :: ssdp01 =       0.378  !(-1.331, 2.087)
real(dl) :: ssdp02 =       2.734  !(-2.905, 8.374)
real(dl) :: ssdp03 =       3.701  !(-2.205, 9.607)
real(dl) :: ssdp10 =       51.22  !(-35.78, 138.2)
real(dl) :: ssdp11 =       304.3  !(-573.5, 1182)
real(dl) :: ssdp12 =       380.4  !(-2388, 3149)
real(dl) :: ssdp13 =       28.88  !(-2715, 2772)
!Goodness of fit:
!SSE: 5.117e+05
!R-square: 0.9998
!Adjusted R-square: 0.9998
!RMSE: 41.86


!--------------- VECTOR --------------
!n>=0
real(dl) :: vv1 =      0.5533
real(dl) :: vv2 =     -0.2252
real(dl) :: vv3 =     0.04957
real(dl) :: vv4 =    -0.00431
real(dl) :: vv5 =    -0.06265
real(dl) :: vv6 =      0.1975
real(dl) :: vv7 =    -0.02277
real(dl) :: vv8 =    0.001899

! This must be replaced
!-1.5 < n < 0

!Linear model Poly45:
!f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y
!+ p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y + p22*x^2*y^2
!+ p13*x*y^3 + p04*y^4 + p41*x^4*y + p32*x^3*y^2 + p23*x^2*y^3
!+ p14*x*y^4 + p05*y^5
!Coefficients (with 95% confidence bounds):
real(dl) :: V00 =       183.9  !(-456.6, 824.5)
real(dl) :: V10 =      -3.996  !(-109.3, 101.3)
real(dl) :: V01 =       178.2  !(-430.5, 786.8)
real(dl) :: V20 =       -3.53  !(-26.79, 19.73)
real(dl) :: V11 =     -0.5511  !(-78.78, 77.68)
real(dl) :: V02 =       66.16  !(-164.9, 297.2)
real(dl) :: V30 =      0.3934  !(-4.924, 5.711)
real(dl) :: V21 =      -3.349  !(-15.73, 9.036)
real(dl) :: V12 =      0.7013  !(-21.19, 22.6)
real(dl) :: V03 =        12.4  !(-31.41, 56.2)
real(dl) :: V40 =      -2.622  !(-3.631, -1.613)
real(dl) :: V31 =       1.472  !(-0.2385, 3.182)
real(dl) :: V22 =      -1.108  !(-3.364, 1.149)
real(dl) :: V13 =      0.2548  !(-2.482, 2.992)
real(dl) :: V04 =        1.15  !(-2.997, 5.297)
real(dl) :: V41 =     -0.7025  !(-0.8977, -0.5073)
real(dl) :: V32 =      0.3475  !(0.1969, 0.498)
real(dl) :: V23 =     -0.1292  !(-0.2696, 0.01115)
real(dl) :: V14 =     0.02458  !(-0.1043, 0.1535)
real(dl) :: V05 =     0.04212  !(-0.1147, 0.199)

!Goodness of fit:
!SSE: 0.0125
!R-square: 1
!Adjusted R-square: 1
!RMSE: 0.008155



!------------ TENSOR ---------------
!n>0
real(dl) :: tt1 =      0.5532
real(dl) :: tt2 =     -0.2252
real(dl) :: tt3 =     0.04958
real(dl) :: tt4 =    -0.00431
real(dl) :: tt5 =     0.08194
real(dl) :: tt6 =      0.2398
real(dl) :: tt7 =     -0.0129
real(dl) :: tt8 =   0.0006638


!-1.5<n<0
!Linear model Poly55:
!     f(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y 
!                    + p12*x*y^2 + p03*y^3 + p40*x^4 + p31*x^3*y + p22*x^2*y^2 
!                    + p13*x*y^3 + p04*y^4 + p50*x^5 + p41*x^4*y + p32*x^3*y^2 
!                    + p23*x^2*y^3 + p14*x*y^4 + p05*y^5
!Coefficients (with 95% confidence bounds):
real(dl) ::        p00 =      -1.374  !(-107.2, 104.4)
real(dl) ::        p10 =       -3.13  !(-19.44, 13.18)
real(dl) ::        p01 =       2.485  !(-98.03, 103)
real(dl) ::        p20 =      -6.108  !(-9.6, -2.616)
real(dl) ::        p11 =     -0.5194  !(-12.63, 11.59)
real(dl) ::        p02 =     -0.1435  !(-38.29, 38)
real(dl) ::        p30 =       -7.28  !(-8.098, -6.462)
real(dl) ::        p21 =      -1.188  !(-3.031, 0.6541)
real(dl) ::        p12 =    -0.01076  !(-3.398, 3.377)
real(dl) ::        p03 =    -0.02336  !(-7.254, 7.208)
real(dl) ::        p40 =      -3.398  !(-3.598, -3.198)
real(dl) ::        p31 =      -1.197  !(-1.454, -0.9396)
real(dl) ::        p22 =   -0.006893  !(-0.3404, 0.3266)
real(dl) ::        p13 =   1.793e-05  !(-0.4233, 0.4233)
real(dl) ::        p04 =   -0.001978  !(-0.6865, 0.6825)
real(dl) ::        p50 =     -0.4739  !(-0.508, -0.4398)
real(dl) ::        p41 =     -0.4924  !(-0.5187, -0.466)
real(dl) ::        p32 =     0.01338  !(-0.009523, 0.03628)
real(dl) ::        p23 =   -0.001639  !(-0.02234, 0.01906)
real(dl) ::        p14 =   0.0001047  !(-0.01983, 0.02004)
real(dl) ::        p05 =  -6.991e-05  !(-0.02596, 0.02582)
!Goodness of fit:
!  SSE: 0.0004587
!  R-square: 1
!  Adjusted R-square: 1
!  RMSE: 0.001447

! END OF COEFFICIENTS FOR FITTING FUNCTIONS
!==================================================================

! psconst needed only for scalar compensated modes.
real(dl), parameter :: psconst(5) = (/0.25, 0.25, 0.25, 2., 2./3. /)
! kDissip and CorrType can be moved inside the P.S parameters.
!real(dl) :: kDissip
!integer :: CorrType = 0 ! In order to run without magnetic fields
!real(dl) :: Magnetic_Index
!ALEX.


    !Make things visible as neccessary...

    public InitialPowerParams, InitialPower_ReadParams, InitializePowers, ScalarPower, TensorPower
    public nnmax,Power_Descript, Power_Name, SetDefPowerParams
!ALEX:
!public kDissip, CorrType, Magnetic_Index
!ALEX.

    contains


    subroutine SetDefPowerParams(AP)
    Type (InitialPowerParams) :: AP

    AP%nn     = 1 !number of initial power spectra
    AP%an     = 1 !scalar spectral index
    AP%n_run   = 0 !running of scalar spectral index
    AP%n_runrun   = 0 !running of running of scalar spectral index
    AP%ant    = 0 !tensor spectra index
    AP%nt_run   = 0 !running of tensor spectral index
    AP%rat    = 1
    AP%k_0_scalar = 0.05
    AP%k_0_tensor = 0.05
    AP%ScalarPowerAmp = 1
    AP%TensorPowerAmp = 1
    AP%tensor_parameterization = tensor_param_indeptilt
    end subroutine SetDefPowerParams

    subroutine InitializePowers(AParamSet,acurv)
    Type (InitialPowerParams) :: AParamSet
    !Called before computing final Cls in cmbmain.f90
    !Could read spectra from disk here, do other processing, etc.

    real(dl) acurv

    if (AParamSet%nn > nnmax) then
        write (*,*) 'To use ',AParamSet%nn,'power spectra you need to increase'
        write (*,*) 'nnmax in power_tilt.f90, currently ',nnmax
    end if
    P = AParamSet

    curv=acurv

    !Write implementation specific code here...

    end subroutine InitializePowers


    function ScalarPower(k,ix)

    !"ix" gives the index of the power to return for this k
    !ScalarPower = const for scale invariant spectrum
    !The normalization is defined so that for adiabatic perturbations the gradient of the 3-Ricci
    !scalar on co-moving hypersurfaces receives power
    ! < |D_a R^{(3)}|^2 > = int dk/k 16 k^6/S^6 (1-3K/k^2)^2 ScalarPower(k)
    !In other words ScalarPower is the power spectrum of the conserved curvature perturbation given by
    !-chi = \Phi + 2/3*\Omega^{-1} \frac{H^{-1}\Phi' - \Psi}{1+w}
    !(w=p/\rho), so < |\chi(x)|^2 > = \int dk/k ScalarPower(k).
    !Near the end of inflation chi is equal to 3/2 Psi.
    !Here nu^2 = (k^2 + curv)/|curv|

    !This power spectrum is also used for isocurvature modes where
    !< |\Delta(x)|^2 > = \int dk/k ScalarPower(k)
    !For the isocurvture velocity mode ScalarPower is the power in the neutrino heat flux.

    real(dl) ScalarPower,k, lnrat, kTilde
    integer ix

! MagCAMB: adding the parameters for fitting the functions
!integer  !This gives me the type of correlation function that has to be computed
real(dl) :: t, t2, t3, t4, t5!, t6, t7, t8, t9, t10, 
real(dl) :: lnt,lnt2,lnt3,lnt4,lnt5
real(dl) :: n,n2, n3, n4, n5
if(P%CorrType /= 0) then
    t = k/P%kDissip
    t2 = t*t
    t3 = t2*t
    t4 = t3*t
    t5 = t4*t
    lnt = log(t)
    lnt2 = lnt*lnt
    lnt3 = lnt2*lnt
    lnt4 = lnt3*lnt
    lnt5 = lnt4*lnt
    !n = Magnetic_Index
    n= P%an(ix)
    n2 = n*n
    n3 = n2*n
    n4 = n3*n
    n5 = n4*n
end if

! INFLATIONARY POWER SPECTRUM
if(P%CorrType==0) then
    lnrat = log(k/P%k_0_scalar)
    ScalarPower=P%ScalarPowerAmp(ix)*exp(lnrat*( P%an(ix)-1 + lnrat*(P%n_run(ix)/2 + P%n_runrun(ix)/6*lnrat)))

! SCALAR MAGNETIC
!Compensated
else if(P%CorrType==1) then!Delta-Delta Fit
    if(n .ge. 0) then !n>=0
        if (t .le. 0.5) then
            !write(*,*) "CorrType = ", CorrType
            !write(*,*) "t = ", t
            ScalarPower = P%ScalarPowerAmp(ix) * (P%kDissip/ (2*pi))**(2.d0*n+6)*&
                    (psconst(1)*(4/(2*n+3)*t3-t4+(sdd1+sdd2*n+sdd3*n2+sdd4*n3)*t5)&
                    -psconst(2)*((1/4)*t4+(-8/15+sdp1*n)*t5))
            !write(*,*) "P(k) = ", ScalarPower
            !write(*,*) "Sc. Amp = ", P%ScalarPowerAmp(ix)
            !ScalarPower = 0.d0
        else
            write(*,*) "t = ", t, "Fitting functions not computed for t>0.5"
            write(*,*) "Set P(k)=0"
            ScalarPower = 0.d0
        end if
    else ! -1.5<n<0
        if(t .le. 0.5) then !t<0.5
            ScalarPower = P%ScalarPowerAmp(ix) * (P%kDissip/ (2*pi))**(2.d0*n+6)*&
            (psconst(1)*exp(ssdd00 + ssdd10*n + ssdd01*lnt + ssdd20*n2 + ssdd11*n*lnt &
            + ssdd02*lnt2 + ssdd30*n3 + ssdd21*n2*lnt + ssdd12*n*lnt2 + ssdd03*lnt3 + &
            ssdd31*n3*lnt + ssdd22*n2*lnt2 + ssdd13*n*lnt3 + ssdd04*lnt4 + ssdd32*n3*lnt2&
            + ssdd23*n2*lnt3 + ssdd14*n*lnt4 + ssdd05*lnt5) - psconst(2)*(t**(2*n+6)*&
            ((ssdp00+ssdp01*n+ssdp02*n2+ssdp03*n3)*t**(-2*n-2 )+(ssdp10+ssdp11*n+ssdp12*n2+ssdp13*n3)*t**(-2*n-1))))


            !exp(ssdp00 + &
            !ssdp10*n + ssdp01*lnt + ssdp20*n2 + ssdp11*n*lnt + ssdp02*lnt2 + ssdp30*n3 + ssdp21*n2*lnt &
            !+ ssdp12*n*lnt2 + ssdp03*lnt3 + ssdp40*n4 + ssdp31*n3*lnt + ssdp22*n2*lnt2 &
            !+ ssdp13*n*lnt3 + ssdp04*lnt4 + ssdp50*n5 + ssdp41*n4*lnt + ssdp32*n3*lnt2 &
            !+ ssdp23*n2*lnt3 + ssdp14*n*lnt4 + ssdp05*lnt5))
            !ScalarPower = 0.d0
        else !t>0.5
            write(*,*) "t = ", t, "Fitting functions not computed for t>0.5"
            write(*,*) "Set P(k)=0"
            ScalarPower = 0.d0
        end if
    end if
else if(P%CorrType==2) then!Pi - Pi
    if(n .ge. 0) then !n>=0
        if (t .le. 0.5) then
            !write(*,*) "CorrType = ", CorrType
            !write(*,*) "t = ", t
            ScalarPower = P%ScalarPowerAmp(ix) * (P%kDissip/ (2*pi))**(2.d0*n+6)*&
            (psconst(3)*(28/(5*(2*n+3))*t3-t4+(spp1+spp2*n+spp3*n2+spp4*n3)*t5)&
            -psconst(2)*((1/4)*t4+(-8/15+sdp1*n)*t5))
            !write(*,*) "P(k) = ", ScalarPower
            !write(*,*) "Sc. Amp = ", P%ScalarPowerAmp(ix)
            !ScalarPower = 0.d0
        else
            write(*,*) "t = ", t, "Fitting functions not computed for t>0.5"
            write(*,*) "Set P(k)=0"
            ScalarPower = 0.d0
        end if
    else ! -1.5<n<0
        if(t .le. 0.5) then !t<0.5
            ScalarPower = P%ScalarPowerAmp(ix) * (P%kDissip/ (2*pi))**(2.d0*n+6)*&
            (psconst(3)*exp(sspp00 + sspp10*n + sspp01*lnt + sspp20*n2 + sspp11*n*lnt &
            + sspp02*lnt2 + sspp30*n3 + sspp21*n2*lnt + sspp12*n*lnt2 + sspp03*lnt3 + &
            sspp31*n3*lnt + sspp22*n2*lnt2 + sspp13*n*lnt3 + sspp04*lnt4 + sspp32*n3*lnt2&
            + sspp23*n2*lnt3 + sspp14*n*lnt4 + sspp05*lnt5) -psconst(2)*(t**(2*n+6)*&
            ((ssdp00+ssdp01*n+ssdp02*n2+ssdp03*n3)*t**(-2*n-2 )+(ssdp10+ssdp11*n+ssdp12*n2+ssdp13*n3)*t**(-2*n-1))))


!psconst(2)*exp(ssdp00 + &
!            ssdp10*n + ssdp01*lnt + ssdp20*n2 + ssdp11*n*lnt + ssdp02*lnt2 + ssdp30*n3 + ssdp21*n2*lnt &
!            + ssdp12*n*lnt2 + ssdp03*lnt3 + ssdp40*n4 + ssdp31*n3*lnt + ssdp22*n2*lnt2 &
!            + ssdp13*n*lnt3 + ssdp04*lnt4 + ssdp50*n5 + ssdp41*n4*lnt + ssdp32*n3*lnt2 &
!            + ssdp23*n2*lnt3 + ssdp14*n*lnt4 + ssdp05*lnt5))
            !ScalarPower=0.d0
        else !t>0.5
            write(*,*) "t = ", t, "Fitting functions not computed for t>0.5"
            write(*,*) "Set P(k)=0"
            ScalarPower = 0.d0
        end if
    end if
else if(P%CorrType==3) then!Delta - Pi
    if(n .ge. 0) then !n>=0
        if (t .le. 0.5) then
            !write(*,*) "CorrType = ", CorrType
            !write(*,*) "t = ", t
            ScalarPower = P%ScalarPowerAmp(ix) * (P%kDissip/ (2*pi))**(2.d0*n+6)*(&
            (1/4)*t4+(-8/15+sdp1*n)*t5)
            !ScalarPower = 0.d0
            !write(*,*) "P(k) = ", ScalarPower
            !write(*,*) "Sc. Amp = ", P%ScalarPowerAmp(ix)
        else
            write(*,*) "t = ", t, "Fitting functions not computed for t>0.5"
            write(*,*) "Set P(k)=0"
            ScalarPower = 0.d0
        end if
    else !-1.5<n<0
        if(t .le. 0.5) then !t<0.5
            ScalarPower = P%ScalarPowerAmp(ix) * (P%kDissip/ (2*pi))**(2.d0*n+6)*&
            t**(2*n+6)*((ssdp00+ssdp01*n+ssdp02*n2+ssdp03*n3)*t**(-2*n-2 )+&
            (ssdp10+ssdp11*n+ssdp12*n2+ssdp13*n3)*t**(-2*n-1))

!            exp(ssdp00 + ssdp10*n + ssdp01*lnt + ssdp20*n2 + ssdp11*n*lnt + &
!            ssdp02*lnt2 + ssdp30*n3 + ssdp21*n2*lnt + ssdp12*n*lnt2 + ssdp03*lnt3 &
!            + ssdp40*n4 + ssdp31*n3*lnt + ssdp22*n2*lnt2 + ssdp13*n*lnt3 + ssdp04*lnt4&
!            + ssdp50*n5 + ssdp41*n4*lnt + ssdp32*n3*lnt2 + ssdp23*n2*lnt3 + ssdp14*n*lnt4 + ssdp05*lnt5)
            !write(*,*) "Amp. = ", P%ScalarPowerAmp(ix)
            !write(*,*) "Delta-Pi P(k) = ", ScalarPower
        else !t>0.5
            write(*,*) "t = ", t, "Fitting functions not computed for t>0.5"
            write(*,*) "Set P(k)=0"
            ScalarPower = 0.d0
        end if
    end if
! Scalar Passive modes
else if(P%CorrType==5) then
    !N>0
    if(n .ge. 0) then
        ScalarPower = P%ScalarPowerAmp(ix) * (P%kDissip/ (2*pi))**(2.d0*n+6)*&
            (28/(5*(2*n+3))*t3-t4+(spp1+spp2*n+spp3*n2+spp4*n3)*t5)
    else !-1.5<n<0
        if(t .le. 0.5) then !t<0.5
            ScalarPower = P%ScalarPowerAmp(ix) * (P%kDissip/ (2*pi))**(2.d0*n+6)*&
            exp(sspp00 + sspp10*n + sspp01*lnt + sspp20*n2 + sspp11*n*lnt &
            + sspp02*lnt2 + sspp30*n3 + sspp21*n2*lnt + sspp12*n*lnt2 + sspp03*lnt3 + &
            sspp31*n3*lnt + sspp22*n2*lnt2 + sspp13*n*lnt3 + sspp04*lnt4 + sspp32*n3*lnt2&
            + sspp23*n2*lnt3 + sspp14*n*lnt4 + sspp05*lnt5)
        else !t>0.5
            write(*,*) "t = ", t, "Fitting functions not computed for t>0.5"
            write(*,*) "Set P(k)=0"
            ScalarPower = 0.d0
        end if
    end if

! VECTOR MAGNETIC
! Only Compensated Modes.
else if(P%CorrType==4) then ! Pi-Pi
    if (n .ge. 0) then
        !write(*,*) "using exact p.s. Amplitude: ", P%ScalarPowerAmp(ix)

        if(t<=0.5) then
            !FITTING FUNCTION
            ScalarPower = P%ScalarPowerAmp(ix) * (P%kDissip/ (2*pi))**(2.d0*n+6)*&
                    ((vv1+vv2*n+vv3*n2+vv4*n3)*t3-&
                    5.d0/12.d0*t4+&
                    (vv5+vv6*n+vv7*n2+vv8*n3)*t5)
        else
            write(*,*) "t = ", t, "Fitting functions not computed for t>0.5"
            write(*,*) "Set P(k)=0"
            ScalarPower = 0.d0
        end if
    else if ((n < 0) .and.( n .ge. -1.5)) then!  -1.5 < n < 0
        if (t <= 0.5 ) then
            ScalarPower =P%ScalarPowerAmp(ix) * (P%kDissip/ (2*pi))**(2.d0*n+6)*exp(&
            V00 + V10*n + V01*lnt + V20*n2 + V11*n*lnt + V02*lnt2 + V30*n3 + V21*n2*lnt &
            + V12*n*lnt2 + V03*lnt3 + V40*n4 + V31*n3*lnt + V22*n2*lnt2 &
            + V13*n*lnt3 + V04*lnt4 + V41*n4*lnt + V32*n3*lnt2 + V23*n2*lnt3 &
            + V14*n*lnt4 + V05*lnt5 )
        else
            write(*,*) "t = ", t, " t>0.5 , P(k) = 0"
            ScalarPower = 0.d0
        end if
    ! This was for testing.....
    !else if (n<-1.5 .and. n> -1.75) then
    !    if(t<=0.5) then
    !        ScalarPower = P%ScalarPowerAmp(ix) * (kDissip/ (2*pi))**(2.d0*n+6)* t**(2.d0*n+6)*&
    !                (p00 + p10*n + p01*t + p20*n2 + p11*n*t + p02*t2 + p30*n3 + p21*n2*t &
    !                + p12*n*t2 + p03*t3 + p40*n4 + p31*n3*t + p22*n2*t2 &
    !                + p13*n*t3 + p04*t4 + p50*n5 + p41*n4*t + p32*n3*t2 &
    !                + p23*n2*t3 + p14*n*t4 + p05*t5)
!if(ScalarPower < 0) write(*,*) "PS < 0!!!!", ScalarPower
        !else
        !    write(*,*) "t = ", t, " t>0.5 , P(k) = 0"
        !    ScalarPower = 0.d0
        !end if
    end if
else
    !ERROR!!!
    write(*,*) "Power Spectrum not set..."
    ! I need to put a flag to kill the process.
end if

    !         ScalarPower = ScalarPower * (1 + 0.1*cos( lnrat*30 ) )

    end function ScalarPower


    function TensorPower(k,ix)

    !TensorPower= const for scale invariant spectrum
    !The normalization is defined so that
    ! < h_{ij}(x) h^{ij}(x) > = \sum_nu nu /(nu^2-1) (nu^2-4)/nu^2 TensorPower(k)
    !for a closed model
    ! < h_{ij}(x) h^{ij}(x) > = int d nu /(nu^2+1) (nu^2+4)/nu^2 TensorPower(k)
    !for an open model
    !"ix" gives the index of the power spectrum to return
    !Here nu^2 = (k^2 + 3*curv)/|curv|

    real(dl) TensorPower,k
    real(dl), parameter :: PiByTwo=3.14159265d0/2._dl
    integer ix
    real(dl) lnrat, k_dep
    !New stuff for magnetic fields
    real(dl) :: t, t2, t3, t4, t5!, t6, t7, t8, t9, t10,
    real(dl) :: lnt,lnt2,lnt3,lnt4,lnt5
    real(dl) :: n,n2, n3, n4, n5
if(P%CorrType /= 0) then
t = k/P%kDissip
    !write(*,*) "t = ", t
    t2 = t*t
    t3 = t2*t
    t4 = t3*t
    t5 = t4*t
    lnt = log(t)
    lnt2 = lnt*lnt
    lnt3 = lnt2*lnt
    lnt4 = lnt3*lnt
    lnt5 = lnt4*lnt
    !n = Magnetic_Index
    n = P%ant(ix)
    n2 = n*n
    n3 = n2*n
    n4 = n3*n
    n5 = n*n4
end if

!INFLATIONARY 
if (P%CorrType==0) then
    lnrat = log(k/P%k_0_tensor)
    k_dep = exp(lnrat*(P%ant(ix) + P%nt_run(ix)/2*lnrat))
    if (P%tensor_parameterization==tensor_param_indeptilt) then
        TensorPower = P%rat(ix)*P%ScalarPowerAmp(ix)*k_dep
    else if (P%tensor_parameterization==tensor_param_rpivot) then
        TensorPower = P%rat(ix)*ScalarPower(P%k_0_tensor,ix) * k_dep
    else if (P%tensor_parameterization==tensor_param_AT) then
        TensorPower = P%TensorPowerAmp(ix) * k_dep
    end if
    if (curv < 0) TensorPower=TensorPower*tanh(PiByTwo*sqrt(-k**2/curv-3))

!Tensor Passive
else if (P%CorrType == 5) then
    !n>0
    if (n .ge. 0) then
        if(t .le. 0.5) then
        TensorPower = P%TensorPowerAmp(ix)* (P%kDissip/ (2*pi))**(2.d0*n+6)*&
        ((tt1+tt2*n+tt3*n2+tt4*n3)*t3-&
        7.d0/12.d0*t4+&
        (tt5+tt6*n+tt7*n2+tt8*n3)*t5)
        else
            write(*,*) "t = ", t, "Fitting functions not computed for t>0.5"
            write(*,*) "Set P(k)=0"
            TensorPower = 0.d0
        end if
    else if(n<0 .and. n .ge. -1.5_dl) then
        if(t .le. 0.5) then
        !write(*,*)"computing tensor power spectrum"
        !TensorPower = P%TensorPowerAmp(ix)* (P%kDissip/ (2*pi))**(2.d0*n+6)*&
        !exp(T00 + T10*n + T01*lnt + T20*n2 + T11*n*lnt + T02*lnt2 + T30*n3 + T21*n2*lnt&
        !+ T12*n*lnt2 + T03*lnt3 + T31*n3*lnt + T22*n2*lnt2 + T13*n*lnt3&
        !+ T04*lnt4 + T32*n3*lnt2 + T23*n2*lnt3 + T14*n*lnt4 + T05*lnt5)

TensorPower = P%TensorPowerAmp(ix)* (P%kDissip/ (2*pi))**(2.d0*n+6)*&
exp(p00 + p10*n + p01*lnt + p20*n2 + p11*n*lnt + p02*lnt2 + p30*n3 + p21*n2*lnt& 
                    + p12*n*lnt2 + p03*lnt3 + p40*n4 + p31*n3*lnt + p22*n2*lnt2& 
                    + p13*n*lnt3 + p04*lnt4 + p50*n5 + p41*n4*lnt + p32*n3*lnt2& 
                    + p23*n2*lnt3 + p14*n*lnt4 + p05*lnt5)


        !write(*,*) "Tensor P(k) = ", TensorPower
        else
            write(*,*) "t = ", t, "Fitting functions not computed for t>0.5"
            write(*,*) "Set P(k)=0"
            TensorPower = 0.d0
        end if
    end if
else
!ERROR
end if
    end function TensorPower

    !Get parameters describing parameterisation (for FITS file)
    !Does not support running extensions
    function Power_Descript(in, Scal, Tens, Keys, Vals)
    character(LEN=8), intent(out) :: Keys(*)
    real(dl), intent(out) :: Vals(*)
    integer, intent(IN) :: in
    logical, intent(IN) :: Scal, Tens
    integer num, Power_Descript
    num=0
    if (Scal) then
        num=num+1
        Keys(num) = 'n_s'
        Vals(num) = P%an(in)
        num=num+1
        Keys(num) = 'n_run'
        Vals(num) = P%n_run(in)
        num=num+1
        Keys(num) = 's_pivot'
        Vals(num) = P%k_0_scalar
    end if
    if (Tens) then
        num=num+1
        Keys(num) = 'n_t'
        Vals(num) = P%ant(in)
        num=num+1
        Keys(num) = 't_pivot'
        Vals(num) = P%k_0_tensor
        if (Scal) then
            num=num+1
            Keys(num) = 'p_ratio'
            Vals(num) = P%rat(in)
        end if
    end if
    Power_Descript = num

    end  function Power_Descript

    subroutine InitialPower_ReadParams(InitPower, Ini, WantTensors)
    use IniFile
    Type(InitialPowerParams) :: InitPower
    Type(TIniFile) :: Ini
    logical, intent(in) :: WantTensors
    integer i

    InitPower%k_0_scalar = Ini_Read_Double_File(Ini,'pivot_scalar',InitPower%k_0_scalar)
    InitPower%k_0_tensor = Ini_Read_Double_File(Ini,'pivot_tensor',InitPower%k_0_tensor)
    InitPower%nn = Ini_Read_Int_File(Ini,'initial_power_num',1)
    if (InitPower%nn>nnmax) call MpiStop('Too many initial power spectra - increase nnmax in InitialPower')
    if (WantTensors) then
        InitPower%tensor_parameterization =  Ini_Read_Int_File(Ini, 'tensor_parameterization',tensor_param_indeptilt)
        if (InitPower%tensor_parameterization < tensor_param_indeptilt .or. &
            & InitPower%tensor_parameterization > tensor_param_AT) &
            & call MpiStop('InitialPower: unknown tensor_parameterization')
    end if
    InitPower%rat(:) = 1
    do i=1, InitPower%nn
        InitPower%an(i) = Ini_Read_Double_Array_File(Ini,'scalar_spectral_index', i)
        InitPower%n_run(i) = Ini_Read_Double_Array_File(Ini,'scalar_nrun',i,0._dl)
        InitPower%n_runrun(i) = Ini_Read_Double_Array_File(Ini,'scalar_nrunrun',i,0._dl)

        if (WantTensors) then
            InitPower%ant(i) = Ini_Read_Double_Array_File(Ini,'tensor_spectral_index',i)
            InitPower%nt_run(i) = Ini_Read_Double_Array_File(Ini,'tensor_nrun',i,0._dl)
            if (InitPower%tensor_parameterization == tensor_param_AT) then
                InitPower%TensorPowerAmp(i) = Ini_Read_Double_Array_File(Ini,'tensor_amp',i)
            else
                InitPower%rat(i) = Ini_Read_Double_Array_File(Ini,'initial_ratio',i)
            end if
        end if

        InitPower%ScalarPowerAmp(i) = Ini_Read_Double_Array_File(Ini,'scalar_amp',i,1._dl)
        !Always need this as may want to set tensor amplitude even if scalars not computed
    end do

    end  subroutine InitialPower_ReadParams


    end module InitialPower

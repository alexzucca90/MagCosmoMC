module richsub

  use precision

contains

  ! General routine for cubic spline interpolation (see NR)
  real(dl) function spline_val(x, xv, yv, y2, n)

    integer, intent(in) :: n
    real(dl), intent(in) :: x
    real(dl), intent(in) :: xv(n), yv(n), y2(n)

    integer :: kh,kl,kn
    real(dl) :: h,a,b,c,d

    ! Extrapolate if value is above of below interval
    if(x < xv(1)) then
       h = xv(2) - xv(1)
       a = (yv(2) - yv(1)) / h
       spline_val = (a - h * y2(2) / 6) * (x - xv(1)) + yv(1)
    else if(x > xv(n)) then
       h = xv(n) - xv(n-1)
       a = (yv(n) - yv(n-1)) / h
       spline_val = (a + h * y2(n-1) / 6) * (x - xv(n)) + yv(n)
    else
       ! Bisection to find correct interval
       kh = n
       kl = 1
       do while(kh - kl > 1)
          kn = (kh + kl) / 2
          if(xv(kn) > x) then
             kh = kn
          else
             kl = kn
          end if
       end do

       ! Set up constants (a la NR)
       h = xv(kh) - xv(kl)

       a = (xv(kh) - x) / h
       b = (x - xv(kl)) / h
       c = (a**3 - a)* h**2 / 6
       d = (b**3 - b)* h**2 / 6

       spline_val = (a*yv(kl) + b*yv(kh) + c*y2(kl) + d*y2(kh))
       
    end if
  end function spline_val

  ! A romberg integrator in reasonable f95 style.
  real(dl) function romberg_integrator(f, a, b, tol, maxit) 
    
    interface
       real(dl) function f(x)
         use precision
         real(dl), intent(in) :: x
       end function f
    end interface

    real(dl), intent(in) :: a, b, tol
    integer, intent(in) :: maxit

    integer, parameter :: minit = 5

    real(dl) :: r(0:maxit,0:maxit)
    real(dl) :: h0, h, t
    integer :: i, j
    logical :: notfinished

    notfinished = .true.
    
    h0 = b-a

    r(0,0) = 0.5 * h0 * (f(a) + f(b));

    do i=1,maxit
       h = h0 / 2**i
       t = 0._dl

       do j=1,2**(i-1)
          t = t + f(a + (2*j-1)*h)
       end do
       t = t * h
       
       r(i,0) = 0.5 * r(i-1,0) + t
       
       do j=1,i
          r(i,j) = r(i,j-1) + (r(i,j-1) - r(i-1,j-1)) / (4._dl**j - 1)
       end do

       if(i > minit .and. abs(1._dl - (r(i,i-1) / r(i,i))) < tol) then
          romberg_integrator = r(i,i)
          notfinished = .false.
          !print *, "Finshed"
          exit
       end if
    end do
    
    if(notfinished) then
       print *, "Did not converge."
       romberg_integrator = r(maxit,maxit)
    end if
    
  end function romberg_integrator


  ! A romberg integrator in reasonable f95 style.
  ! Equivalent to rombinr_obj in CAMB, passesin an extra real(dl) argument
  real(dl) function romberg_integrator_obj(y, f, a, b, tol, maxit) 
    
    real(dl), intent(in) :: y

    interface
       real(dl) function f(y, x)
         use precision
         real(dl), intent(in) :: y
         real(dl), intent(in) :: x
       end function f
    end interface

    real(dl), intent(in) :: a, b, tol
    integer, intent(in) :: maxit

    integer, parameter :: minit = 5

    real(dl) :: r(0:maxit,0:maxit)
    real(dl) :: h0, h, t, err
    integer :: i, j
    logical :: notfinished = .true.
    
    h0 = b-a

    r(0,0) = 0.5 * h0 * (f(y, a) + f(y, b));

    do i=1,maxit
       h = h0 / 2**i
       t = 0._dl

       do j=1,2**(i-1)
          t = t + f(y, a + (2*j-1)*h)
       end do
       t = t * h
       
       r(i,0) = 0.5 * r(i-1,0) + t
       
       do j=1,i
          r(i,j) = r(i,j-1) + (r(i,j-1) - r(i-1,j-1)) / (4._dl**j - 1)
       end do

       err = abs( (r(i,i) - r(i,i-1)) / r(i,i))
       if(i > minit .and. err < tol) then
          romberg_integrator_obj = r(i,i)
          notfinished = .false.
          !print *, "Finshed"
          exit
       end if
    end do
    
    if(notfinished) then
       print *, "Did not converge. Error: ", err, " Value: ", r(maxit,maxit) 
       romberg_integrator_obj = r(maxit,maxit)
    end if
    
  end function romberg_integrator_obj

  ! A romberg integrator in reasonable f95 style.
  ! Equivalent to rombinr_obj in CAMB, passes in an extra real(dl) argument
  real(dl) function romberg_integrator_arr(y, f, a, b, tol, maxit) 
    
    real(dl), intent(in), dimension(:) :: y

    interface
       real(dl) function f(y, x)
         use precision
         real(dl), intent(in), dimension(:) :: y
         real(dl), intent(in) :: x
       end function f
    end interface

    real(dl), intent(in) :: a, b, tol
    integer, intent(in) :: maxit

    integer, parameter :: minit = 5

    real(dl) :: r(0:maxit,0:maxit)
    real(dl) :: h0, h, t, err
    integer :: i, j
    logical :: notfinished = .true.
    
    h0 = b-a

    r(0,0) = 0.5 * h0 * (f(y, a) + f(y, b));

    do i=1,maxit
       h = h0 / 2**i
       t = 0._dl

       do j=1,2**(i-1)
          t = t + f(y, a + (2*j-1)*h)
       end do
       t = t * h
       
       r(i,0) = 0.5 * r(i-1,0) + t
       
       do j=1,i
          r(i,j) = r(i,j-1) + (r(i,j-1) - r(i-1,j-1)) / (4._dl**j - 1)
       end do

       err = abs( (r(i,i) - r(i,i-1)) / r(i,i))
       if(i > minit .and. err < tol) then
          romberg_integrator_arr= r(i,i)
          notfinished = .false.
          !print *, "Finshed"
          exit
       end if
    end do
    
    if(notfinished) then
       print *, "Did not converge. Error: ", err, " Value: ", r(maxit,maxit) 
       romberg_integrator_arr = r(maxit,maxit)
    end if
    
  end function romberg_integrator_arr



  
  real(dl) function trapezoidrule(f, h, n)
    ! Trapezoid integration
    real(dl), intent(in) :: f(n)
    real(dl), intent(in) :: h
    integer, intent(in) :: n
    
    integer :: i
    real(dl) :: t

    t = 0._dl
    
    do i=2,n-1
       t = t+f(i)
    end do
    t = t + (f(1) + f(n)) / 2._dl
    trapezoidrule = t*h
  end function trapezoidrule


  
  real(dl) function trapezoidrule_irregular(f, x, n)
    ! Trapezoid integration
    real(dl), intent(in) :: f(n)
    real(dl), intent(in) :: x(n)
    integer, intent(in) :: n
    
    integer :: i
    real(dl) :: t

    t = 0._dl
    
    do i=2,n
       t = t+ 0.5_dl * (f(i)+f(i-1)) * (x(i)-x(i-1))
    end do

    trapezoidrule_irregular = t
  end function trapezoidrule_irregular



end module richsub



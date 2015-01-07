module priors

      use utils1
      private betain, xinbta, dierfc

contains

!=======================================================================
! Prior distribution functions: r is a uniform deviate from the unit
!
!

! Uniform[0:1]  ->  Delta[x1]

function DeltaFunctionPrior(r,x1,x2)

      implicit none

      double precision r,x1,x2,DeltaFunctionPrior

      DeltaFunctionPrior=x1
      
end function DeltaFunctionPrior

!=======================================================================
! Uniform[0:1]  ->  Uniform[x1:x2]

function UniformPrior(r,x1,x2)

            implicit none

            double precision r,x1,x2,UniformPrior

            UniformPrior=x1+r*(x2-x1)

end function UniformPrior

!=======================================================================
! Uniform[0:1]  ->  LogUniform[x1:x2]

function LogPrior(r,x1,x2)

            implicit none

            double precision r,x1,x2,LogPrior
            double precision lx1,lx2

            if (r.le.0.0d0) then
                  LogPrior=-1.0d32
            else
                  lx1=dlog10(x1)
                  lx2=dlog10(x2)
                  LogPrior=10.d0**(lx1+r*(lx2-lx1))
            endif
      
end function LogPrior

!=======================================================================
! Uniform[0:1]  ->  Exponential[0:inf]

function ExpPrior(r)

            implicit none

            double precision r,ExpPrior
            
            ExpPrior=-log(1.d0-r)

end function ExpPrior

!=======================================================================
! Uniform[0:1]  ->  Jeffreys[x1:x2]
function JeffreysPrior(r, x1, x2)

            implicit none

            double precision r,x1,x2,JeffreysPrior

            JeffreysPrior=x1*(x2/x1)**r
            
end function JeffreysPrior

!=======================================================================
! Uniform[0:1]  ->  Mod. Jeffreys[0:x2] with break at x1
function ModJeffreysPrior(r, x1, x2)

            implicit none

            double precision r, x1, x2, ModJeffreysPrior

            ModJeffreysPrior=x1*((1.d0+x2/x1)**r - 1.d0)

end function ModJeffreysPrior

!=======================================================================
! Uniform[0:1]  ->  Sin[x1:x2]  (angles in degrees):

function SinPrior(r,x1,x2)

            implicit none

            double precision r,x1,x2,SinPrior
            real cx1,cx2,deg2rad
            parameter(deg2rad=0.017453292)

            cx1=cos(x1*deg2rad)
            cx2=cos(x2*deg2rad)
            SinPrior=1.d0*acos(cx1+r*(cx2-cx1))

end function SinPrior

!=======================================================================
! Uniform[0:1]  ->  Cauchy[mean=x0,FWHM=2*gamma]

function CauchyPrior(r,x0,gamma)

            implicit none

            double precision r,x0,gamma,CauchyPrior
            real Pi
            parameter(Pi=3.141592654)

            CauchyPrior=x0+gamma*tan(Pi*(r-0.5))

end function CauchyPrior

!=======================================================================
! Uniform[0:1]  ->  Gaussian[mean=mu,variance=sigma**2]

function GaussianPrior(r,mu,sigma)

            implicit none

            double precision r,mu,sigma,GaussianPrior
            double precision SqrtTwo
            parameter(SqrtTwo=1.414213562d0)

            if (r.le.1.0d-16.or.(1.0d0-r).le.1.0d-16) then
            GaussianPrior=-1.0d32
            else
            GaussianPrior=mu+sigma*SqrtTwo*dierfc(2.d0*(1.d0-r))
            endif
      
end function GaussianPrior

!=======================================================================

! Uniform[0:1]  ->  LogNormal[mode=a,width parameter=sigma]

function LogNormalPrior(r,a,sigma)

            implicit none

            double precision r,a,sigma,LogNormalPrior
            double precision SqrtTwo,bracket
            parameter(SqrtTwo=1.414213562d0)

            bracket=sigma*sigma+sigma*SqrtTwo*dierfc(2.d0*r)
            LogNormalPrior=a*dexp(bracket)

end function LogNormalPrior

!=======================================================================

! Uniform[0:1]  ->  Beta[a,b]

function BetaPrior(r,a,b)

            implicit none

            double precision r, a, b, BetaPrior
            double precision beta_log
            integer ifault

            beta_log = lgamma(a) + lgamma(b) - lgamma(a+b)
            BetaPrior = xinbta(a, b, beta_log, r, ifault )
            if (ifault /= 0) print *, ifault

end function BetaPrior


!=======================================================================       
! Inverse of complimentary error function in double precision

function dierfc(y)

      implicit none
            double precision y,dierfc
            double precision qa,qb,qc,qd,q0,q1,q2,q3,q4,pa,pb,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18
            double precision p19,p20,p21,p22,x,z,w,u,s,t
            double precision infinity
            parameter (infinity=5.0d0)
            parameter (qa=9.16461398268964d-01, &
            qb=2.31729200323405d-01, &
            qc=4.88826640273108d-01, &
            qd=1.24610454613712d-01, &
            q0=4.99999303439796d-01, &
            q1=1.16065025341614d-01, &
            q2=1.50689047360223d-01, &
            q3=2.69999308670029d-01, &
            q4=-7.28846765585675d-02)
            parameter (pa=3.97886080735226000d0, &
            pb=1.20782237635245222d-01, &
            p0=2.44044510593190935d-01, &
            p1=4.34397492331430115d-01, &
            p2=6.86265948274097816d-01, &
            p3=9.56464974744799006d-01, &
            p4=1.16374581931560831d0, &
            p5=1.21448730779995237d0, &
            p6=1.05375024970847138d0, &
            p7=7.13657635868730364d-01, &
            p8=3.16847638520135944d-01, &
            p9=1.47297938331485121d-02, &
            p10=-1.05872177941595488d-01, &
            p11=-7.43424357241784861d-02)
            parameter (p12=2.20995927012179067d-03, &
            p13=3.46494207789099922d-02, &
            p14=1.42961988697898018d-02, &
            p15=-1.18598117047771104d-02, &
            p16=-1.12749169332504870d-02, &
            p17=3.39721910367775861d-03, &
            p18=6.85649426074558612d-03, &
            p19=-7.71708358954120939d-04, &
            p20=-3.51287146129100025d-03, &
            p21=1.05739299623423047d-04, &
            p22=1.12648096188977922d-03)
            if (y==0.0) then
            dierfc=infinity
            return
            endif  
            z=y
            if (y .gt. 1) z=2-y
            w=qa-log(z)
            u=sqrt(w)
            s=(qc+log(u))/w
            t=1/(u+qb)
            x=u*(1-s*(0.5d0+s*qd))-((((q4*t+q3)*t+q2)*t+q1)*t+q0)*t
            t=pa/(pa+x)
            u=t-0.5d0
            s=(((((((((p22*u+p21)*u+p20)*u+p19)*u+p18)*u+p17)*u+p16)*u+p15)*u+p14)*u+p13)*u+p12
            s=((((((((((((s*u+p11)*u+p10)*u+p9)*u+p8)*u+p7)*u+p6)*u+p5)*u+p4)*u+p3)*u+p2) &
            *u+p1)*u+p0)*t-z*exp(x*x-pb)
            x=x+s*(1+x*s)
            if (y .gt. 1) x=-x
            dierfc=x
      
end function dierfc



function betain ( x, p, q, beta, ifault )
    ! BETAIN computes the incomplete Beta function ratio.
    !  Author: Original FORTRAN77 version by KL Majumder, GP Bhattacharjee.
    !          FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input, real (kind=8) X, the argument, between 0 and 1.
    !    Input, real (kind=8) P, Q, the parameters, which must be positive.
    !    Input, real (kind=8) BETA, the logarithm of the complete beta function.
    !    Output, integer ( kind = 4 ) IFAULT, error flag. 0, no error.
    !    Output, real (kind=8) BETAIN, the value of the incomplete Beta function ratio.
  implicit none

  real (kind=8), parameter :: acu = 0.1D-14
  real (kind=8) ai, beta, betain, cx
  integer ( kind = 4 ) ifault, ns
  logical indx
  real (kind=8) p, pp, psq, q, qq, rx, temp, term, x, xx

  betain = x
  ifault = 0

  !  Check the input arguments.
  if ( p <= 0.0d0 ) then
    ifault = 1
    return
  end if

  if ( q <= 0.0d0 ) then
    ifault = 1
    return
  end if

  if ( x < 0.0d0 .or. 1.0d0 < x ) then
    ifault = 2
    return
  end if

  !  Special cases.
  if ( x == 0.0d0 .or. x == 1.0d0 ) then
    return
  end if

  !  Change tail if necessary and determine S.
  psq = p + q
  cx = 1.0d0 - x

  if ( p < psq * x ) then
    xx = cx
    cx = x
    pp = q
    qq = p
    indx = .true.
  else
    xx = x
    pp = p
    qq = q
    indx = .false.
  end if

  term = 1.0d0
  ai = 1.0d0
  betain = 1.0d0
  ns = int ( qq + cx * psq )

  !  Use Soper's reduction formula.
  rx = xx / cx
  temp = qq - ai
  if ( ns == 0 ) then
    rx = xx
  end if

  do
    term = term * temp * rx / ( pp + ai )
    betain = betain + term
    temp = abs ( term )

    if ( temp <= acu .and. temp <= acu * betain ) then
      betain = betain * exp ( pp * log ( xx ) + ( qq - 1.0d0 ) * log ( cx ) - beta ) / pp

      if ( indx ) then
        betain = 1.0d0 - betain
      end if
      exit
    end if

    ai = ai + 1.0d0
    ns = ns - 1

    if ( 0 <= ns ) then
      temp = qq - ai
      if ( ns == 0 ) then
        rx = xx
      end if
    else
      temp = psq
      psq = psq + 1.0d0
    end if

  end do

  return
end function betain

function xinbta ( p, q, beta, alpha, ifault )
    !  XINBTA computes the inverse of the incomplete Beta function.
    !
    !  Author: Original FORTRAN77 version by GW Cran, KJ Martin, GE Thomas.
    !           FORTRAN90 version by John Burkardt.
    !
    !  Parameters:
    !
    !    Input, real (kind=8) P, Q, the parameters of the incomplete Beta function.
    !
    !    Input, real (kind=8) BETA, the logarithm of the value of the complete Beta function.
    !
    !    Input, real (kind=8) ALPHA, the value of the incomplete Beta function.  
    !                         0 <= ALPHA <= 1.
    !
    !    Output, integer ( kind = 4 ) IFAULT, error flag. 0, no error occurred.
    !
    !    Output, real (kind=8) XINBTA, the argument of the incomplete Beta function 
    !                                  which produces the value ALPHA.
    !
    !  Local Parameters:
    !
    !    Local, real (kind=8) SAE, requests an accuracy of about 10^SAE.
  implicit none

  real (kind=8) a, acu, adj, alpha, beta, fpu, g, h
  integer ( kind = 4 ) iex
  integer ( kind = 4 ) ifault
  logical indx
  real (kind=8) p, pp, prev, q, qq, r, s
  real (kind=8), parameter :: sae = -30.0d0
  real (kind=8) sq, t, tx, w, value, xin, xinbta, y, yprev

  fpu = 10.0d0 ** sae

  ifault = 0
  value = alpha

  !  Test for admissibility of parameters.
  if ( p <= 0.0d0 ) then
    ifault = 1
    stop 1
  end if

  if ( q <= 0.0d0 ) then
    ifault = 1
    stop 1
  end if

  if ( alpha < 0.0d0 .or. 1.0d0 < alpha ) then
    ifault = 2
    stop 1
  end if

  !  Return immediately if the answer is easy to determine.  
  if ( alpha == 0.0d0 ) then
    value = 0.0d0
    xinbta = value
    return
  end if

  if ( alpha == 1.0d0 ) then
    value = 1.0d0
    xinbta = value
    return
  end if

  !  Change tail if necessary.
  if ( 0.5d0 < alpha ) then
    a = 1.0d0 - alpha
    pp = q
    qq = p
    indx = .true.
  else
    a = alpha
    pp = p
    qq = q
    indx = .false.
  end if

  !  Calculate the initial approximation.
  r = sqrt ( - log ( a * a ) )

  y = r - ( 2.30753d0 + 0.27061d0 * r ) / ( 1.0d0 + ( 0.99229d0 + 0.04481d0 * r ) * r )

  if ( 1.0d0 < pp .and. 1.0d0 < qq ) then
    r = ( y * y - 3.0d0 ) / 6.0d0
    s = 1.0d0 / ( pp + pp - 1.0d0 )
    t = 1.0d0 / ( qq + qq - 1.0d0 )
    h = 2.0d0 / ( s + t )
    w = y * sqrt ( h + r ) / h - ( t - s ) * ( r + 5.0d0 / 6.0d0 - 2.0d0 / ( 3.0d0 * h ) )
    value = pp / ( pp + qq * exp ( w + w ) )
  else
    r = qq + qq
    t = 1.0d0 / ( 9.0d0 * qq )
    t = r * ( 1.0d0 - t + y * sqrt ( t ) ) ** 3

    if ( t <= 0.0d0 ) then
      value = 1.0d0 - exp ( ( log ( ( 1.0d0 - a ) * qq ) + beta ) / qq )
    else
      t = ( 4.0d0 * pp + r - 2.0d0 ) / t
      if ( t <= 1.0d0 ) then
        value = exp ( ( log ( a * pp ) + beta ) / pp )
      else
        value = 1.0d0 - 2.0d0 / ( t + 1.0d0 )
      end if
    end if
  end if

  !  Solve for X by a modified Newton-Raphson method, using the function BETAIN.
  r = 1.0d0 - pp
  t = 1.0d0 - qq
  yprev = 0.0d0
  sq = 1.0d0
  prev = 1.0d0

  if ( value < 0.0001d0 ) then
    value = 0.0001d0
  end if

  if ( 0.9999d0 < value ) then
    value = 0.9999d0
  end if

  iex = max ( - 5.0d0 / pp / pp - 1.0d0 / a ** 0.2d0 - 13.0d0, sae )

  acu = 10.0d0 ** iex

  !  Iteration loop.
  do
    y = betain ( value, pp, qq, beta, ifault )
    if ( ifault /= 0 ) then
      write ( *, '(a,i6)' ) '  BETAIN returned IFAULT = ', ifault
      stop 1
    end if

    xin = value
    y = ( y - a ) * exp ( beta + r * log ( xin ) + t * log ( 1.0d0 - xin ) )

    if ( y * yprev <= 0.0d0 ) then
      prev = max ( sq, fpu )
    end if

    g = 1.0d0

    do
    !  Choose damping factor.
      do
        adj = g * y
        sq = adj * adj
        if ( sq < prev ) then
          tx = value - adj
          if ( 0.0d0 <= tx .and. tx <= 1.0d0 ) then
            exit
          end if
        end if
        g = g / 3.0d0
      end do

    !  Check whether current estimate is acceptable.
    !  The change "VALUE = TX" was suggested by Ivan Ukhov.
      if ( prev <= acu .or. y * y <= acu ) then
        value = tx
        if ( indx ) then
          value= 1.0d0 - value
        end if
        xinbta = value
        return
      end if
      if ( tx /= 0.0d0 .and. tx /= 1.0d0 ) then
        exit
      end if
      g = g / 3.0d0
    end do

    if ( tx == value ) then
      exit
    end if

    value = tx
    yprev = y

  end do

  if ( indx ) then
    value = 1.0d0 - value
  end if

  xinbta = value

  return
end function xinbta

!===========================================================)============
end module priors

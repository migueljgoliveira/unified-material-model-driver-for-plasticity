************************************************************************
c     GOTOH BIQUADRATIC YIELD FUNCTION AND DERIVATIVES
c
c       doi:
c
      subroutine ummdp_gotoh ( s,se,dseds,d2seds2,nreq,pryld,ndyld )          
c------------------------------------------------------------- variables
      implicit none
c
      integer,intent(in) :: nreq,ndyld
      real*8 ,intent(in) :: s(3),pryld(ndyld)
c
      real*8,intent(out) :: se
      real*8,intent(out) :: dseds(3)
			real*8,intent(out) :: d2seds2(3,3)
c
      integer i,j,k,m,n
      real*8 phi,q
      real*8 a(9),v(4),t(4)
      real*8 c(4,4),dtds(4,3)
      real*8 d2tds2(4,3,3)
c-----------------------------------------------------------------------
c     a(i)      : coef.s of Gotoh's 4th order function
c     t(i)      : stress^2 vector
c                 ={ sx^2, sx*sy, sy^2, txy^2 }
c     c(i,j)    : matrix to calc. se
c              Gotoh's function se = ( {t}^T*[c]*{t} )^(1/4)
c
c     dtds(i,j) : diff. of t(i) with respect to s(j)
c     d2tds2(i,j,k)
c               : 2nd order diff. of t(i) w.r.t s(j) & s(k)
c-----------------------------------------------------------------------
c
c                                            ---- anisotropic parameters
      do i = 1,9
        a(i) = pryld(1+i)
      end do
c                                                      ---- coef. matrix
      c(1,1) = a(1)
      c(1,2) = a(2)*0.5d0
      c(1,3) = 0.0
      c(1,4) = a(6)*0.5d0
      c(2,2) = a(3)
      c(2,3) = a(4)*0.5d0
      c(2,4) = a(7)*0.5d0
      c(3,3) = a(5)
      c(3,4) = a(8)*0.5d0
      c(4,4) = a(9)
      do i = 2,4
        do j = 1,i-1
          c(i,j) = c(j,i)
        end do
      end do
c                                                    ---- t-vector (s^2)
      t(1) = s(1) * s(1)
      t(2) = s(1) * s(2)
      t(3) = s(2) * s(2)
      t(4) = s(3) * s(3)
c                                                 ---- equivalent stress
      call ummdp_utility_mv  ( v,c,t,4,4 )
      call ummdp_utility_vvs ( phi,t,v,4 )
c
      if ( phi <= 0.0d0 ) phi = 0.0d0
      se = sqrt(sqrt(phi))
c                                            ---- 1st order differential
      if ( nreq >= 1 ) then
        call ummdp_utility_clear2 ( dtds,4,3 )
        dtds(1,1) = s(1) * 2.0d0
        dtds(2,1) = s(2)
        dtds(2,2) = s(1)
        dtds(3,2) = s(2) * 2.0d0
        dtds(4,3) = s(3) * 2.0d0
        call ummdp_utility_clear1 ( v,4 )
        do i = 1,3
          do j = 1,4
            do k = 1,4
              v(i) = v(i) + 2.0d0*t(j)*c(j,k)*dtds(k,i)
            end do
          end do
        end do
        q = 0.25d0 * phi**(-0.75d0)
        do i = 1,3
          dseds(i) = q * v(i)
        end do
      end if
c                                            ---- 2nd order differential
      if ( nreq >= 2 ) then
        call ummdp_utility_clear3 ( d2tds2,4,3,3 )
        d2tds2(1,1,1) = 2.0d0
        d2tds2(2,1,2) = 1.0d0
        d2tds2(2,2,1) = 1.0d0
        d2tds2(3,2,2) = 2.0d0
        d2tds2(4,3,3) = 2.0d0
        call ummdp_utility_clear2 ( d2seds2,3,3 )
        do i = 1,3
          do j = 1,3
            do m = 1,4
              do n = 1,4
                d2seds2(i,j) = d2seds2(i,j)+
     1                       2.0d0*c(m          ,n       )*
     2                        ( dtds(m,i)*dtds(  n  ,j)+
     3                           t(  m)  *d2tds2(n,i,j)  )
              end do
            end do
          end do
        end do
        do i = 1,3
          do j = 1,3
            d2seds2(i,j) = q*(d2seds2(  i,   j)
     1                      -0.75d0*v(i)*v(j)/phi)
          end do
        end do
      end if
c
      return
      end subroutine ummdp_gotoh
c
c
c
c
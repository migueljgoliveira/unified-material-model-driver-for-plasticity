************************************************************************
c     VON MISES YIELD FUNCTION AND DERIVATIVES
c
      subroutine ummdp_mises ( s,se,dseds,d2seds2,nreq )
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: nreq
      real*8 ,intent(in) :: s(6)
c
      real*8,intent(out) :: se
      real*8,intent(out) :: dseds(6)
			real*8,intent(out) :: d2seds2(6,6)
c
      integer i,j
			real*8 phi
      real*8 v(6)
			real*8 c(6,6)
c-----------------------------------------------------------------------
c	
c                                               ---- coefficients matrix
      call ummdp_utility_clear2 ( c,6,6 )
      do i = 1,3
        do j = 1,3
          c(i,j) = -0.5d0
        end do
        c(i,i) = 1.0d0
      end do
      do i = 4,6
        c(i,i) = 3.0d0
      end do
c
      call ummdp_utility_mv  ( v,c,s,6,6 )
      call ummdp_utility_vvs ( phi,s,v,6 )
c                                                 ---- equivalent stress
      se = sqrt(phi)
c                                              ---- 1st order derivative
      if ( nreq >= 1 ) then
        do i = 1,6
          dseds(i) = v(i) / se
        end do
      end if
c                                              ---- 2nd order derivative
      if ( nreq >= 2 ) then
        do i = 1,6
          do j = 1,6
            d2seds2(i,j) = (-v(i)*v(j)/phi+c(i,j)) / se
          end do
        end do
      end if
c
      return
      end subroutine ummdp_mises
c
c
c

************************************************************************
c
c     HU2005 YIELD FUNCTION
c
c       doi: https://doi.org/10.1016/j.ijplas.2004.11.004
c
      subroutine ummdp_yield_hu2005 ( s,se,dseds,d2seds2,nreq,pryld,
     1                                ndyld )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,parameter :: maxa = 100
c
      integer,intent(in) :: nreq,ndyld
      real*8 ,intent(in) :: s(6),pryld(ndyld)
c
      real*8,intent(out) :: se
      real*8,intent(out) :: dseds(6)
			real*8,intent(out) :: d2seds2(6,6)
c 
      integer nd0,n,nterms,it,jy,jx,ndmax
      real*8 a(maxa)
      real*8 ipow(maxa,3)
c-----------------------------------------------------------------------
c
      nd0 = 2
c
      n = 0
      do it = 0,nd0
        n = n + (nd0-it)*2 + 1
      end do
      nterms = n
      if ( maxa < nterms ) then
        write (6,*) 'increase maxa :',maxa,nterms
        call ummdp_exit ( 9000 )
      end if
c
      n = 0
      ipow = 0.0d0
      do it = 0,nd0
        ndmax = nd0*2 - it*2
        do jy = 0,ndmax
          jx = ndmax - jy
          n = n + 1
          ipow(n,1) = jx
          ipow(n,2) = jy
          ipow(n,3) = it
        end do
      end do
c
      a    =  0.0d0
      a(1) =  pryld(1+1)    ! X1 
      a(2) =  pryld(1+2)    ! X2
      a(3) =  pryld(1+3)    ! X3
      a(4) =  pryld(1+4)    ! X4
      a(5) =  pryld(1+5)    ! X5
      a(6) =  pryld(1+7)    ! C1 <-
      a(7) = -pryld(1+9)    ! C3 <-
      a(8) =  pryld(1+8)    ! C2 <-
      a(9) =  pryld(1+6)    ! X7
c
      call ummdp_yield_hy_polytype ( s,se,dseds,d2seds2,nreq,nd0,a,
     1                               ipow,maxa,nterms )
c
      return
      end subroutine ummdp_yield_hu2005
c
c
c
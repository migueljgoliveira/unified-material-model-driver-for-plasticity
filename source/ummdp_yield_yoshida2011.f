************************************************************************
c     YOSHIDA2011 YIELD FUNCTION AND DERIVATIVES
c
c       doi: https://doi.org/10.1016/j.ijplas.2013.01.010
c
      subroutine ummdp_yoshida2011 ( s,se,dseds,d2seds2,nreq,
     1                               pryld,ndyld )
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
      integer nd0,n,it,nterms,ndmax,jy,jx
      integer ipow(maxa,3)
      real*8 a(maxa)
c-----------------------------------------------------------------------
c
      nd0 = 3
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
      ipow = 0
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
      a = 0.0d0
      a(1) =  1.0d0 *          pryld(1+1)           !       c1
      a(2) = -3.0d0 *          pryld(1+2)           !    -3*c2
      a(3) =  6.0d0 *          pryld(1+3)           !     6*c3
      a(4) = -7.0d0 *          pryld(1+4)           !    -7*c4
      a(5) =  6.0d0 *          pryld(1+5)           !     6*c5
      a(6) =- 3.0d0 *          pryld(1+6)           !    -3*c6
      a(7) =  1.0d0 *          pryld(1+7)           !       c7
      a(8) =  1.0d0 * 9.0d0  * pryld(1+8)           !  9   *c8
      a(9) = -2.0d0 * 9.0d0  * pryld(1+9)           !  9*-2*c9
      a(10) = 3.0d0 * 9.0d0  * pryld(1+10)          !  9* 3*c10
      a(11) =-2.0d0 * 9.0d0  * pryld(1+11)          !  9*-2*c11
      a(12) = 1.0d0 * 9.0d0  * pryld(1+12)          !  9   *c12
      a(13) = 1.0d0 * 27.0d0 * pryld(1+13)          ! 27   *c13
      a(14) =-1.0d0 * 27.0d0 * pryld(1+14)          ! 27*-1*c14
      a(15) = 1.0d0 * 27.0d0 * pryld(1+15)          ! 27   *c15
      a(16) = 1.0d0 * 27.0d0 * pryld(1+16)          ! 27   *c16
c
      call ummdp_hy_polytype ( s,se,dseds,d2seds2,nreq,nd0,
     1                          a,ipow,maxa,nterms )
c
      return
      end subroutine ummdp_yoshida2011
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     HU2005 & YOSHIDA2011 STYLE POLYNOMIAL TYPE YIELD FUNCTION
c
      subroutine ummdp_hy_polytype ( s,se,dseds,d2seds2,nreq,nd0,
     1                                a,ipow,maxa,nterms )
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: nreq,nd0,maxa,nterms
      integer,intent(in) :: ipow(maxa,3)
      real*8 ,intent(in) :: s(6),a(maxa)
c
      real*8,intent(out) :: se
      real*8,intent(out) :: dseds(6)
      real*8,intent(out) :: d2seds2(6,6)
c    
      integer i,j,k,n,id,idmax,jdmax,jd,nd
      integer ii(3)
      real*8 dinv,fai,q,dd,ff,fff,ddi,ddj
      real*8 sterm(3),v(6)
c-----------------------------------------------------------------------
c     nd        : order of polynominal nd=2*nd0
c     a(n)      : constants of function
c     ipow(n,i) : power of terms
c-----------------------------------------------------------------------
c
      nd = nd0 * 2
      dinv = 1.0d0 / float(nd)
c
      sterm(1) = s(1) - s(3)                  ! sx-sz
      sterm(2) = s(2) - s(3)                  ! sy-sz
      sterm(3) = s(4)**2 + s(5)**2 + s(6)**2  ! txy^2+tyz^2+tzx^2
c
      fai = 0.0d0
      do n = 1,nterms
        q = a(n)
        do k = 1,3
          if ( ipow(n,k) > 0 ) then
            q = q * sterm(k)**ipow(n,k)
          end if
        end do
        fai = fai + q
      end do
      se = fai ** dinv
      if ( nreq == 0 ) return
c
      v = 0.0d0
      do i = 1,6
        idmax = 1
        if ( i == 3 ) idmax = 2
        do id = 1,idmax
          do n = 1,nterms
            do k = 1,3
              ii(k) = ipow(n,k)
            end do
            select case ( i )
            case ( 1,2 )
              dd = float(ii(i))
              ii(i) = ii(i) - 1
            case ( 3 )
              dd = -1.0d0 * float(ii(id))
              ii(id) = ii(id) - 1
            case default
              dd = 2.0d0 * s(i) * float(ii(3))
              ii(3) = ii(3) - 1
            end select
            q = dd * a(n)
            do k = 1,3
              if ( ii(k) > 0 ) then
                q = q * sterm(k)**ii(k)
              else if ( ii(k) < 0 ) then
                q = 0.0
              end if
            end do
            v(i) = v(i) + q
          end do
        end do
      end do
      ff = dinv * fai**(dinv-1.0d0)
      do i = 1,6
        dseds(i) = ff * v(i)
      end do
      if ( nreq == 1 ) return
c
      fff = dinv * (dinv-1.0d0) * fai**(dinv-2.0d0)
      do i = 1,6
        do j = 1,6
          d2seds2(i,j) = fff * v(i) * v(j)
        end do
      end do
      do i = 1,6
        idmax = 1
        if ( i == 3 ) idmax = 2
        do id = 1,idmax
          do j = 1,6
            jdmax = 1
            if (  j == 3 ) jdmax = 2
            if ( ( j > 3 ) .and. ( i == j ) ) jdmax = 2
            do jd = 1,jdmax
              do n = 1,nterms
                do k = 1,3
                  ii(k) = ipow(n,k)
                end do
                select case ( i )
                case ( 1,2 )
                  ddi = float(ii(i))
                  ii(i) = ii(i) - 1
                case ( 3 )
                  ddi = -1.0d0 * float(ii(id))
                  ii(id) = ii(id) - 1
                case default
                  ddi = 2.0d0 * s(i) * float(ii(3))
                  ii(3) = ii(3) - 1
                end select 
                select case ( j )
                case ( 1,2 )
                  ddj = float(ii(j))
                  ii(j) = ii(j) - 1
                case ( 3 )
                  ddj = -1.0d0 * float(ii(jd))
                  ii(jd) = ii(jd) - 1
                case default
                  if ( jd == 1 ) then
                    ddj = 2.0d0 * s(j) * float(ii(3))
                    ii(3) = ii(3) - 1
                  else
                    ddi = 2.0d0 * float(ipow(n,3))
                    ddj = 1.0d0
                  end if
                end select
                q = a(n) * ddi * ddj
                do k = 1,3
                  if ( ii(k) > 0 ) then
                    q = q * sterm(k)**ii(k)
                  else if ( ii(k) < 0 ) then
                    q = 0.0d0
                  end if
                end do
                d2seds2(i,j) = d2seds2(i,j) + ff*q
              end do
            end do
          end do
        end do
      end do
c
      return
      end subroutine ummdp_hy_polytype
c
c
c
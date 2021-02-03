c----------------------------------------------------------(yoshida2011)
c     F.Yoshida (2011,2013) yield function and its dfferentials
c
c     NUMISHEET 2011 Proceedings, AIP Conf. Proc.1383 (2011), pp.807-814
c
c     "A user-friendly 3D yield function to describe anisotropy of 
c       steel sheets ",IJP,v.45(2013), pp.1119-139. )
c
      subroutine jancae_yoshida_2011 ( s,se,dseds,d2seds2,nreq,
     &                                 pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (maxa=100)
      dimension s(6),dseds(6),d2seds2(6,6),pryld(ndyld)
      dimension a(maxa),ipow(maxa,3)
c
c
      nd0 = 3
c
      n = 0
      do it = 0,nd0
        n = n + (nd0-it)*2 + 1
      end do
      nterms = n
      if ( maxa .lt. nterms ) then
        write (6,*) 'increase maxa :',maxa,nterms
        call jancae_exit ( 9000 )
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
      a = 0.0
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
      call jancae_hy_polytype ( s,se,dseds,d2seds2,nreq,
     &                          nd0,a,ipow,maxa,nterms )
c
      return
      end
c
c
c
c----------------------------------------------------------(yoshida2011)
c     Weilong Hu & F.Yoshida style polynominal type yield function
c
      subroutine jancae_hy_polytype ( s,se,dseds,d2seds2,nreq,
     &                                nd0,a,ipow,maxa,nterms )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension s(6),dseds(6),d2seds2(6,6)
      dimension a(maxa),ipow(maxa,3)
      dimension sterm(3),ii(3),v(6)
c
c     nd        : order of polynominal nd=2*nd0
c     a(n)      : constants of function
c     ipow(n,i) : power of terms
c
      nd = nd0 * 2
      dinv = 1.0d0 / float(nd)
c
      sterm(1) = s(1) - s(3)                  ! sx-sz
      sterm(2) = s(2) - s(3)                  ! sy-sz
      sterm(3) = s(4)**2 + s(5)**2 + s(6)**2  ! txy^2+tyz^2+tzx^2
c
      fai = 0.0
      do n = 1,nterms
        q = a(n)
        do k = 1,3
          if ( ipow(n,k) .gt. 0 ) then
            q = q * sterm(k)**ipow(n,k)
          end if
        end do
        fai = fai + q
      end do
      se = fai ** dinv
      if ( nreq .eq. 0 ) return
c
      v = 0.0
      do i = 1,6
        idmax = 1
        if ( i .eq. 3 ) idmax = 2
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
              if ( ii(k) .gt. 0 ) then
                q = q * sterm(k)**ii(k)
              else if ( ii(k) .lt. 0 ) then
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
      if ( nreq .eq. 1 ) return
c
      fff = dinv * (dinv-1.0d0) * fai**(dinv-2.0d0)
      do i = 1,6
        do j = 1,6
          d2seds2(i,j) = fff * v(i) * v(j)
        end do
      end do
      do i = 1,6
        idmax = 1
        if ( i .eq. 3 ) idmax = 2
        do id = 1,idmax
          do j = 1,6
            jdmax = 1
            if (  j .eq. 3 ) jdmax = 2
            if ( ( j .gt. 3 ) .and. ( i .eq. j ) ) jdmax = 2
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
                  if ( jd .eq. 1 ) then
                    ddj = 2.0d0 * s(j) * float(ii(3))
                    ii(3) = ii(3) - 1
                  else
                    ddi = 2.0d0 * float(ipow(n,3))
                    ddj = 1.0d0
                  end if
                end select
                q = a(n) * ddi * ddj
                do k = 1,3
                  if ( ii(k) .gt. 0 ) then
                    q = q * sterm(k)**ii(k)
                  else if ( ii(k) .lt. 0 ) then
                    q = 0.0
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
      end
c
c
c
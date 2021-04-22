c----------------------------------------------------------------(mises)
c     von Mises yield function and its dfferentials
c     ( 1913 )
c
      subroutine jancae_mises ( s,se,dseds,d2seds2,nreq )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension s(6),dseds(6),d2seds2(6,6)
      dimension c(6,6),v(6)
c                                                      ---- coef. matrix
      call jancae_clear2 ( c,6,6 )
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
      call jancae_mv  ( v,c,s,6,6 )
      call jancae_vvs ( phi,s,v,6 )
c                                                 ---- equivalent stress
      se = sqrt(phi)
c                                            ---- 1st order differential
      if ( nreq >= 1 ) then
        do i = 1,6
          dseds(i) = v(i) / se
        end do
      end if
c                                            ---- 2nd order differential
      if ( nreq >= 2 ) then
        do i = 1,6
          do j = 1,6
            d2seds2(i,j) = (-v(i)*v(j)/phi+c(i,j)) / se
          end do
        end do
      end if
c
      return
      end
c
c
c

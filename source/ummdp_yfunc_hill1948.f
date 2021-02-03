c-------------------------------------------------------------(hill1948)
c     Hill 1948 yield function and its dfferentials
c     (  Proc. Roy. Soc. A193(1948) p281-297 )
c
c     ( flow curve must be defined in uniaxial sx vs ex )
c
      subroutine jancae_hill1948 ( s,se,dseds,d2seds2,nreq,
     &                             pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension s(6),dseds(6),d2seds2(6,6),pryld(ndyld)
      dimension c(6,6),v(6)
c
c     c(i,j)    : matrix to calc. se
c     Hill's 1948 function se = ({s}^T*[c]*{s})^(1/2)
c
c                                            ---- anisotropic parameters
c                              pf means "F", pg means "G"....
c                         F=G=H=1 and L=M=N=3 means von Mises
      pf = pryld(1+1)
      pg = pryld(1+2)
      ph = pryld(1+3)
      pl = pryld(1+4)
      pm = pryld(1+5)
      pn = pryld(1+6)
c                                                      ---- coef. matrix
      call jancae_clear2 ( c,6,6 )
      c(1,1) = pg + ph
      c(1,2) = -ph
      c(1,3) = -pg
      c(2,1) = -ph
      c(2,2) = pf + ph
      c(2,3) = -pf
      c(3,1) = -pg
      c(3,2) = -pf
      c(3,3) = pf + pg
      c(4,4) = 2.0d0 * pn
      c(5,5) = 2.0d0 * pm
      c(6,6) = 2.0d0 * pl
      do i = 1,6
        do j = 1,6
          c(i,j) = c(i,j) / (pg+ph)
        end do
      end do
c
      call jancae_mv  ( v,c,s,6,6 )
      call jancae_vvs ( phi,s,v,6 )
c                                                 ---- equivalent stress
      if ( phi .le. 0.0 ) phi = 0.0
      se = sqrt(phi)
c                                            ---- 1st order differential
      if ( nreq .ge. 1 ) then
        do i = 1,6
          dseds(i) = v(i) / se
        end do
      end if
c                                            ---- 2nd order differential
      if ( nreq .ge. 2 ) then
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

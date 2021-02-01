c----------------------------------------------------------------(gotoh)
c     Gotoh biquadratic yield function and its dfferentials
c     (J.JSTP vol.19 no.208 p377-385 (1978-5) )
c
      subroutine jancae_gotoh ( s,se,dseds,d2seds2,nreq,
     &                          pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension s(3),dseds(3),d2seds2(3,3),pryld(ndyld)
      dimension a(9),c(4,4),v(4),
     &          t(4),dtds(4,3),d2tds2(4,3,3)
c
c     a(i)      : coef.s of Gotoh's 4th order function
c     t(i)      : stress^2 vector 
c                 ={ sx^2, sx*sy, sy^2, txy^2 }
c     c(i,j)    : matrix to calc. se
c              Gotoh's function se = ( {t}^T*[c]*{t} )^(1/4)
c
c     dtds(i,j) : diff. of t(i) with respect to s(j)
c     d2tds2(i,j,k)
c               : 2nd order diff. of t(i) w.r.t s(j) & s(k)
c
c                                            ---- anisotropic parameters
      do i=1,9
        a(i)=pryld(1+i)
      enddo
c                                                      ---- coef. matrix
      c(1,1)=a(1)
      c(1,2)=a(2)*0.5d0
      c(1,3)=0.0
      c(1,4)=a(6)*0.5d0
      c(2,2)=a(3)
      c(2,3)=a(4)*0.5d0
      c(2,4)=a(7)*0.5d0
      c(3,3)=a(5)
      c(3,4)=a(8)*0.5d0
      c(4,4)=a(9)
      do i=2,4
        do j=1,i-1
          c(i,j)=c(j,i)
        enddo
      enddo
c                                                    ---- t-vector (s^2)
      t(1)=s(1)*s(1)
      t(2)=s(1)*s(2)
      t(3)=s(2)*s(2)
      t(4)=s(3)*s(3)
c                                                 ---- equivalent stress
      call jancae_mv  ( v,c,t,4,4 )
      call jancae_vvs ( phi,t,v,4 )
c
      if ( phi.le.0.0 ) phi=0.0
      se=sqrt(sqrt(phi))
c                                            ---- 1st order differential
      if ( nreq.ge.1 ) then
        call jancae_clear2 ( dtds,4,3 )
        dtds(1,1)=s(1)*2.0d0
        dtds(2,1)=s(2)
        dtds(2,2)=s(1)
        dtds(3,2)=s(2)*2.0d0
        dtds(4,3)=s(3)*2.0d0
        call jancae_clear1 ( v,4 )
        do i=1,3
          do j=1,4
            do k=1,4
              v(i)=v(i)+2.0d0*t(j)*c(j,k)*dtds(k,i)
            enddo
          enddo
        enddo
        q=0.25d0*phi**(-0.75d0)
        do i=1,3
          dseds(i)=q*v(i)
        enddo
      endif
c                                            ---- 2nd order differential
      if ( nreq.ge.2 ) then
        call jancae_clear3 ( d2tds2,4,3,3 )
        d2tds2(1,1,1)=2.0d0
        d2tds2(2,1,2)=1.0d0
        d2tds2(2,2,1)=1.0d0
        d2tds2(3,2,2)=2.0d0
        d2tds2(4,3,3)=2.0d0
        call jancae_clear2 ( d2seds2,3,3 )
        do i=1,3
          do j=1,3
            do m=1,4
              do n=1,4
                d2seds2(i,j)=d2seds2(i,j)+
     &                       2.0d0*c(m          ,n       )*
     &                        ( dtds(m,i)*dtds(  n  ,j)+
     &                           t(  m)  *d2tds2(n,i,j)  )
              enddo
            enddo
          enddo
        enddo
        do i=1,3
          do j=1,3
            d2seds2(i,j)=q*(d2seds2(  i,   j)
     &                      -0.75d0*v(i)*v(j)/phi)
          enddo
        enddo
      endif
c
      return
      end
c
c
c
c
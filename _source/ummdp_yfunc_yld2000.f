c--------------------------------------------------------------(yld2000)
c     Barlat YLD2000 yield function and its dfferentials
c     ( IJP v.19(203) p1297-1319. )
c

      subroutine jancae_yld2000 ( s,se,dseds,d2seds2,nreq,
     &                            pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension s(3),dseds(3),d2seds2(3,3),pryld(ndyld)
      dimension a(8),am(2,3,3),x(2,2),y(2,3),phi(2),
     &          dsedphi(2),dphidx(2,2),
     &          dxdy(2,2,3),dyds(2,3,3),
     &          d2sedphi2(2,2),d2phidx2(2,2,2),
     &          d2xdy2(2,2,3,3)
c
c
c     variables  : symbols in Barlat's paper
c
c     s(i)       : Sigma (i=1~2)
c     x(1,i)     : X'i   (i=1~2)
c     x(2,i)     : X"i   (i=1~2)
c     y(1,i)     : X'xx,X'yy,X'xy (i=1~3)
c     y(2,i)     : X"xx,X"yy,X"xy (i=1~3)
c     phi(1)     : phi'
c     phi(2)     : phi"
c     se         : equivalent stress =(phi'+ph")^(1/M)
c     am(1,i,j)  : liner transf. matrix for sigma to X'xx
c     am(2,i,j)  : liner transf. matrix for sigma to X"xx
c     a(i)       : anisotropic parameter a1~a8
c
c       1st index means number of dash
c       2nd index means suffix
c
c                                            ---- anisotropic parameters
      do i = 1,8
        a(i) = pryld(i+1)
      enddo
      em = pryld(9+1)
c                                         ---- set linear transf. matrix
      call jancae_yld2000_2d_am ( a,am )
c                                                 ---- equivalent stress
      call jancae_yld2000_2d_xyphi ( s,em,am,x,y,phi )
      q = phi(1) + phi(2)
      if ( q .le. 0.0 ) q = 0.0
      se = (0.5d0*q) ** (1.0d0/em)
c                                            ---- 1st order differential
      if ( nreq .ge. 1 ) then
        call jancae_yld2000_2d_ds1 ( em,am,x,y,phi,
     &                               dsedphi,dphidx,
     &                               dxdy,dyds,se )
        call jancae_clear1 ( dseds,3 )
        do nd = 1,2
          do m = 1,2
            do k = 1,3
              do i = 1,3
                dseds(i) = dseds(i) + dsedphi(nd)*dphidx(nd,m)*
     &                                dxdy(nd,m,k)*dyds(nd,k,i)
              enddo
            enddo
          enddo
        enddo
      endif
c                                            ---- 2nd order differential
      if ( nreq .ge. 2 ) then
        call jancae_yld2000_2d_ds2 ( phi,x,y,em,
     &                               d2sedphi2,d2phidx2,
     &                               d2xdy2,se )
        call jancae_clear2 ( d2seds2,3,3 )
        do i = 1,3
        do j = 1,3
          do nd1 = 1,2
          do nd2 = 1,2
            do k = 1,2
            do l = 1,2
              do m = 1,3
              do n = 1,3
                d2seds2(i,j) = d2seds2(i,j) + d2sedphi2(nd1,nd2)*
     &                          dphidx(nd1,k)*
     &                          dxdy(nd1,k,m)*dyds(nd1,m,i)*
     &                          dphidx(nd2,l)*
     &                          dxdy(nd2,l,n)*dyds(nd2,n,j)
              enddo
              enddo
            enddo
            enddo
          enddo
          enddo
          do nd = 1,2
            do k = 1,2
            do l = 1,2
              do m = 1,3
              do n = 1,3
                d2seds2(i,j) = d2seds2(i,j) + dsedphi(nd)*
     &                          d2phidx2(nd,k,l)*
     &                          dxdy(nd,k,m)*dyds(nd,m,i)*
     &                          dxdy(nd,l,n)*dyds(nd,n,j)
              enddo
              enddo
            enddo
            enddo
          enddo
          do nd = 1,2
            do k = 1,2
              do m = 1,3
              do n = 1,3
                d2seds2(i,j) = d2seds2(i,j) + dsedphi(nd)*
     &                          dphidx(nd,k)*d2xdy2(nd,k,m,n)*
     &                          dyds(nd,m,i)*dyds(nd,n,j)
              enddo
              enddo
            enddo
          enddo
        enddo
        enddo
      endif
c
      return
      end
c
c
c
c--------------------------------------------------------------(yld2000)
c     set barlat-yld2k linear transformation matrix am
c
      subroutine jancae_yld2000_2d_am ( a,am )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(8),am(2,3,3)
c
c                                      ---- linear transformation matrix
      am(1,1,1) =  2.0d0*a(1)
      am(1,1,2) = -1.0d0*a(1)
      am(1,1,3) =  0.0
      am(1,2,1) = -1.0d0*a(2)
      am(1,2,2) =  2.0d0*a(2)
      am(1,2,3) =  0.0
      am(1,3,1) =  0.0
      am(1,3,2) =  0.0
      am(1,3,3) =  3.0d0*a(7)
c
      am(2,1,1) = -2.0d0*a(3) + 2.0d0*a(4) + 8.0d0*a(5) - 2.0d0*a(6)
      am(2,1,2) =        a(3) - 4.0d0*a(4) - 4.0d0*a(5) + 4.0d0*a(6)
      am(2,1,3) =  0.0
      am(2,2,1) =  4.0d0*a(3) - 4.0d0*a(4) - 4.0d0*a(5) +       a(6)
      am(2,2,2) = -2.0d0*a(3) + 8.0d0*a(4) + 2.0d0*a(5) - 2.0d0*a(6)
      am(2,2,3) =  0.0
      am(2,3,1) =  0.0
      am(2,3,2) =  0.0
      am(2,3,3) =  9.0d0*a(8)
c
      do i = 1,3
        do j = 1,3
          am(1,i,j) = am(1,i,j) / 3.0d0
          am(2,i,j) = am(2,i,j) / 9.0d0
        enddo
      enddo
c
      return
      end
c
c
c
c--------------------------------------------------------------(yld2000)
c     calc. barlat-yld2k function x,y,phi
c
      subroutine jancae_yld2000_2d_xyphi ( s,em,am,x,y,phi )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension s(3),am(2,3,3),x(2,2),y(2,3),phi(2)
      dimension p(2)
c
      p(1) =  1.0d0
      p(2) = -1.0d0
c                                                       ---- {y}=[am]{s}
      call jancae_clear2 ( y,2,3 )
      do nd = 1,2
        do i = 1,3
          do j = 1,3
            y(nd,i) = y(nd,i) + am(nd,i,j)*s(j)
          enddo
        enddo
      enddo
c                                        ---- {x}=principle value of {y}
      do nd = 1,2
        a = (y(nd,1)-y(nd,2))**2.d0 + 4.0d0*y(nd,3)**2.d0
        a = sqrt(a)
        do i = 1,2
          x(nd,i) = 0.5d0 * (y(nd,1)+y(nd,2)+p(i)*a)
        enddo
      enddo
c                                                 ---- phi(1) and phi(2)
      nd = 1
      phi(nd) = abs(x(nd,1)-x(nd,2))**em
      nd = 2
      phi(nd) = abs(2.0d0*x(nd,2)+x(nd,1))**em +
     &          abs(2.0d0*x(nd,1)+x(nd,2))**em
c
      return
      end
c
c
c
c--------------------------------------------------------------(yld2000)
c     set 1st order differential of parameters
c
      subroutine jancae_yld2000_2d_ds1 ( em,am,x,y,phi,
     &                                   dsedphi,dphidx,
     &                                   dxdy,dyds,se )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension am(2,3,3),x(2,2),y(2,3),phi(2),
     &          dsedphi(2),dphidx(2,2),
     &          dxdy(2,2,3),dyds(2,3,3)
c
      dimension p(2)
c
      eps = 1.0d-16
c
      p(1) =  1.0d0
      p(2) = -1.0d0
      emi = 1.0d0 / em
c                                                          ---- dse/dphi
      q = phi(1) + phi(2)
      if ( q .le. 0.0 ) q = 0.0
      do nd = 1,2
        dsedphi(nd) = (0.5d0**emi) * emi * q**(emi-1.0d0)
      enddo
c                                                           ---- dphi/dx
      nd = 1
      a0 = x(nd,1) - x(nd,2)
      b0 = abs(a0)
      sgn0 = 0
      if ( b0 .ge. eps*se ) sgn0 = a0 / b0
      dphidx(nd,1) =  em * b0**(em-1.0d0) * sgn0
      dphidx(nd,2) = -em * b0**(em-1.0d0) * sgn0
c
      nd = 2
      a1 = 2.0d0*x(nd,1) +       x(nd,2)
      a2 =       x(nd,1) + 2.0d0*x(nd,2)
      b1 = abs(a1)
      b2 = abs(a2)
      sgn1 = 0.0
      sgn2 = 0.0
      if ( b1 .ge. eps*se ) sgn1 = a1 / b1
      if ( b2 .ge. eps*se ) sgn2 = a2 / b2
      dphidx(nd,1) = em*(2.0d0*b1**(em-1.0d0)*sgn1 +
     &                         b2**(em-1.0d0)*sgn2 )
      dphidx(nd,2) = em*(      b1**(em-1.0d0)*sgn1 +
     &                   2.0d0*b2**(em-1.0d0)*sgn2 )
c
      do nd = 1,2
        a = (y(nd,1)-y(nd,2))*(y(nd,1)-y(nd,2)) + 4.0d0*y(nd,3)*y(nd,3)
        a = sqrt(a)
        if ( a .gt. eps*se ) then
          do j = 1,2
            dxdy(nd,j,1) = 0.5d0 * (1.0d0+p(j)*(y(nd,1)-y(nd,2))/a)
            dxdy(nd,j,2) = 0.5d0 * (1.0d0-p(j)*(y(nd,1)-y(nd,2))/a)
            dxdy(nd,j,3) = 2.0d0 *        p(j)* y(nd,3)         /a
          enddo
        else
          do j = 1,2
            dxdy(nd,j,1) = 0.5d0 * (1.0d0+0.0)
            dxdy(nd,j,2) = 0.5d0 * (1.0d0-0.0)
            dxdy(nd,j,3) = 2.0d0 *        0.0
          enddo
        endif
      enddo
c
      do nd = 1,2
        do i = 1,3
          do j = 1,3
            dyds(nd,i,j) = am(nd,i,j)
          enddo
        enddo
      enddo
c
      return
      end
c
c
c
c--------------------------------------------------------------(yld2000)
c     set 2nd order differential of parameters
c
      subroutine jancae_yld2000_2d_ds2 ( phi,x,y,em,
     &                                   d2sedphi2,d2phidx2,
     &                                   d2xdy2,se )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension phi(2),x(2,2),y(2,3),
     &          d2sedphi2(2,2),d2phidx2(2,2,2),
     &          d2xdy2(2,2,3,3)
      dimension p(2)
c
      eps = 1.0d-16
c
      p(1) =  1.0d0
      p(2) = -1.0d0
      emi = 1.0d0 / em
c                                                        ---- d2se/dphi2
      q = phi(1) + phi(2)
      if ( q .le. 0.0 ) q = 0.0
      do nd1 = 1,2
        do nd2 = 1,2
          a = 0.5d0**emi * emi * (emi-1.0d0) * q**(emi-2.0d0)
          d2sedphi2(nd1,nd2) = a
        enddo
      enddo
c                                                         ---- d2phi/dx2
      nd = 1
      do i = 1,2
        do j = 1,2
          a = (em-1.0d0) * em * (abs(x(nd,1)-x(nd,2)))**(em-2.0d0)
          if ( i .ne. j ) a = -a
          d2phidx2(nd,i,j) = a
        enddo
      enddo
      nd = 2
      do i = 1,2
        do j = 1,2
          if ( i .eq. j ) then
            if ( i .eq. 1 ) then
              a = (em-1.0d0) * em*
     &            (4.0d0*(abs(2.0d0*x(nd,1)+      x(nd,2)))**(em-2.0d0)+
     &                   (abs(      x(nd,1)+2.0d0*x(nd,2)))**(em-2.0d0))
            else
              a = (em-1.0d0)*em*
     &            (      (abs(2.0d0*x(nd,1)+      x(nd,2)))**(em-2.0d0)+
     &            4.0d0*(abs(      x(nd,1)+2.0d0*x(nd,2)))**(em-2.0d0))
            endif
          else
            a = (em-1.0d0) * em * 
     &          (2.0d0*(abs(2.0d0*x(nd,1)+      x(nd,2)))**(em-2.0d0)+
     &           2.0d0*(abs(      x(nd,1)+2.0d0*x(nd,2)))**(em-2.0d0) )
          endif
          d2phidx2(nd,i,j) = a
        enddo
      enddo
c                                                           ---- d2x/dy2
      do nd = 1,2
        a = (y(nd,1)-y(nd,2))*(y(nd,1)-y(nd,2)) +
     &      4.0d0*   y(nd,3) *         y(nd,3)
        if ( a .gt. eps*se ) then
          a = 1.0d0 / sqrt(a**3)
          do m = 1,2
            do i = 1,3
              do j = 1,3
                ij = i*10+j
                if ( ( ij .eq. 11 ) .or. ( ij .eq. 22 ) ) then
                  q = y(nd,3) * y(nd,3)
                else if ( ij .eq. 33 ) then
                  q = (y(nd,1)-y(nd,2)) * (y(nd,1)-y(nd,2))
                else if ( ( ij .eq. 12 ) .or. ( ij .eq. 21 ) ) then
                  q = -y(nd,3) * y(nd,3)
                else if ( ( ij .eq. 23 ) .or. ( ij .eq. 32 ) ) then
                  q = y(nd,3) * (y(nd,1)-y(nd,2))
                else
                  q = -y(nd,3) * (y(nd,1)-y(nd,2))
                endif
                d2xdy2(nd,m,i,j) = 2.0d0 * a * p(m) * q
              enddo
            enddo
          enddo
        else
          do m = 1,2
            do i = 1,3
              do j = 1,3
                d2xdy2(nd,m,i,j) = 0.0
              enddo
            enddo
          enddo
        endif
      enddo
c
      return
      end
c
c
c

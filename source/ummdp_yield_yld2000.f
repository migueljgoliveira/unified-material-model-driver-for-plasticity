************************************************************************
c     YLD2000-2D YIELD FUNCTION AND DERIVATIVES
c
c       doi: https://doi.org/10.1016/S0749-6419(02)00019-0
c
      subroutine ummdp_yld2000 ( s,se,dseds,d2seds2,nreq,pryld,ndyld )     
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: nreq,ndyld
      real*8 ,intent(in) :: s(3),pryld(ndyld)
c
      real*8,intent(out) :: se
      real*8,intent(out) :: dseds(3)
      real*8,intent(out) :: d2seds2(3,3)
c
      integer i,j,k,l,m,n,nd,nd1,nd2
      real*8 em,q
      real*8 a(8),phi(2),dsedphi(2),dphidx(2,2)
      real*8 x(2,2),y(2,3),d2sedphi2(2,2)
      real*8 am(2,3,3),dxdy(2,2,3),dyds(2,3,3),d2phidx2(2,2,2)
      real*8 d2xdy2(2,2,3,3)
c-----------------------------------------------------------------------         
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
      end do
      em = pryld(9+1)
c                                  ---- set linear transformation matrix
      call ummdp_yld2000_2d_am ( a,am )
c                                                 ---- equivalent stress
      call ummdp_yld2000_2d_xyphi ( s,em,am,x,y,phi )
      q = phi(1) + phi(2)
      if ( q <= 0.0 ) q = 0.0
      se = (0.5d0*q) ** (1.0d0/em)
c                                              ---- 1st order derivative
      if ( nreq >= 1 ) then
        call ummdp_yld2000_2d_ds1 ( em,am,x,y,phi,dsedphi,dphidx,dxdy,
     1                              dyds,se )
        call ummdp_utility_clear1 ( dseds,3 )
        do nd = 1,2
          do m = 1,2
            do k = 1,3
              do i = 1,3
                dseds(i) = dseds(i) + (dsedphi(nd)*dphidx(nd,m)
     1                                 * dxdy(nd,m,k)*dyds(nd,k,i))
              end do
            end do
          end do
        end do
      end if
c                                              ---- 2nd order derivative
      if ( nreq >= 2 ) then
        call ummdp_yld2000_2d_ds2 ( phi,x,y,em,d2sedphi2,d2phidx2,
     1                              d2xdy2,se )                   
        call ummdp_utility_clear2 ( d2seds2,3,3 )
        do i = 1,3
        do j = 1,3
          do nd1 = 1,2
          do nd2 = 1,2
            do k = 1,2
            do l = 1,2
              do m = 1,3
              do n = 1,3
                d2seds2(i,j) = d2seds2(i,j) 
     1                         + d2sedphi2(nd1,nd2)*dphidx(nd1,k)                     
     2                           * (dxdy(nd1,k,m)*dyds(nd1,m,i)
     3                              * dphidx(nd2,l)*dxdy(nd2,l,n)
     4                              * dyds(nd2,n,j))
              end do
              end do
            end do
            end do
          end do
          end do
          do nd = 1,2
            do k = 1,2
            do l = 1,2
              do m = 1,3
              do n = 1,3
                d2seds2(i,j) = d2seds2(i,j) 
     1                          + (dsedphi(nd)*d2phidx2(nd,k,l)
     2                             * dxdy(nd,k,m)*dyds(nd,m,i)
     3                             * dxdy(nd,l,n)*dyds(nd,n,j))
              end do
              end do
            end do
            end do
          end do
          do nd = 1,2
            do k = 1,2
              do m = 1,3
              do n = 1,3
                d2seds2(i,j) = d2seds2(i,j) 
     1                         + (dsedphi(nd)*dphidx(nd,k)
     2                            * d2xdy2(nd,k,m,n)*dyds(nd,m,i)
     3                            * dyds(nd,n,j))
              end do
              end do
            end do
          end do
        end do
        end do
      end if
c
      return
      end subroutine ummdp_yld2000
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     SET LINEAR TRANSFORMATION MATRIX
c
      subroutine ummdp_yld2000_2d_am ( a,am )
c-----------------------------------------------------------------------
      implicit none
c
      real*8,intent(in) :: a(8)
c
      real*8,intent(out) :: am(2,3,3)
c
      integer i,j
c-----------------------------------------------------------------------
c
c                                      ---- linear transformation matrix
      am(1,1,1) =  2.0d0*a(1)
      am(1,1,2) = -1.0d0*a(1)
      am(1,1,3) =  0.0d0
      am(1,2,1) = -1.0d0*a(2)
      am(1,2,2) =  2.0d0*a(2)
      am(1,2,3) =  0.0d0
      am(1,3,1) =  0.0d0
      am(1,3,2) =  0.0d0
      am(1,3,3) =  3.0d0*a(7)
c
      am(2,1,1) = -2.0d0*a(3) + 2.0d0*a(4) + 8.0d0*a(5) - 2.0d0*a(6)
      am(2,1,2) =        a(3) - 4.0d0*a(4) - 4.0d0*a(5) + 4.0d0*a(6)
      am(2,1,3) =  0.0d0
      am(2,2,1) =  4.0d0*a(3) - 4.0d0*a(4) - 4.0d0*a(5) +       a(6)
      am(2,2,2) = -2.0d0*a(3) + 8.0d0*a(4) + 2.0d0*a(5) - 2.0d0*a(6)
      am(2,2,3) =  0.0d0
      am(2,3,1) =  0.0d0
      am(2,3,2) =  0.0d0
      am(2,3,3) =  9.0d0*a(8)
c
      do i = 1,3
        do j = 1,3
          am(1,i,j) = am(1,i,j) / 3.0d0
          am(2,i,j) = am(2,i,j) / 9.0d0
        end do
      end do
c
      return
      end subroutine ummdp_yld2000_2d_am
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     calc. barlat-yld2k function x,y,phi
c
      subroutine ummdp_yld2000_2d_xyphi ( s,em,am,x,y,phi )
c-----------------------------------------------------------------------
      implicit none
c
      real*8,intent(in) :: em
      real*8,intent(in) :: s(3)
      real*8,intent(in) :: am(2,3,3)
c
      real*8,intent(out) :: phi(2)
      real*8,intent(out) :: x(2,2),y(2,3)
      
c
      integer i,j,nd
      real*8 a
      real*8 p(2)
c-----------------------------------------------------------------------
c
      p(1) =  1.0d0
      p(2) = -1.0d0
c                                                       ---- {y}=[am]{s}
      call ummdp_utility_clear2 ( y,2,3 )
      do nd = 1,2
        do i = 1,3
          do j = 1,3
            y(nd,i) = y(nd,i) + am(nd,i,j)*s(j)
          end do
        end do
      end do
c                                        ---- {x}=principle value of {y}
      do nd = 1,2
        a = (y(nd,1)-y(nd,2))**2.d0 + 4.0d0*y(nd,3)**2.d0
        a = sqrt(a)
        do i = 1,2
          x(nd,i) = 0.5d0 * (y(nd,1)+y(nd,2)+p(i)*a)
        end do
      end do
c                                                 ---- phi(1) and phi(2)
      nd = 1
      phi(nd) = abs(x(nd,1)-x(nd,2))**em
      nd = 2
      phi(nd) = abs(2.0d0*x(nd,2)+x(nd,1))**em
     1          + abs(2.0d0*x(nd,1)+x(nd,2))**em
c
      return
      end subroutine ummdp_yld2000_2d_xyphi
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     SET 1ST ORDER DERIVATIVE OF PARAMETERS
c
      subroutine ummdp_yld2000_2d_ds1 ( em,am,x,y,phi,dsedphi,dphidx,
     1                                  dxdy,dyds,se )       
c-----------------------------------------------------------------------
      implicit none
c
      real*8 se,em
      real*8 phi(2),dsedphi(2)
      real*8 x(2,2),y(2,3),dphidx(2,2)
      real*8 am(2,3,3),dxdy(2,2,3),dyds(2,3,3)
c
      integer i,j,nd
      real*8 eps,emi,q,a,a0,a1,a2,b0,b1,b2,sgn0,sgn1,sgn2
      real*8 p(2)
c-----------------------------------------------------------------------
c
      eps = 1.0d-16
c
      p(1) =  1.0d0
      p(2) = -1.0d0
      emi = 1.0d0 / em
c                                                          ---- dse/dphi
      q = phi(1) + phi(2)
      if ( q <= 0.0 ) q = 0.0
      do nd = 1,2
        dsedphi(nd) = (0.5d0**emi) * emi * q**(emi-1.0d0)
      end do
c                                                           ---- dphi/dx
      nd = 1
      a0 = x(nd,1) - x(nd,2)
      b0 = abs(a0)
      sgn0 = 0
      if ( b0 >= eps*se ) sgn0 = a0 / b0
      dphidx(nd,1) =  em * b0**(em-1.0d0) * sgn0
      dphidx(nd,2) = -em * b0**(em-1.0d0) * sgn0
c
      nd = 2
      a1 = 2.0d0*x(nd,1) +       x(nd,2)
      a2 =       x(nd,1) + 2.0d0*x(nd,2)
      b1 = abs(a1)
      b2 = abs(a2)
      sgn1 = 0.0d0
      sgn2 = 0.0d0
      if ( b1 >= eps*se ) sgn1 = a1 / b1
      if ( b2 >= eps*se ) sgn2 = a2 / b2
      dphidx(nd,1) = em*( 2.0d0*b1**(em-1.0d0)*sgn1
     1                    + b2**(em-1.0d0)*sgn2 )
      dphidx(nd,2) = em*( b1**(em-1.0d0)*sgn1
     1                    + 2.0d0*b2**(em-1.0d0)*sgn2 )
c
      do nd = 1,2
        a = (y(nd,1)-y(nd,2))*(y(nd,1)-y(nd,2)) + 4.0d0*y(nd,3)*y(nd,3)
        a = sqrt(a)
        if ( a > eps*se ) then
          do j = 1,2
            dxdy(nd,j,1) = 0.5d0 * (1.0d0+p(j)*(y(nd,1)-y(nd,2))/a)
            dxdy(nd,j,2) = 0.5d0 * (1.0d0-p(j)*(y(nd,1)-y(nd,2))/a)
            dxdy(nd,j,3) = 2.0d0 *        p(j)* y(nd,3)         /a
          end do
        else
          do j = 1,2
            dxdy(nd,j,1) = 0.5d0 * (1.0d0+0.0d0)
            dxdy(nd,j,2) = 0.5d0 * (1.0d0-0.0d0)
            dxdy(nd,j,3) = 2.0d0 *        0.0d0
          end do
        end if
      end do
c
      do nd = 1,2
        do i = 1,3
          do j = 1,3
            dyds(nd,i,j) = am(nd,i,j)
          end do
        end do
      end do
c
      return
      end subroutine ummdp_yld2000_2d_ds1
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     SET 2ND ORDER DERIVATIVE OF PARAMETERS
c
      subroutine ummdp_yld2000_2d_ds2 ( phi,x,y,em,d2sedphi2,d2phidx2,
     1                                  d2xdy2,se )
c-----------------------------------------------------------------------
      implicit none
c
      real*8 em,se
      real*8 phi(2)
      real*8 x(2,2),y(2,3),d2sedphi2(2,2)
      real*8 d2phidx2(2,2,2)
      real*8 d2xdy2(2,2,3,3)
c
      integer i,j,m,nd,nd1,nd2,ij
      real*8 eps,emi,q,a
      real*8 p(2)
c-----------------------------------------------------------------------
c
      eps = 1.0d-16
c
      p(1) =  1.0d0
      p(2) = -1.0d0
      emi = 1.0d0 / em
c                                                        ---- d2se/dphi2
      q = phi(1) + phi(2)
      if ( q <= 0.0d0 ) q = 0.0d0
      do nd1 = 1,2
        do nd2 = 1,2
          a = 0.5d0**emi * emi * (emi-1.0d0) * q**(emi-2.0d0)
          d2sedphi2(nd1,nd2) = a
        end do
      end do
c                                                         ---- d2phi/dx2
      nd = 1
      do i = 1,2
        do j = 1,2
          a = (em-1.0d0) * em * (abs(x(nd,1)-x(nd,2)))**(em-2.0d0)
          if ( i /= j ) a = -a
          d2phidx2(nd,i,j) = a
        end do
      end do
      nd = 2
      do i = 1,2
        do j = 1,2
          if ( i == j ) then
            if ( i == 1 ) then
              a = (em-1.0d0) * em
     1            * (4.0d0*(abs(2.0d0*x(nd,1)+x(nd,2)))**(em-2.0d0)
     2               + (abs(x(nd,1)+2.0d0*x(nd,2)))**(em-2.0d0))
            else
              a = (em-1.0d0) * em
     1             * ((abs(2.0d0*x(nd,1)+x(nd,2)))**(em-2.0d0)
     2                + 4.0d0*(abs(x(nd,1)+2.0d0*x(nd,2)))**(em-2.0d0))
            end if
          else
            a = (em-1.0d0) * em
     1           * (2.0d0*(abs(2.0d0*x(nd,1)+x(nd,2)))**(em-2.0d0)
     2              + 2.0d0*(abs(x(nd,1)+2.0d0*x(nd,2)))**(em-2.0d0))
          end if
          d2phidx2(nd,i,j) = a
        end do
      end do
c                                                           ---- d2x/dy2
      do nd = 1,2
        a = (y(nd,1)-y(nd,2))*(y(nd,1)-y(nd,2)) + 4.0d0*y(nd,3)*y(nd,3)       
        if ( a > eps*se ) then
          a = 1.0d0 / sqrt(a**3)
          do m = 1,2
            do i = 1,3
              do j = 1,3
                ij = i*10+j
                if ( (ij == 11) .or. (ij == 22) ) then
                  q = y(nd,3) * y(nd,3)
                else if ( ij == 33 ) then
                  q = (y(nd,1)-y(nd,2)) * (y(nd,1)-y(nd,2))
                else if ( (ij == 12) .or. (ij == 21) ) then
                  q = -y(nd,3) * y(nd,3)
                else if ( (ij == 23) .or. (ij == 32) ) then
                  q = y(nd,3) * (y(nd,1)-y(nd,2))
                else
                  q = -y(nd,3) * (y(nd,1)-y(nd,2))
                end if
                d2xdy2(nd,m,i,j) = 2.0d0 * a * p(m) * q
              end do
            end do
          end do
        else
          do m = 1,2
            do i = 1,3
              do j = 1,3
                d2xdy2(nd,m,i,j) = 0.0d0
              end do
            end do
          end do
        end if
      end do
c
      return
      end subroutine ummdp_yld2000_2d_ds2
c
c
c

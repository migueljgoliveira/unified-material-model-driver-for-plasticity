************************************************************************
c     BBC2005 YIELD FUNCTION AND DERIVATIVES
c
c       doi: 
c
      subroutine ummdp_bbc2005 ( s,se,dseds,d2seds2,nreq,pryld,ndyld )
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
      integer k,i,j,d,e,mm,fact,ii,jj
      real*8 a,b,L,M,N,P,Q,R,nn,fir,las,oo,pp,phi,Al,AA,BB,CC,DD,kk,
     1       lth1_2,lth12,lth2_3,lth23,d2phi_a,d2phi_b,dsedphi,
     2       d2sedphi2,d2pdg2,k_1,k_2,d2pdp2,k_a,k_b,k_c
      real*8 th(3),lth(3),dphidlth(3)
      real*8 dlthds(3,3),d2phidlth2(3,3),d2gds2(3,3),d2lds2(3,3),
     1       d2pds2(3,3)
      real*8 d2lthds2(3,3,3)
c-----------------------------------------------------------------------
c     local variables : symbols in BBC2005 Yield Functions document
c
c     In this subroutine, intermidiate function th(1), th(2) and th(3)
c     of the original paper are represented with Theta(i).
c       th(1) = L*s(1)*M*s(2)
c       th(2) = sqrt((N*s(1)-P*s(2))**2+s(3)**2)
c       th(3) = sqrt((Q*s(1)-R*s(2))**2+s(3)**2)
c     k,a,b,L,M,N,P,Q,R    : material parameter
c     phi                  : phi
c     Al                   : A
c     lth(1) = th(1)**2
c     lth(2) = th(2)**2
c     lth(3) = th(3)**2
c     dlthds(3,3)            :d(lth(i))/d(sigma_x,x,txy)
c     dphidlth(3)            :d(phi)/d(lth(i))
c     dsephi                 :d(se)/d(phi)
c     d2phidlth2(3,3)        :d2(phi)/dg(lth(i))2 ,i =1,2,3
c     d2lthds2(i,3,3)        :d2(lth(i))/d(sigma_x,y,txy)2, i=1,2,3
c
c-----------------------------------------------------------------------
c
c                                            ---- anisotropic parameters
      k = pryld(1+1)
      a = pryld(1+2)
      b = pryld(1+3)
      L = pryld(1+4)
      M = pryld(1+5)
      N = pryld(1+6)
      P = pryld(1+7)
      Q = pryld(1+8)
      R = pryld(1+9)
c
c                                                 ---- equivalent stress
      th(1) = L*s(1)+M*s(2)
      th(2) = sqrt((N*s(1)-P*s(2))**2+s(3)**2)
      th(3) = sqrt((Q*s(1)-R*s(2))**2+s(3)**2)
c
      AA = th(2)**2+th(1)**2
      BB = 2*th(2)*th(1)
      CC = th(2)**2+th(3)**2
      DD = 2*th(2)*th(3)
c
      mm = int(k/2)
      nn = 0.0d0
      oo = 0.0d0
c
      do i = 0,mm
        nn = nn + fact(k)/(fact(k-2*i)*fact(2*i))
     1       *AA**(k-2*i)*BB**(2*i)
        oo = oo + fact(k)/(fact(k-2*i)*fact(2*i))
     1       *CC**(k-2*i)*DD**(2*i)
      end do
c
      fir = 2.0d0 * a * nn
      las = 2.0d0 * b * oo
c
      phi = fir + las
c
      kk = dble(k)
c
      Al = (a*(N+L)**(2*k)+a*(N-L)**(2*k)
     1   +b*(N+Q)**(2*k)+b*(N-Q)**(2*k))**(1.0d0/(2.0d0*kk))
c
      se = phi**(1/(2*kk)) / Al
c
      dseds(:) = 0
      d2seds2(:,:) = 0
c
c                                           ----  1st order differential
      if ( nreq >= 1 ) then
c
        dsedphi = (1.0d0/(2.0d0*kk))*(phi**(1.0d0/(2.0d0*kk)-1.0d0)/Al)
c
        lth(1) = (L*s(1)+M*s(2))**2
        lth(2) = (N*s(1)-P*s(2))**2+s(3)**2
        lth(3) = (Q*s(1)-R*s(2))**2+s(3)**2
c
c                    ----to avoid division by zero or zero of zero power
c 
        if ( lth(1) < 1e-15 * se**2) then
          lth(1) = 1e-15 * se**2
        end if
c
        if ( lth(2) < 1e-15 * se**2) then
          lth(2) = 1e-15 * se**2
        end if
c
        if ( lth(3) < 1e-15 * se**2) then
          lth(3) = 1e-15 * se**2
        end if
c
        lth1_2 = lth(1)+lth(2)
        lth12  = lth(1)*lth(2)
        lth2_3 = lth(2)+lth(3)
        lth23  = lth(2)*lth(3)
c
        if (lth1_2 < 1e-15 * se**2) then
          lth1_2 = 1e-15 * se**2
        end if
c
        if (lth2_3 < 1e-15 * se**2) then
          lth2_3 = 1e-15 * se**2
        end if
c
        dphidlth(:) = 0.0d0
c
        do i = 0,mm  
          dphidlth(1) = dphidlth(1) + fact(k)/(fact(k-2*i)*fact(2*i))*
     1     (2*a*(i*4**i*lth(2)**i*lth(1)**(i-1)*lth1_2
     2     **(k-2*i)+(k-2*i)*4**i*lth12**i*lth1_2**(-2*i+k-1)))
c
          dphidlth(2) = dphidlth(2) + fact(k)/(fact(k-2*i)*fact(2*i))*
     1    (2*a*(i*4**i*lth(2)**(i-1)*lth(1)**i*lth1_2
     2    **(k-2*i)+(k-2*i)*4**i*lth12**i*
     3    lth1_2**(-2*i+k-1))
     4    +2*b*(i*4**i*lth(2)**(i-1)*lth(3)**i*lth2_3
     5    **(k-2*i)+(k-2*i)*4**i*lth23**i*
     6    lth2_3**(-2*i+k-1)))
c
          dphidlth(3) = dphidlth(3) + fact(k)/(fact(k-2*i)*fact(2*i))*
     1     (2*b*(i*4**i*lth(2)**i*lth(3)**(i-1)*lth2_3
     2     **(k-2*i)+(k-2*i)*4**i*lth23**i*lth2_3**(-2*i+k-1)))
        end do
c
        dlthds(1,1) = 2 * L * (M*s(2)+L*s(1))
        dlthds(1,2) = 2 * M * (M*s(2)+L*s(1))
        dlthds(1,3) = 0.0d0
        dlthds(2,1) = 2 * N * (N*s(1)-P*s(2))
        dlthds(2,2) = -2 * P * (N*s(1)-P*s(2))
        dlthds(2,3) = 2 * s(3)
        dlthds(3,1) = 2 * Q * (Q*s(1)-R*s(2))
        dlthds(3,2) = -2 * R * (Q*s(1)-R*s(2))
        dlthds(3,3) = 2 * s(3)
c
        do i = 1,3
          dseds(i) = 0.0d0
          do j = 1,3
            dseds(i) = dseds(i) + dsedphi*dphidlth(j)*dlthds(j,i)
          end do
c         write (150,*) s(1), s(2), s(3), i, dseds(i)
        end do
      end if
c
c                                            ---- 2nd order differential
c
      if ( nreq >= 2 ) then
c
        d2sedphi2 = (1/kk/2.0d00-1)*phi**(1/kk/2.0d0-2)/
     1            (kk*Al)/2.0d0
c
        d2phidlth2(:,:) = 0.0d0
c
        do i = 0,m
c
          k_a = (i-1) * i * 4**i
          k_b = (k-2*i) * i * 4**i
          k_c = (k-2*i-1) * (k-2*i) * 4**i
c
          d2phidlth2(1,1) = d2phidlth2(1,1) + 
     1     fact(k)/(fact(k-2*i)*fact(2*i))
     2     *2*a*(k_a*lth(1)**(i-2)*lth(2)**i*lth1_2**(k-2*i)
     3     +2*(k_b*lth(2)**i*lth(1)**(i-1)*lth1_2**(k-2*i-1))
     4     +k_c*lth12**i*lth1_2**(k-2*i-2))
c
          d2phidlth2(1,2) = d2phidlth2(1,2) + 
     1     fact(k)/(fact(k-2*i)*fact(2*i))
     2     *2*a*(i**2*4**i*lth(1)**(i-1)*lth(2)**(i-1)*lth1_2**(k-2*i)
     3     +k_b*lth(1)**i*lth(2)**(i-1)*lth1_2**(k-2*i-1)
     4     +k_b*lth(1)**(i-1)*lth(2)**i*lth1_2**(k-2*i-1)
     5     +(k-2*i-1)*(k-2*i)*4**i*lth12**i*lth1_2**(k-2*i-2))
c
          d2phidlth2(2,2) = d2phidlth2(2,2) + 
     1     fact(k)/(fact(k-2*i)*fact(2*i))
     2     *(2*b*(k_a*lth(2)**(i-2)*lth(3)**i*lth2_3**(k-2*i)
     3     +2*k_b*lth(3)**i*lth(2)**(i-1)*lth2_3**(k-2*i-1)
     4     +k_c*lth23**i*lth2_3**(k-2*i-2))
     5     +2*a*(k_a*lth(2)**(i-2)*lth(1)**i*lth1_2**(k-2*i)
     6     +2*k_b*lth(1)**i*lth(2)**(i-1)*lth1_2**(k-2*i-1)
     7     +k_c*lth12**i*lth1_2**(k-2*i-2)))
c
          d2phidlth2(2,3) = d2phidlth2(2,3) +
     1     fact(k)/(fact(k-2*i)*fact(2*i))
     2     *2*b*(i**2*4**i*lth(2)**(i-1)*lth(3)**(i-1)*lth2_3**(k-2*i)
     3     +k_b*lth(2)**(i-1)*lth(3)**i*lth2_3**(k-2*i-1)
     4     +k_b*lth(2)**i*lth(3)**(i-1)*lth2_3**(k-2*i-1)
     5     +k_c*lth23**i*lth2_3**(k-2*i-2))
c
          d2phidlth2(3,3) = d2phidlth2(3,3) +
     1     fact(k)/(fact(k-2*i)*fact(2*i))
     2     *2*b*(k_a*lth(2)**i*lth(3)**(i-2)*lth2_3**(k-2*i)
     3     +2*k_b*lth(2)**i*lth(3)**(i-1)*lth2_3**(k-2*i-1)
     4     +k_c*lth23**i*lth2_3**(k-2*i-2))
c
        end do
c
        d2phidlth2(2,1) = d2phidlth2(1,2)
        d2phidlth2(3,2) = d2phidlth2(2,3)
c
        d2lthds2(:,:,:) = 0.0d0
c
        d2lthds2(1,1,1) = 2.0d0 * L**2
        d2lthds2(1,1,2) = 2.0d0 * L * M
        d2lthds2(1,2,1) = d2lthds2(1,1,2)
        d2lthds2(1,2,2) = 2.0d0 * M**2
c
        d2lthds2(2,1,1) = 2.0d0 * N**2
        d2lthds2(2,1,2) = -2.0d0 * N * P
        d2lthds2(2,2,1) = d2lthds2(2,1,2)
        d2lthds2(2,2,2) = 2.0d0 * P**2
        d2lthds2(2,3,3) = 2.0d0
c
        d2lthds2(3,1,1) = 2.0d0 * Q**2
        d2lthds2(3,1,2) = -2.0d0 * Q * R
        d2lthds2(3,2,1) = d2lthds2(3,1,2)
        d2lthds2(3,2,2) = 2.0d0 * R**2
        d2lthds2(3,3,3) = 2.0d0
c
        do i=1,3
          do j=1,3
            d2seds2(i,j)=0.0d0
            do ii=1,3
              do jj=1,3
                d2seds2(i,j) = d2seds2(i,j)
     1                         + d2sedphi2*dphidlth(jj)*dlthds(jj,j)
     2                           * dphidlth(ii)*dlthds(ii,i)
     3                         + dsedphi*d2phidlth2(ii,jj)*dlthds(jj,j)
     4                           * dlthds(ii,i)  
              end do
              d2seds2(i,j)= d2seds2(i,j) 
     1                      + dsedphi*dphidlth(ii)*d2lthds2(ii,i,j)                    
            end do
          end do
        end do
c
      end if
c
      return
      end subroutine ummdp_bbc2005
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     SOLVE FACTORIAL
c
      integer function fact(n) result(m)
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: n
c
      integer i
c-----------------------------------------------------------------------
c
      m = 1
      if ( n>=1 ) then
          do i = 1,n
            m = m * i
          end do
      end if
c
      return
      end function fact
c
c
c
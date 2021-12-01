************************************************************************
c     HILL 1990 YIELD FUNCTION AND DERIVATIVES
c
c       doi: https://doi.org/10.1016/0022-5096(90)90006-P
c
      subroutine ummdp_yield_hill1990 ( s,se,dseds,d2seds2,nreq,pryld,
     1                                  ndyld )                    
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
      integer i,j
      real*8 a,b,tau,sigb,am,syini,sigbtm,alarge,x1,x2,x3,x4,fai1,fai2,
     1       fai3,fai4,fyild,wrk,wrk1,wrk2,wrk3,wrk4
      real*8 v(3),a1(3),dfds1(3),dfds2(3),dfds3(3),dfds4(3),dfds_t(3),
     1       dxds1(3),dxds2(3),dxds3(3),dxds4(3)
      real*8 c(3,3),a2(3,3),a3(3,3),a4(3,3),d2fds1(3,3),d2fds2(3,3),
     1       d2fds3(3,3),d2fds4(3,3),d2fds_t(3,3),dx1dx1(3,3),
     2       dx2dx2(3,3),dx3dx3(3,3),dx4dx4(3,3),df4df3(3,3),
     3       df3df4(3,3),d2xds1(3,3),d2xds2(3,3),d2xds3(3,3),d2xds4(3,3)
      character*32 text
c-----------------------------------------------------------------------
c                                                         ---- variables
c
c     s(3) : stress
c     pryld(ndyld) : material parameter for yield function Hill's 1990
c     pryld(1+1) = a
c     pryld(1+2) = b
c     pryld(1+3) = tau
c     pryld(1+4) = sigb
c     pryld(1+5) = m (=>am)
c
c     se           : equivalent stress
c     dseds(3)     : differential coefficient of first order
c     d2seds2(3,3) : differential coefficient of second order
c
c     a1(3)   : vector to calc for equivalent stress
c     a2(3,3) : matrix to calc for equivalent stress
c     a3(3,3) : matrix to calc for equivalent stress
c     a4(3,3) : matrix to calc for equivalent stress
c
c     dfds1(3) : df1/ds =fai1 1st order derivative by stress
c     dfds2(3) : df2/ds =fai2 1st order derivative by stress
c     dfds3(3) : df3/ds =fai3 1st order derivative by stress
c     dfds4(3) : df4/ds =fai4 1st order derivative by stress
c
c
c     dxds1(3) : dx1/ds =x1 1st order derivative by stress
c     dxds2(3) : dx2/ds =x2 1st order derivative by stress
c     dxds3(3) : dx3/ds =x3 1st order derivative by stress
c     dxds4(3) : dx4/ds =x4 1st order derivative by stress
c
c     d2fds1(3,3)  : d2f1/ds2 =fai1 2nd order derivative by stress
c     d2fds2(3,3)  : d2f2/ds2 =fai2 2nd order derivative by stress
c     d2fds3(3,3)  : d2f3/ds2 =fai3 2nd order derivative by stress
c     d2fds4(3,3)  : d2f4/ds2 =fai4 2nd order derivative by stress
c     d2fds_t(3,3) : d2f/ds2 =fai  2nd order derivative by stress
c-----------------------------------------------------------------------
c
c                                       ---- define a1-matrix, a2-matrix
      data a1/ 1.0d0 , 1.0d0 , 0.0d0 /,
     1     a2/ 1.0d0 ,-1.0d0 , 0.0d0 ,
     2        -1.0d0 , 1.0d0 , 0.d00 ,
     3         0.0d0 , 0.0d0 , 4.0d0 /,
     4     a3/ 1.0d0 , 0.0d0 , 0.0d0 ,
     5         0.0d0 , 1.0d0 , 0.0d0 ,
     6         0.0d0 , 0.0d0 , 2.0d0 /
c                                        ---- set anisotropic parameters
      a    = pryld(1+1)
      b    = pryld(1+2)
      tau  = pryld(1+3)
      sigb = pryld(1+4)
      am   = pryld(1+5)
c
      syini = 1.0d0
      sigbtm = (sigb/tau)**am
      alarge = 1.0d0 + sigbtm - 2.0d0*a + b
c
c                               ---- coef. matrix of material parameters
c                                ---- define a4-matrix consists of a & b
      a4 = 0.0d0
      a4(1,1) = -2.0d0*a + b
      a4(2,2) = 2.0d0*a + b
      a4(1,2) = - b
      a4(2,1) = - b
c                                                              ---- fai1
      x1 = s(1) + s(2)
      fai1 = abs(x1)**am
c                                                              ---- fai2
      call ummdp_utility_mv  ( v,a2,s,3,3 )
      call ummdp_utility_vvs ( x2,s,v,3 )
      fai2 = sigbtm * (x2)**(am/2.0d0)
c                                                              ---- fai3
      call ummdp_utility_mv  ( v,a3,s,3,3 )
      call ummdp_utility_vvs ( x3,s,v,3 )
      fai3 = (x3)**(am/2.0d0-1.0d0)
c                                                              ---- fai4
      call ummdp_utility_mv  ( v,a4,s,3,3 )
      call ummdp_utility_vvs ( x4,s,v,3 )
      fai4 = x4
c                                             ---- yield fuction : fyild
      fyild = fai1 + fai2 + fai3*fai4
c
c                                                 ---- equivalent stress
      se = (fyild/alarge)**(1.0d0/am)
c
      if ( nreq <= 1 ) return

c               ---- 1st order differential coefficient of yield fuction
c     dfdsi(i) : diff. of fai-i(i) with respect to s(j)
c     dxdsi(i) : diff. of x-i(i) with respect to s(j)
c                         1st order differential of x_number_i
c
c                                                          ---- dfai1/ds
      dxds1(1) = 1.0
      dxds1(2) = 1.0
      dxds1(3) = 0.0
c
      wrk = am * (abs(x1)**(am-2))*x1
      do i = 1,3
        dfds1(i) = wrk * dxds1(i)
      end do
c                                                          ---- dfai2/ds
      wrk = sigbtm * (am/2.0) * (x2)**(am/2.0-1.0)
      call ummdp_utility_mv( dxds2,a2,s,3,3 )
      do i = 1,3
        dxds2(i) = 2.0 * dxds2(i)
        dfds2(i) = wrk * dxds2(i)
      end do
c                                                          ---- dfai3/ds
      wrk = (am/2.0-1.0) * (x3)**(am/2.0-2.0)
      call ummdp_utility_mv( dxds3,a3,s,3,3 )
      do i = 1,3
        dxds3(i) = 2.0 * dxds3(i)
        dfds3(i) = wrk * dxds3(i)
      end do
c                                                          ---- dfai4/ds
      call ummdp_utility_mv( dxds4,a4,s,3,3 )
      do i = 1,3
        dxds4(i) = 2.0 * dxds4(i)
        dfds4(i) = dxds4(i)
      end do
c
c        ---- 1st order differential coefficient of yield fuction result
c                                     ---- dfai/ds()= result = dfds_t(i)
      do i = 1,3
        dfds_t(i) = dfds1(i) + dfds2(i) + dfds3(i)*fai4 + fai3*dfds4(i)
      end do
c           ---- 1st order differential coefficient of equivalent stress
      wrk = (abs(fyild/alarge))**(1.0/am-1.0) / (am*alarge)
      do i = 1,3
        dseds(i) = wrk * dfds_t(i)
      end do
c
c
      if ( nreq <= 2 ) return
c
c            --- 2st order differential coefficient of equivalent stress
c                                                   with respect to s(j)
c-----------------------------------------------------------------------
c     dfds_t(3,3) : 2nd order differ of fai by s(j) & s(k)
c     df2ds1(3,3) : 2nd order diff. of  s(j) & s(k)
c-----------------------------------------------------------------------
c
c                                                        ---- d2fai1/ds2
      wrk = am * (am-1.0) * (abs(x1))**(am-2.0)
      do i = 1,3
        do j = 1,3
          d2fds1(i,j) = wrk * dxds1(i) * dxds1(j)
        end do
      end do
c                                                        ---- d2fai2/ds2
      wrk1 = sigbtm * (am/2.0)
      if ( abs(x2) < 1e-10 ) x2 = 1e-10
      wrk2 = (am/2.0-1.0) * (x2**(am/2.0-2.0))
      wrk3 = x2**(am/2.0-1.0)
      wrk2 = wrk1 * wrk2
      wrk3 = wrk1 * wrk3
c             ---- make [ dx2 * dx2(t) ] & [d2x/ds2] & make [d2fai2/ds2]
      do i = 1,3
        do j = 1,3
          dx2dx2(i,j) = dxds2(i) * dxds2(j)
          d2xds2(i,j) = 2.0 * a2(j,i)
          d2fds2(i,j) = wrk2 * dx2dx2(i,j) + wrk3 * d2xds2(i,j)
        end do
      end do
c                                   ---- d2fai3/ds2   make   d2fds3(i,j)
      wrk1 = am/2.0 - 1.0
      wrk2 = (am/2.0-2.0) * (x3**(am/2.0-3.0))
      wrk3 = x3**(am/2.0-2.0)
      wrk2 = wrk1 * wrk2
      wrk3 = wrk1 * wrk3
c                                   ---- [d2x3/ds2] &  make [d2fai3/ds2]
      do i = 1,3
        do j = 1,3
          dx3dx3(i,j) = dxds3(i) * dxds3(j)
          d2xds3(i,j) = 2.0 * a3(j,i)
          d2fds3(i,j) = wrk2 * dx3dx3(i,j) + wrk3 * d2xds3(i,j)
        end do
      end do
c                                                 ---- [d2fai3/ds2]*fai4
      do i = 1,3
        do j = 1,3
            d2fds3(i,j) = d2fds3(i,j) * fai4
          end do
        end do
c                                          ---- [dfai4/ds]*[dfai3/ds](T)
        do i = 1,3
          do j = 1,3
            df4df3(i,j) = dfds4(i) * dfds3(j)
          end do
        end do
c                                                        ---- d2fai4/ds2
c                                                 ---- make [d2fai3/ds2]
        do i = 1,3
          do j = 1,3
            d2fds4(i,j) = 2.0 * a4(i,j)
          end do
        end do
c
c        ---- 2nd order differential coefficient of yield fuction result
c                                  ---- d2fai/ds2()= result = d2fds_t(i)
c
        do i = 1,3
          do j = 1,3
            d2fds_t(i,j) = d2fds1(i,j) + d2fds2(i,j) + 
     1                     d2fds3(i,j)*fai4 + df4df3(i,j) + 
     2                     df4df3(j,i) + fai3*d2fds4(i,j)
          end do
        end do
c
c           ---- 2nd order differential coefficient of equivalent stress
c                                                              by stress
        wrk1 = 1.0/(am*alarge)
        wrk2 = (1.0/am-1.0)/alarge
        wrk3 = (fyild/alarge)**(1.0/am-2.0)
        wrk4 = (fyild/alarge)**(1.0/am-1.0)
        wrk2 =  wrk1 * wrk2 * wrk3
        wrk4 =  wrk1 * wrk4
c
        do i = 1,3
          do j = 1,3
            d2seds2(i,j) = wrk2 * dfds_t(i) * dfds_t(j) + 
     1                     wrk4 * d2fds_t(i,j)
          end do
        end do
c
      return
      end subroutine ummdp_yield_hill1990
c
c
c

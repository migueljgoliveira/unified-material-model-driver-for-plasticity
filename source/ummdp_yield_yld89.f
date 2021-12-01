************************************************************************
c     YLD89 YIELD FUNCTION AND DERIVATIVES
c
c       doi: 
c
      subroutine ummdp_yield_yld89 ( s,se,dseds,d2seds2,nreq,pryld,
     1                               ndyld )
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
      real*8 delta,se0,ta1,ta2,ta3,tb1,tb2,tb3,asy
      real*8 s0(3),dsedsz(3)
      real*8 d2seds2z(3,3)
c-----------------------------------------------------------------------
c
      call ummdp_yield_yld89_branch ( s,se,dseds,d2seds2,nreq,
     1                                pryld,ndyld )
c
      if ( nreq <= 1 ) return
c
c                                         ---- numerical differentiation
      delta = 1.0d-6 * se
c
      do j = 1,3
        s0(j) = s(j)
      end do
c
      do i = 1,3
        s0(i) = s(i) + delta
        call ummdp_yield_yld89_branch ( s0,se0,dsedsz,d2seds2z,1,
     1                                  pryld,ndyld )
        ta1 = dsedsz(1)
        ta2 = dsedsz(2)
        ta3 = dsedsz(3)
c
        s0(i) = s(i) - delta
        call ummdp_yield_yld89_branch ( s0,se0,dsedsz,d2seds2z,1,
     2                                  pryld,ndyld )
        tb1 = dsedsz(1)
        tb2 = dsedsz(2)
        tb3 = dsedsz(3)
c
        s0(i) = s(i)   ! re-set s0(i) for next loop
c
        d2seds2(1,i) = (ta1-tb1) / (2.0d0*delta)
        d2seds2(2,i) = (ta2-tb2) / (2.0d0*delta)
        d2seds2(3,i) = (ta3-tb3) / (2.0d0*delta)
      end do
c
c                                         ---- d2seds2 must be symmetric
      do i = 1,3
        do j = i+1,3
          asy = 0.5d0*(d2seds2(i,j) + d2seds2(j,i))
          d2seds2(i,j) = asy
          d2seds2(j,i) = asy
        end do
      end do
c
      return
      end subroutine ummdp_yield_yld89
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     BRANCH OF YLD89
c
      subroutine ummdp_yield_yld89_branch ( s,se,dseds,d2seds2,nreq,
     1                                      pryld,ndyld )
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
      integer i,j,k,l
      real*8 pM,a,h,p,pK1,pK3,pK2,pK22,f1,f2,f3,f,DpK22,dsedf,d2sedf2,
     1       sc1,sc2,sc3
      real*8 dfdK(2),df1dK(2),df2dK(2),df3dK(2),u(2),v(2),w(2),x(2)
      real*8 dKds(2,3),d2fdK2(2,2),dfdKdfdK(2,2),d2K1ds2(3,3),
     1       d2K2ds2(3,3),d2f1dK2(2,2),d2f2dK2(2,2),d2f3dK2(2,2)
c-----------------------------------------------------------------------
c
c                                                        ---- parameters
      pM = pryld(1+1)
      a  = pryld(1+2)
      h  = pryld(1+3)
      p  = pryld(1+4)                                   
c                                                         ---- clear end
c
c     K1
c
      pK1 = (s(1)+h*s(2)) * 0.5d0

c     K2
c     K22=K2^2
      pK3 = (s(1)-h*s(2)) * 0.5d0
      pK22 = pK3*pK3 + p*p*s(3)*s(3)
      pK2 = pK22 ** 0.5d0
c
      f1 = a * abs(pK1+pK2)**pM
      f2 = a * abs(pK1-pK2)**pM
      f3 = (2.0d0-a) * abs(2.0d0*pK2)**pM
c
      f = f1 + f2 + f3
c
      se = (0.5d0*f)**(1.0d0/pM)
c
c                                                             ---- dseds
      if ( nreq >= 1 ) then
c                                                         ---- dKds(2,3)
        dKds(1,1) = 0.5d0
        dKds(1,2) = h * 0.5d0
        dKds(1,3) = 0.0d0
c
        if ( pK2 == 0.0d0 ) then
          dKds(2,1) = 0.0d0
          dKds(2,2) = 0.0d0
          dKds(2,3) = 0.0d0
          if ( s(3) == 0.0d0 ) dKds(2,3) = p * p
        else
          dKds(2,1) = pK3 / (2.0d0*pK2)
          dKds(2,2) = -h*pK3 / (2.0d0*pK2)
          dKds(2,3) = p*p*s(3) / pK2
        end if
c                                                  ---- d2K1ds2(3,3)=0.0
        d2K1ds2 = 0.0d0
c                                                      ---- d2K2ds2(3,3)
        if ( pK2 == 0.0d0 ) then
          DpK22 = 1.0d-32
c
          d2K2ds2(1,1) = (DpK22**(-0.5d0)-pK3*DpK22**(-1.5d0)) / 4.0d0
          d2K2ds2(1,2) = -h * d2K2ds2(1,1)
          d2K2ds2(1,3) = -p*p*s(3)*pK3/2.0d0*DpK22**(-1.5d0)
c
          d2K2ds2(2,1) = d2K2ds2(1,2)
          d2K2ds2(2,2) = h * h * d2K2ds2(1,1)
          d2K2ds2(2,3) = -h * d2K2ds2(1,3)
c
          d2K2ds2(3,1) = d2K2ds2(1,3)
          d2K2ds2(3,2) = d2K2ds2(2,3)
          d2K2ds2(3,3) = p*p*(DpK22**(-0.5d0) - 
     1                   p*p*s(3)*s(3)*DpK22**(-1.5d0))
        else
          d2K2ds2(1,1) = (pK22**(-0.5d0)-pK3*pK22**(-1.5d0))/4.0d0
          d2K2ds2(1,2) = -h*d2K2ds2(1,1)
          d2K2ds2(1,3) = -p*p*s(3)*pK3/2.0d0*pK22**(-1.5d0)
c
          d2K2ds2(2,1) = d2K2ds2(1,2)
          d2K2ds2(2,2) = h * h * d2K2ds2(1,1)
          d2K2ds2(2,3) = -h * d2K2ds2(1,3)
c
          d2K2ds2(3,1) = d2K2ds2(1,3)
          d2K2ds2(3,2) = d2K2ds2(2,3)
          d2K2ds2(3,3) = p*p
     1                   * (pK22**(-0.5d0)-p*p*s(3)*s(3)*pK22**(-1.5d0))
        end if
c                                                              ---- dfdK
c                                                             ---- df1dK
        do i = 1,2
          df1dK(i) = a*pM*(pK1+pK2)*abs(pK1+pK2)**(pM-2.0d0)
        end do
c                                                             ---- df2dK
        df2dK(1) = a*pM*(pK1-pK2)*abs(pK1-pK2)**(pM-2.0d0)
        df2dK(2) = -a*pM*(pK1-pK2)*abs(pK1-pK2)**(pM-2.0d0)
c                                                             ---- df3dK
        df3dK(1) = 0.0d0
        df3dK(2) = 2.0d0*pM*(2.0d0-a)*(2.0d0*pK2) * 
     1              abs(2.0d0*pK2)**(pM-2.0d0)
c                                            ---- dfdK=df1dK+df2dK+df3dK
        dfdK = 0.0d0
        do i = 1,2
          dfdK(i) = df1dK(i)+df2dK(i)+df3dK(i)
        end do
c                                                             ---- dsedf
        dsedf = 0.0d0
        dsedf = (0.5d0*f) ** ((1.0d0-pM)/pM)
        dsedf = dsedf / (2.0d0*pM)
c                                                             ---- dseds
        dseds = 0.0d0
        do i = 1,3
          do j = 1,2
            dseds(i) = dseds(i)+dfdK(j)*dKds(j,i)
          end do
          dseds(i) = dseds(i)*dsedf
        end do
      end if
c                                                           ---- d2seds2
      if ( nreq >= 2 ) then
c                                                            ---- d2fdK2
c                                                           ---- d2f1dK2
        do i = 1,2
          do j = 1,2
            d2f1dK2(i,j) = a*pM*(pM-1.0d0)*abs(pK1+pK2)**(pM-2.0d0)
          end do
        end do
c                                                           ---- d2f2dK2
        d2f2dK2(1,1) = a*pM*(pM-1.0d0)*abs(pK1-pK2)**(pM-2.0d0)
        d2f2dK2(1,2) = -a*pM*(pM-1.0d0)*abs(pK1-pK2)**(pM-2.0d0)
        d2f2dK2(2,1) = d2f2dK2(1,2)
        d2f2dK2(2,2) = d2f2dK2(1,1)
c                                                           ---- d2f3dK2
        d2f3dK2(1,1) = 0.0d0
        d2f3dK2(1,2) = 0.0d0
        d2f3dK2(2,1) = 0.0d0
        d2f3dK2(2,2) = 4.0d0*pM*(pM-1.0d0)*(2.0d0-a)
     1                                 *abs(2.0d0*pK2)**(pM-2.0d0)
c
c                                    ---- d2fdK2=d2f1dK2+d2f2dK2+d2f3dK2
        do i = 1,2
          do j = 1,2
            d2fdK2(i,j) = d2f1dK2(i,j)+d2f2dK2(i,j)+d2f3dK2(i,j)
          end do
        end do
c                                                     ---- dfdKdfdK(2,2)
        do i = 1,2
          do j = 1,2
            dfdKdfdK(i,j) = dfdK(i)*dfdK(j)
          end do
        end do
c                                                           ---- d2sedf2
        d2sedf2 = (0.5d0*f)**((1.0d0-2.0d0*pM)/pM)
        d2sedf2 = d2sedf2*(1.0d0-pM) / (4.0d0*pM*pM)
c
        do i = 1,3
          do j = 1,3
            w = 0.0d0
            x = 0.0d0
            do k = 1,2
              do l = 1,2
                w(k) = w(k) + dfdKdfdK(k,l)*dKds(l,j)
                x(k) = x(k) + d2fdK2(k,l)*dKds(l,j)
              end do
            end do
            sc1 = 0.0d0
            sc2 = 0.0d0
            sc3 = 0.0d0
            do k = 1,2
              sc1 = sc1 + dKds(k,i)*w(k)
              sc2 = sc2 + dKds(k,i)*x(k)
            end do
            sc3 = dfdK(1)*d2K1ds2(j,i) + dfdK(2)*d2K2ds2(j,i)
            d2seds2(i,j) = d2sedf2*sc1 + dsedf*sc2 + d2sedf2*sc3
          end do
        end do
      end if
c
      return
      end subroutine ummdp_yield_yld89_branch
c
c
c

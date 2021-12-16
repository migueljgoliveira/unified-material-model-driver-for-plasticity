************************************************************************
*
*     KINEMATIC HARDENING LAWS
*
************************************************************************
c
c      0 : No Kinematic Hardening
c      1 : Prager
c      2 : Ziegler
c      3 : Armstrong & Frederick (1966)
c      4 : Chaboche (1979)
c      5 : Chaboche (1979) - Ziegler
c      6 : Yoshida-Uemori
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     KINEMATIC HARDENING LAW
c
      subroutine ummdp_kinematic ( vk,dvkdp,dvkds,dvkdx,dvkdxt,p,s,x,xt,
     1                             nttl,nnrm,nshr,mxpbs,npbs,prkin,
     2                             ndkin,pryld,ndyld )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer nttl,nnrm,nshr,mxpbs,npbs,ndkin,ndyld
      real*8 p
      real*8 prkin(ndkin),pryld(ndyld),s(nttl),xt(nttl)
      real*8 vk(npbs,nttl),dvkdp(npbs,nttl),x(mxpbs,nttl)
      real*8 dvkds(npbs,nttl,nttl),dvkdxt(npbs,nttl,nttl)
      real*8 dvkdx(npbs,npbs,nttl,nttl)
c
      integer i,j,k,l
      integer ntkin
c-----------------------------------------------------------------------
c
      ntkin = nint(prkin(1))
c                                                        ---- initialize
      do i = 1,npbs
        do j = 1,nttl
          vk(i,j) = 0.0
          dvkdp(i,j) = 0.0
          do k = 1,nttl
            dvkds(i,j,k) = 0.0
            dvkdxt(i,j,k) = 0.0
            do l = 1,npbs
              dvkdx(i,l,j,k) = 0.0
            end do
          end do
        end do
      end do
c
      select case ( ntkin )
c
      case ( 0 )                                ! No Kinematic Hardening
        return
c
      case ( 1 )                                                ! Prager
        call ummdp_kinematic_prager ( vk,dvkdp,dvkds,dvkdx,dvkdxt,p,s,x,
     1                                xt,nttl,nnrm,nshr,mxpbs,npbs,
     2                                prkin,ndkin,pryld,ndyld )
c
      case ( 2 )                                               ! Ziegler
        call ummdp_kinematic_ziegler ( vk,dvkdp,dvkds,dvkdx,dvkdxt,p,s,
     1                                 x,xt,nttl,nnrm,nshr,mxpbs,npbs,
     2                                 prkin,ndkin,pryld,ndyld )
c
      case ( 3 )                                 ! Armstrong & Frederick
        call ummdp_kinematic_armstrong ( vk,dvkdp,dvkds,dvkdx,dvkdxt,p,
     1                                   s,x,xt,nttl,nnrm,nshr,mxpbs,
     2                                   npbs,prkin,ndkin,pryld,ndyld )
c
      case ( 4 )                                            ! Chaboche I
        call ummdp_kinematic_chabocheI ( vk,dvkdp,dvkds,dvkdx,dvkdxt,p,
     1                                  s,x,xt, nttl,nnrm,nshr,mxpbs,
     2                                  npbs,prkin,ndkin,pryld,ndyld )
c
      case ( 5 )                                           ! Chaboche II
        call ummdp_kinematic_chabocheII ( vk,dvkdp,dvkds,dvkdx,
     1                                          dvkdxt,p,s,x,xt,nttl,
     2                                          nnrm,nshr,mxpbs,npbs,
     3                                          prkin,ndkin,pryld,
     4                                          ndyld )
c
      case ( 6 )                                        ! Yoshida-Uemori
        call ummdp_kinematic_yoshida_uemori ( vk,dvkdp,dvkds,dvkdx,
     1                                        dvkdxt,p,s,x,xt,nttl,nnrm,
     2                                        nshr,mxpbs,npbs,prkin,
     3                                        ndkin,pryld,ndyld )
c
      case default
        write (6,*) 'still not be supported. ntkin=',ntkin
        call ummdp_exit ( 204 )
      end select
c
      return
      end subroutine ummdp_kinematic
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     CALCULATE 1ST AND 2ND DERIVATIVES FOR KINEMATIC HARDENING LAWS
c
      subroutine ummdp_kinematic_dseds ( eta,seta,dseds,d2seds2,nttl,
     1                                   nnrm,nshr,pryld,ndyld )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer nttl,nnrm,nshr,ndyld
      real*8 seta
      real*8 eta(nttl),dseds(nttl),pryld(ndyld)
      real*8 d2seds2(nttl,nttl)
c
      integer i,j
      real*8 em1,em2,en11,en12,en13,en21,en22,en23
c-----------------------------------------------------------------------
c
c                    ---- dseds and d2seds2 for plastic strain increment
      call ummdp_yield ( seta,dseds,d2seds2,2,eta,nttl,nnrm,nshr,
     1                   pryld,ndyld )
c
c                   ---- engineering shear strain -> tensor shear strain
      do i = nnrm+1,nttl
        dseds(i) = 0.5d0*dseds(i)
        do j = 1,nttl
          d2seds2(i,j) = 0.5d0*d2seds2(i,j)
        end do
      end do
c                                          ---- for plane stress problem
      if ( nnrm == 2 ) then
        em1 = dseds(1)
        em2 = dseds(2)
        dseds(1) = dseds(1) + em1 + em2
        dseds(2) = dseds(2) + em2 + em1
        en11 = d2seds2(1,1)
        en12 = d2seds2(1,2)
        en13 = d2seds2(1,3)
        en21 = d2seds2(2,1)
        en22 = d2seds2(2,2)
        en23 = d2seds2(2,3)
        d2seds2(1,1) = d2seds2(1,1) + en11 + en21
        d2seds2(1,2) = d2seds2(1,2) + en12 + en22
        d2seds2(1,3) = d2seds2(1,3) + en13 + en23
        d2seds2(2,1) = d2seds2(2,1) + en21 + en11
        d2seds2(2,2) = d2seds2(2,2) + en22 + en12
        d2seds2(2,3) = d2seds2(2,3) + en23 + en13
      end if
c
      return
      end subroutine ummdp_kinematic_dseds
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     PRAGER KINEMATIC HARDENING LAW
c
      subroutine ummdp_kinematic_prager ( vk,dvkdp,dvkds,dvkdx,dvkdxt,
     1                                    p,s,x,xt,nttl,nnrm,nshr,mxpbs,
     2                                    npbs,prkin,ndkin,pryld,ndyld )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer nttl,nnrm,nshr,mxpbs,npbs,ndkin,ndyld
      real*8 p
      real*8 s(nttl),xt(nttl),prkin(ndkin),pryld(ndyld)
      real*8 vk(npbs,nttl),dvkdp(npbs,nttl),x(mxpbs,nttl)
      real*8 dvkds(npbs,nttl,nttl),dvkdxt(npbs,nttl,nttl)
      real*8 dvkdx(npbs,npbs,nttl,nttl)
c
      integer i,j,n
      real*8 c,seta,dcdp
      real*8 eta(nttl),dseds(nttl)
      real*8 d2seds2(nttl,nttl)
c-----------------------------------------------------------------------
c
      c = prkin(2)/3.0d0*2.0d0
c
      do i = 1,nttl
        eta(i) = s(i) - xt(i)
      end do
c
      call ummdp_kinematic_dseds ( eta,seta,dseds,d2seds2,nttl,nnrm,
     1                             nshr,pryld,ndyld )
c
      n = 1
c
      do i = 1,nttl
        vk(n,i) = c * dseds(i)
      end do
c
      dcdp = 0.0d0
      do i = 1,nttl
        dvkdp(n,i) = dcdp * dseds(i)
      end do
c
      do i = 1,nttl
        do j = 1,nttl
          dvkds(n,i,j) = c * d2seds2(i,j)
          dvkdx(n,n,i,j) = -c * d2seds2(i,j)
          dvkdxt(n,i,j) = 0.0
        end do
      end do
c
      return
      end subroutine ummdp_kinematic_prager
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     ZIEGLER KINEMATIC HARDENING LAW
c
      subroutine ummdp_kinematic_ziegler ( vk,dvkdp,dvkds,dvkdx,dvkdxt,
     1                                     p,s,x,xt,nttl,nnrm,nshr,
     2                                     mxpbs,npbs,prkin,ndkin,
     3                                     pryld,ndyld )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer nttl,nnrm,nshr,mxpbs,npbs,ndkin,ndyld
      real*8 p
      real*8 s(nttl),xt(nttl),prkin(ndkin),pryld(ndyld)
      real*8 vk(npbs,nttl),dvkdp(npbs,nttl),x(mxpbs,nttl)
      real*8 dvkds(npbs,nttl,nttl),dvkdxt(npbs,nttl,nttl)
      real*8 dvkdx(npbs,npbs,nttl,nttl)
c
      integer i,j
      integer n
      real*8 c,dcdp
      real*8 eta(nttl),dseds(nttl)
      real*8 d2seds2(nttl,nttl),am(nttl,nttl)
c-----------------------------------------------------------------------
c
      c = prkin(2)
c
      do i = 1,nttl
        eta(i) = s(i) - xt(i)
      end do
c
      n = 1
      do i = 1,nttl
        vk(n,i) = c * eta(i)
      end do
c
      dcdp = 0.0
      do i = 1,nttl
        dvkdp(n,i) = dcdp * eta(i)
      end do
c
      call ummdp_utility_setunitm ( am,nttl )
      do i = 1,nttl
        do j = 1,nttl
          dvkds(n,i,j) = c * am(i,j)
          dvkdx(n,n,i,j) = -c * am(i,j)
          dvkdxt(n,i,j) = 0.0
        end do
      end do
c
      return
      end subroutine ummdp_kinematic_ziegler
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     ARMSTRONG-FREDERICK KINEMATIC HARDENING LAW
c
      subroutine ummdp_kinematic_armstrong ( vk,dvkdp,dvkds,dvkdx,
     1                                       dvkdxt,p,s,x,xt,nttl,nnrm,
     2                                       nshr,mxpbs,npbs,prkin,
     3                                       ndkin,pryld,ndyld )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer nttl,nnrm,nshr,mxpbs,npbs,ndkin,ndyld
      real*8 p
      real*8 s(nttl),xt(nttl),prkin(ndkin),pryld(ndyld)
      real*8 vk(npbs,nttl),dvkdp(npbs,nttl),x(mxpbs,nttl)
      real*8 dvkds(npbs,nttl,nttl),dvkdxt(npbs,nttl,nttl)
      real*8 dvkdx(npbs,npbs,nttl,nttl)
c
      integer i,j
      integer n
      real*8 c,g,seta,dcdp,dgdp
      real*8 eta(nttl),dseds(nttl)
      real*8 d2seds2(nttl,nttl),am(nttl,nttl)
c-----------------------------------------------------------------------
c
      c = prkin(1+1)/3.0d0*2.0d0
      g = prkin(1+2)
c
      do i = 1,nttl
        eta(i) = s(i) - xt(i)
      end do
c
      call ummdp_kinematic_dseds ( eta,seta,dseds,d2seds2,nttl,nnrm,
     1                             nshr,pryld,ndyld )
c
      n = 1
      do i = 1,nttl
        vk(n,i) = c*dseds(i) - g*xt(i)
      end do
c
      dcdp = 0.0d0
      dgdp = 0.0d0
      do i = 1,nttl
        dvkdp(n,i) = dcdp*dseds(i) - dgdp*xt(i)
      end do
c
      call ummdp_utility_setunitm ( am,nttl )
      do i = 1,nttl
        do j = 1,nttl
          dvkds(n,i,j) = c * d2seds2(i,j)
          dvkdx(n,n,i,j) = -c*d2seds2(i,j) - g*am(i,j)
          dvkdxt(n,i,j) = 0.0
        end do
      end do
c
      return
      end subroutine ummdp_kinematic_armstrong
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     CHABOCHE I KINEMATIC HARDENING LAW
c
      subroutine ummdp_kinematic_chabocheI ( vk,dvkdp,dvkds,dvkdx,
     1                                       dvkdxt,p,s,x,xt,nttl,nnrm,
     2                                       nshr,mxpbs,npbs,prkin,
     3                                       ndkin,pryld,ndyld )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer nttl,nnrm,nshr,mxpbs,npbs,ndkin,ndyld
      real*8 p
      real*8 s(nttl),xt(nttl),prkin(ndkin),pryld(ndyld)
      real*8 vk(npbs,nttl),dvkdp(npbs,nttl),x(mxpbs,nttl)
      real*8 dvkds(npbs,nttl,nttl),dvkdxt(npbs,nttl,nttl)
      real*8 dvkdx(npbs,npbs,nttl,nttl)
c
      integer i,j,n
      integer n0
      real*8 seta,c,g,dcdp,dgdp
      real*8 eta(nttl),dseds(nttl)
      real*8 d2seds2(nttl,nttl),am(nttl,nttl)
c-----------------------------------------------------------------------
c
      do i = 1,nttl
        eta(i) = s(i) - xt(i)
      end do
c
      call ummdp_kinematic_dseds ( eta,seta,dseds,d2seds2,nttl,nnrm,
     1                             nshr,pryld,ndyld )
c
      call ummdp_utility_setunitm ( am,nttl )
      do n = 1,npbs
        n0 = 1 + (n-1)*2
        c = prkin(1+n0+1)/3.0d0*2.0d0
        g = prkin(1+n0+2)
        do i = 1,nttl
          vk(n,i) = c*dseds(i) - g*x(n,i)
        end do
        dcdp = 0.0d0
        dgdp = 0.0d0
        do i = 1,nttl
          dvkdp(n,i) = dcdp*dseds(i) - dgdp*x(n,i)
        end do
        do i = 1,nttl
          do j = 1,nttl
            dvkds(n,i,j) = c * d2seds2(i,j)
            dvkdx(n,n,i,j) = -g * am(i,j)
            dvkdxt(n,i,j) = -c * d2seds2(i,j)
          end do
        end do
      end do
c
      return
      end subroutine ummdp_kinematic_chabocheI
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     CHABOCHE II KINEMATIC HARDENING LAW
c
      subroutine ummdp_kinematic_chabocheII ( vk,dvkdp,dvkds,dvkdx,
     1                                        dvkdxt,p,s,x,xt,nttl,nnrm,
     2                                        nshr,mxpbs,npbs,prkin,
     3                                        ndkin,pryld,ndyld )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer nttl,nnrm,nshr,mxpbs,npbs,ndkin,ndyld
      real*8 p
      real*8 s(nttl),xt(nttl),prkin(ndkin),pryld(ndyld)
      real*8 vk(npbs,nttl),dvkdp(npbs,nttl),x(mxpbs,nttl)
      real*8 dvkds(npbs,nttl,nttl),dvkdxt(npbs,nttl,nttl)
      real*8 dvkdx(npbs,npbs,nttl,nttl)
c
      integer i,j,n
      integer n0
      real*8 seta,c,g,dcdp,dgdp
      real*8 eta(nttl),dseds(nttl)
      real*8 d2seds2(nttl,nttl),am(nttl,nttl)
c-----------------------------------------------------------------------
c
      do i = 1,nttl
        eta(i) = s(i) - xt(i)
      end do
c
      call ummdp_kinematic_dseds ( eta,seta,dseds,d2seds2,nttl,nnrm,
     1                             nshr,pryld,ndyld )
c
      call ummdp_utility_setunitm ( am,nttl )
      do n = 1,npbs
        n0 = 1 + (n-1)*2
        c = prkin(1+n0+1)
        g = prkin(1+n0+2)
        do i = 1,nttl
          vk(n,i) = (c/seta)*eta(i) - g*x(n,i)
        end do
        dcdp = 0.0d0
        dgdp = 0.0d0
        do i = 1,nttl
          dvkdp(n,i) = (dcdp/seta)*eta(i) - dgdp*x(n,i)
        end do
        do i = 1,nttl
          do j = 1,nttl
            dvkds(n,i,j) = (c/seta) * am(i,j)
            dvkdx(n,n,i,j) = -g * am(i,j)
            dvkdxt(n,i,j) = -(c/seta) * am(i,j)
          end do
        end do
      end do
c
      return
      end subroutine ummdp_kinematic_chabocheII
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     YOSHIDA-UEMORI KINEMATIC HARDENING LAW
c
      subroutine ummdp_kinematic_yoshida_uemori ( vk,dvkdp,dvkds,dvkdx,
     1                                            dvkdxt,p,s,x,xt,nttl,
     2                                            nnrm,nshr,mxpbs,npbs,
     3                                            prkin,ndkin,pryld,
     4                                            ndyld )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer nttl,nnrm,nshr,mxpbs,npbs,ndkin,ndyld
      real*8 p
      real*8 s(nttl),xt(nttl),prkin(ndkin),pryld(ndyld)
      real*8 vk(npbs,nttl),dvkdp(npbs,nttl),x(mxpbs,nttl)
      real*8 dvkds(npbs,nttl,nttl),dvkdxt(npbs,nttl,nttl)
      real*8 dvkdx(npbs,npbs,nttl,nttl)
c
      integer i,j
      integer n
      real*8 pc,py,pa,pk,pb,seta
      real*8 eta(nttl),dseds(nttl)
      real*8 d2seds2(nttl,nttl),am(nttl,nttl)
c-----------------------------------------------------------------------
c
      pc = prkin(1+1)
      py = prkin(1+2)
      pa = prkin(1+3)
      pk = prkin(1+4)
      pb = prkin(1+5)
c
      do i = 1,nttl
        eta(i) = s(i) - xt(i)
      end do
c
      call ummdp_kinematic_dseds ( eta,seta,dseds,d2seds2,nnrm,nshr,
     1                             nttl,pryld,ndyld )
      call ummdp_utility_setunitm ( am,nttl )
c
      n = 1
      do i = 1,nttl
        vk(n,i) = pc*((pa/py)*eta(i) - sqrt(pa/seta)*x(n,i))
      end do
      do i = 1,nttl
        dvkdp(n,i) = 0.0
      end do
      do i = 1,nttl
        do j = 1,nttl
          dvkds(n,i,j) = pc*pa/py*am(i,j)
          dvkdxt(n,i,j) = -pc*pa/py*am(i,j)
          dvkdx(n,n,i,j) = pc*sqrt(pa)*
     1                    ( -am(i,j)/sqrt(seta)
     2                       + x(n,i)*dseds(j)/(2.0d0*seta**(1.5d0)) )
        end do
      end do
c
      n = 2
      do i = 1,nttl
        vk(n,i) = pk*(2.0d0/3.0d0*pb*dseds(i) - x(n,i))
      end do
      do i = 1,nttl
        dvkdp(n,i) = 0.0
      end do
      do i = 1,nttl
        do j = 1,nttl
          dvkds(n,i,j) = 2.0d0/3.0d0*pb + pk*d2seds2(i,j)
          dvkdxt(n,i,j) = -2.0d0/3.0d0*pb + pk*d2seds2(i,j)
          dvkdx(n,n,i,j) = -pk * am(i,j)
        end do
      end do
c
      return
      end subroutine ummdp_kinematic_yoshida_uemori
c
c
c
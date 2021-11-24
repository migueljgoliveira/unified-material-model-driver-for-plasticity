************************************************************************
c     BBC2008 YIELD FUNCTION AND DERIVATIVES
c
      subroutine jancae_bbc2008 (s,se,dseds,d2seds2,nreq,
     1                           pryld,ndyld)
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
      integer nds,ndk
c-----------------------------------------------------------------------
c
      nds = nint(pryld(2))
      ndk = nint(pryld(3))
c
      call jancae_bbc2008_core ( s,se,dseds,d2seds2,nreq,
     1                           pryld,ndyld,nds,ndk )
c
      return
      end subroutine jancae_bbc2008
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine jancae_bbc2008_core ( s,se,dseds,d2seds2,nreq,
     1                                 pryld,ndyld,sp,kp )
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: nreq,ndyld,sp,kp
      real*8 ,intent(in) :: s(3),pryld(ndyld)
c
      real*8,intent(out) :: se
      real*8,intent(out) :: dseds(3)
      real*8,intent(out) :: d2seds2(3,3)
c
      integer csp,m,eta
      real*8 se1,wp,phiL,phiM,phiN,phiL_m,phiM_kp_m,phiN_m,se2k,
     1       jancae_bbc2008_get_se 
      real*8 wpi(2),dFds(3),dphiLds(3),dphiMds(3),dphiNds(3)
      real*8 d2Fds2(3,3),d2phiLds2(3,3),d2phiMds2(3,3),d2phiNds2(3,3)
      real*8 Lp(sp,3,3),Mp(sp,3,3),Np(sp,3,3)
      real*8 kCm(0:kp)
c-----------------------------------------------------------------------
c     se1: The factor to normalize equivalent stress calculated by 
c          original formula.
c     kCm(): Combination variables (caution, this array starts from zero.)
c          kCm(0)  = 2kp_C_0
c          kCm(1)  = 2kp_C_2
c          ...
c          kCm(kp) = 2kp_C_2kp
c
c     Variables expressed by "xp" in this code are used as "x" in the 
c       BBC2008 paper.
c   
c     If sp and kp don't have appropriate values, we can't setup Lp, Mp, 
c        Np tensors and kCm.
c   
c     pryld(1  ) = -6 (identification code for bbc2008)
c     pryld(1+1) = sp: s
c     pryld(1+2) = kp: k
c               wp: w = (1.5)^(1.0d0/s)
c   
c     pryld(3+1)~pryld(3+sp*8) are used for matrix Lp, Mp and Np 
c     i 1~sp
c     n0=1+2+(i-1)*8
c     pryld(n0+1) = l_1^(i) -> Lp(i,1~3,1~3)
c     pryld(n0+2) = l_2^(i) -> Lp(i,1~3,1~3)
c     pryld(n0+3) = m_1^(i) -> Mp(i,1~3,1~3)
c     pryld(n0+4) = m_2^(i) -> Mp(i,1~3,1~3)
c     pryld(n0+5) = m_3^(i) -> Mp(i,1~3,1~3)
c     pryld(n0+6) = n_1^(i) -> Np(i,1~3,1~3)
c     pryld(n0+7) = n_2^(i) -> Np(i,1~3,1~3)
c     pryld(n0+8) = n_3^(i) -> Np(i,1~3,1~3)   
c
c     Temporary variables.
c      wpi(1) : wp^(i-1)
c      wpi(2) : wp^(sp-i)
c
c      dFds(i) : dF/ds(i), see (x.y.2f)
c      d2Fds2(i,j) : d2F/(ds(i)ds(j)), see (x.y.2g)
c
c      phiX, (X=L,M,N) : see (x.y.1o)
c      dphiXds, (X=L,M,N): see (x.y.2d)
c      d2phiXds2, (X=L,M,N): see (x.y.2e)
c
c      phiL_m, phiN_m: phi_L^m, phi_N^m
c      phiM_kp_m : phi_M^(kp-m)
c
c     Here, (x.y.10), (x.y.2d) etc are equation numbering used in a 
c       text / a paper which will be released by JANCAE.
c
c-----------------------------------------------------------------------
c     ----------------
c        parameters
c     ----------------
c
      call jancae_bbc2008_setup ( pryld,ndyld,sp,kp,wp,
     1                            Lp,Mp,Np,kCm,se1)
c
c ----------------
c  se section
c ----------------
c
c                             ---- The unit of this se is (stress)^(2kp)
      se = jancae_bbc2008_get_se ( sp,kp,wp,s,Lp,Mp,Np,kCm)
c
c      ---- see eq.(x.y.2b) and eq.(x.y.2) for se2k and se, respectively
      se2k = se / se1
      se = se2k**(1.0d0 / (2.0d0*kp))
c
c   ---- If a main routine requests this subroutine to calculate se only
      if ( nreq == 0 ) then
        return
      end if
c
c
c --------------------------
c  dseds & d2seds2 section
c --------------------------
      call jancae_clear1 ( dFds,3 )
      call jancae_clear2 ( d2Fds2,3,3 )

c                              ---- long-long-long do loops starts here.
      do csp = 1,sp
c
        call jancae_bbc2008_get_w_phi ( wpi,phiL,phiM,phiN,
     1                                  csp,sp,wp,Lp,Mp,Np,s )
c
        do m = 0,kp
c
c     phiM^m, phiL^(kp-m) and phiN^m terms sometimes become 0**0.
c     To get consistency with the yield function and its differentials
c       disscussed by Banabic et al., 0**0 is required to be 1.
c
          phiL_m = 1.0d0
          phiN_m = 1.0d0
          if ( m /= 0 ) then
            phiL_m = phiL**m
            phiN_m = phiN**m
          end if
c
          phiM_kp_m = 1.0d0
          if ( (kp-m) /= 0 ) then
            phiM_kp_m = phiM**(kp-m)
          end if
c
c
          call jancae_bbc2008_get_dphiXds ( dphiLds,Lp,s,csp,m,sp )
          call jancae_bbc2008_get_dphiXds ( dphiMds,Mp,s,csp,(kp-m),sp )
          call jancae_bbc2008_get_dphiXds ( dphiNds,Np,s,csp,m,sp )
c
c                                             ---- <dseds>, see (x.y.2f)
          dFds(1:3) = dFds(1:3) + kCm(m) * 
     1                    (  wpi(1)*(  dphiMds(1:3) * phiL_m
     2                               + dphiLds(1:3) * phiM_kp_m )
     3                     + wpi(2)*(  dphiMds(1:3) * phiN_m
     4                               + dphiNds(1:3) * phiM_kp_m ))
c
c
c                     ---- <d2seds2>, see (x.y.2g), d2F/ds(eta)ds(gamma)
          if ( nreq ==2 ) then
c
            call jancae_bbc2008_get_d2phiXds2 (d2phiLds2,Lp,s,csp,m,sp)
            call jancae_bbc2008_get_d2phiXds2 (d2phiMds2,Mp,s,csp,
     1                                                       (kp-m),sp)
            call jancae_bbc2008_get_d2phiXds2 (d2phiNds2,Np,s,csp,m,sp)
c
            do eta = 1,3
              d2Fds2(eta,1:3) = d2Fds2(eta,1:3) + kCm(m)
     1             * (  wpi(1) * (  d2phiMds2(eta,1:3) * phiL_m
     2                            + dphiMds(1:3) * dphiLds(eta)
     3                            + dphiLds(1:3) * dphiMds(eta)
     4                            + d2phiLds2(eta,1:3) * phiM_kp_m )
     5                + wpi(2) * (  d2phiMds2(eta,1:3) * phiN_m
     6                            + dphiMds(1:3) * dphiNds(eta)
     7                            + dphiNds(1:3) * dphiMds(eta)
     8                            + d2phiNds2(eta,1:3) * phiM_kp_m ))
            end do
c
          end if
c
c                                                ---- end of m=0,kp loop
        end do
c
c                                                ---- end of i=1,sp loop
      end do
c
c
c                            ---- < dseds >, see (x.y.2f), se2k = se^2kp
      dseds(1:3) =  dFds(1:3) * se / (se1 * 2.0d0 * kp * se2k)

c                  ---- < d2seds2 >, see (x.y.2g), d2se/ds(eta)ds(gamma)
      if ( nreq == 2 ) then
        do eta = 1,3
          d2seds2(eta,1:3) = 
     1            d2Fds2(eta,1:3) * se / (se1 * 2.0d0 * kp * se2k)
     2            - (2.0d0*kp-1.0d0) * dseds(eta) * dseds(1:3) / se
        end do
      end if
c
      return
      end subroutine jancae_bbc2008_core
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     jancae_bbc2008_get_w_phi ()
c     A subroutine to get w^(i-1), w^(s-i) and phiX variables
c
      subroutine jancae_bbc2008_get_w_phi ( wpi,phiL,phiM,phiN,csp,sp,
     1                                      wp,Lp,Mp,Np,s )
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: csp, sp
      real*8 ,intent(in) :: wp
      real*8 ,intent(in) :: s(3)
      real*8 ,intent(in) :: Lp(sp,3,3),Mp(sp,3,3),Np(sp,3,3)
c
      real*8,intent(out) :: phiM, phiL, phiN
      real*8,intent(out) :: wpi(2)
c
      real*8 jancae_bbc2008_get_phiX
c-----------------------------------------------------------------------
c
      wpi(1) = wp**(csp-1)
      wpi(2) = wp**(sp-csp)
      phiL = jancae_bbc2008_get_phiX (Lp, s, csp , sp)
      phiM = jancae_bbc2008_get_phiX (Mp, s, csp , sp)
      phiN = jancae_bbc2008_get_phiX (Np, s, csp , sp)
c
      return
      end subroutine jancae_bbc2008_get_w_phi
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CALCULATE EQUIVALENT STRESS
c
      real*8 function jancae_bbc2008_get_se ( sp,kp,wp,s,Lp,Mp,Np,kCm )
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: sp,kp
      real*8 ,intent(in) :: wp
      real*8 ,intent(in) :: s(3)
      real*8 ,intent(in) :: Lp(sp,3,3),Mp(sp,3,3),Np(sp,3,3)
      real*8 ,intent(in) :: kCm(0:kp)
c
      integer csp,m
      real*8 phiM,phiL,phiN,phiL_m,phiM_kp_m,phiN_m
      real*8 wpi(2)
c-----------------------------------------------------------------------
c
      jancae_bbc2008_get_se = 0.0d0
c
      do csp = 1,sp
c
        call jancae_bbc2008_get_w_phi ( wpi,phiL,phiM,phiN,csp,sp,wp,
     1                                  Lp,Mp,Np,s )
c
        do m = 0,kp
c
          phiL_m = 1.0d0
          phiN_m = 1.0d0
          if ( m /= 0 ) then
            phiL_m = phiL**m
            phiN_m = phiN**m
          end if
c
          phiM_kp_m = 1.0d0
          if ( (kp-m) /= 0 ) then
            phiM_kp_m = phiM**(kp-m)
          end if
c
          jancae_bbc2008_get_se = 
     1    jancae_bbc2008_get_se 
     2     + kCm(m) * phiM_kp_m * ( wpi(1) * phiL_m + wpi(2) * phiN_m )
c
        end do
c
      end do
c
      return
      end function jancae_bbc2008_get_se
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     jancae_bbc2008_get_phiX (Xp, s, csp , sp)
c     A function to calculate s(a)*X(a,b)*s(b) (summation convention)
c
      real*8 function jancae_bbc2008_get_phiX ( Xp,s,csp,sp )
c-----------------------------------------------------------------------
      implicit none
c
      integer,parameter :: nc = 3
c
      integer,intent(in) :: csp,sp
      real*8 ,intent(in) :: s(nc)
      real*8 ,intent(in) :: Xp(sp,nc,nc)
c
      integer i
      real*8 v(nc)
      real*8 XXp(nc,nc)
c-----------------------------------------------------------------------
c     Xp: coefficient tensor
c     s : stress vector
c     csp: CURRENT sp value, 1<=csp<=sp
c
c     local variables
c     nc: the number of components.
c     XXp: = Xp(csp, nc, nc)
c-----------------------------------------------------------------------
c
c                                  ---- convert 3rd tensor to 2nd tensor
      do i = 1,nc
        XXp(i,1:nc) = Xp(csp,i,1:nc)
      end do
c
      call jancae_mv ( v,XXp,s,nc,nc)
      call jancae_vvs (jancae_bbc2008_get_phiX,v,s,nc)
c
      return
      end function jancae_bbc2008_get_phiX
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     jancae_bbc2008_get_dphiXds (dphiXds, Xp, s, csp, lambda,sp)
c     A subroutine to calculate d(phiX^(lambda))/ds.
c     It returns dphiXds(nc).
c
      subroutine jancae_bbc2008_get_dphiXds (dphiXds,Xp,s,csp,lambda,sp)
c-----------------------------------------------------------------------
      implicit none
c
      integer,parameter :: nc = 3
c
      integer,intent(in) :: csp,lambda,sp
      real*8 ,intent(in) :: s(nc)
      real*8 ,intent(in) :: Xp(sp,nc,nc)
c
      real*8,intent(out) :: dphiXds(nc)
c
      integer i
      real*8 phi
      real*8 v(nc)
      real*8 XXp(nc,nc)
c-----------------------------------------------------------------------
c     dphiXds: d(phiX^(lambda))/ds
c     Xp: coefficient tensor
c     s : stress vector
c     csp: CURRENT sp value, 1<=csp<=sp. 
c         Usually, it equals do-loop iterator "i".
c     lambda: see above.
c
c     local variables
c     nc: the number of components.
c     XXp: = Xp(csp, nc, nc)
c-----------------------------------------------------------------------
      call jancae_clear1(dphiXds,nc)
c
c                              ---- If lambda is 0, return dphiXds = {0}
      if ( lambda == 0) then
        return
      end if
c
c                                  ---- convert 3rd tensor to 2nd tensor
      do i = 1,nc
        XXp(i,1:nc) = Xp(csp, i, 1:nc)
      end do
c
c     In the bbc2008 section of the document "User subroutines for 
c       Metalic Plasticity model?", expression (x.y.2d) has 
c       "chi(gamma, beta)*s(beta)" (summation convention). 
c     In this routine, corresponding term is v(nc) obtained 
c       from jancae_mv().
c
      call jancae_mv (v, XXp, s, nc, nc)
c
      if ( lambda == 1 ) then
        dphiXds(1:nc) = 2.0d0 * v(1:nc)
      else
        call jancae_vvs (phi, v, s, nc)
        dphiXds(1:nc) = 2.0d0 * lambda * phi**(lambda-1) * v(1:nc)
      end if
c
      return
      end subroutine jancae_bbc2008_get_dphiXds
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     jancae_bbc2008_get_d2phiXds2 (d2phiXds2, Xp, s, csp, lambda,sp)
c     A subroutine to calculate d2(phiX^(lambda))/(dsds').
c     It returns d2phiXdsds(nc,nc).
c
      subroutine jancae_bbc2008_get_d2phiXds2 ( d2phiXds2,Xp,s,csp,
     1                                          lambda,sp )
c-----------------------------------------------------------------------
      implicit none
c
      integer,parameter :: nc = 3
c
      integer,intent(in) :: csp,lambda,sp
      real*8 ,intent(in) :: s(nc)
      real*8 ,intent(in) :: Xp(sp,nc,nc)
c
      real*8,intent(out) :: d2phiXds2(nc,nc)
c
      integer i
      real*8 phi,phi_lambda2
      real*8 v(nc)
      real*8 XXp(nc,nc)
c-----------------------------------------------------------------------
c     d2phiXdsds: d2(phiX^(lambda))/(dsds')
c     Xp: coefficient tensor
c     s : stress vector
c     csp: CURRENT sp value, 1<=csp<=sp. 
c         Usually, it equals do-loop iterator "i".
c     lambda: see above.
c
c     local variables
c     nc: the number of components.
c     XXp: = Xp(csp, nc, nc)
c-----------------------------------------------------------------------
c
c                             ---- see eq.(x.y.2e), the case lambda <= 1
      if ( lambda <= 1 ) then
        do i = 1,nc
          d2phiXds2(i,1:nc) = 2.0d0 * lambda * Xp(csp, i, 1:nc)
        end do
        return
      end if
c
c
      do i = 1,nc
        XXp(i,1:nc) = Xp(csp, i, 1:nc)
      end do
c
      call jancae_mv (v, XXp, s, nc, nc)
      call jancae_vvs (phi, v, s, nc)
c
      phi_lambda2 = 1.0d0
      if ( lambda /= 2 ) then
        phi_lambda2 = phi**(lambda-2)
      end if
c
      call jancae_clear2 ( d2phiXds2,nc,nc )
c
c                                            ---- d2phiX/(ds(i)ds(1:nc))
      do i = 1,nc
        d2phiXds2(i, 1:nc) = 2.0d0 * lambda * phi_lambda2 * 
     1   ( 2.0d0 * (lambda - 1) * v(1:nc) * v(i)
     2    + phi * XXp(1:nc, i) )
      end do
c
      return
      end subroutine jancae_bbc2008_get_d2phiXds2
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     setup_bbc2008_parameters()
c     A routine to setup local variables.
c
      subroutine jancae_bbc2008_setup ( pryld,ndyld,sp,kp,wp,
     1                                  Lp,Mp,Np,kCm,se1 )
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: sp,kp,ndyld
      real*8 ,intent(in) :: pryld(ndyld)
c
      real*8,intent(inout) :: wp,se1
      real*8,intent(inout) :: Lp(sp,3,3),Mp(sp,3,3),Np(sp,3,3)
      real*8,intent(inout) :: kCm(0:kp)
c
      integer csp,k,l,m,n
      real*8 jancae_bbc2008_get_se
      real*8 dummy_s(3)
      real*8 Comb(kp*2,0:kp*2)
c-----------------------------------------------------------------------
c
      wp = 1.5d0 ** (1.0d0/sp)
c
c                                        ---- Combination variables, kCm
      kCm(0) = 1.0d0
      kCm(kp) = 1.0d0
c
      do k = 1,2*kp
c                                    ---- caution: Comb has zero origin.
        Comb(k,0) = 1.0d0
        Comb(k,k) = 1.0d0
        do m = 1,k-1
          Comb(k, m) = Comb(k-1, m-1) + Comb(k-1, m)
        end do
c
c     We need Comb(k=2kp,2m) array.
c       Comb(k,0) --> kCm(0)  (which has already been setup)
c       Comb(k,k) --> kCm(kp) (which has also been setup)
c
c     In do m=2,k-2,2 loop, we have
c       Comb(k,2) --> kCm(1)
c       Comb(k,4) --> kCm(2)
c       ...
c       Comb(k,k-2) --> kCm(kp-1)
c
        if ( k == (2*kp) ) then
          n = 1
          do m = 2,k-2,2
            kCm(n) = Comb(k, m)
            n = n + 1
          end do
        end if

      end do
c
c                                                           ---- tensors
      do csp = 1,sp
c                                                      ---- L^(i) tensor
        l = 3 + 8 * (csp - 1)
        Lp(csp,1,1) = pryld(l+1)**2
        Lp(csp,1,2) = pryld(l+1)*pryld(l+2)
        Lp(csp,1,3) = 0.0d0
        Lp(csp,2,2) = pryld(l+2)**2
        Lp(csp,2,3) = 0.0d0
        Lp(csp,3,3) = 0.0d0
        Lp(csp,2,1) = Lp(csp,1,2)
        Lp(csp,3,1) = Lp(csp,1,3)
        Lp(csp,3,2) = Lp(csp,2,3)
c                                                      ---- M^(i) tensor
        m = l + 2
        call jancae_bbc2008_setup_MN_tensors (m,csp,pryld,ndyld,Mp,sp)
c                                                      ---- N^(i) tensor
        n = m + 3
        call jancae_bbc2008_setup_MN_tensors (n,csp,pryld,ndyld,Np,sp)
      end do
c
c
c     equiv. stress in uniaxial stress state.
c     dummy_s = (1.0d0, 0.0d0, 0.0d0)
c     ** The unit of this se1 is (stress)^(2kp)
c
      call jancae_clear1 (dummy_s, 3)
      dummy_s(1) = 1.0d0
      se1 = jancae_bbc2008_get_se (sp, kp, wp, dummy_s, Lp, Mp, Np, kCm)
c
      return
      end subroutine jancae_bbc2008_setup
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     setup_MN_tensors (ic,csp,pryld,ndyld, Xp,sp)
c       ic: initial counter
c       
c       This routine returns Mp or Np tensor.
c       Mp and Np tensors are the same style,
c       thus this subroutine has been created.
c
      subroutine jancae_bbc2008_setup_MN_tensors ( ic,csp,
     1                                             pryld,ndyld,Xp,sp )
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: ic,csp,sp,ndyld
      real*8 ,intent(in) :: pryld(ndyld)
c
      real*8,intent(inout) :: Xp(sp,3,3)
c-----------------------------------------------------------------------
c
      Xp(csp,1,1) = pryld(ic+1)**2
      Xp(csp,1,2) = -pryld(ic+1)*pryld(ic+2)
      Xp(csp,1,3) = 0.0d0
      Xp(csp,2,2) = pryld(ic+2)**2
      Xp(csp,2,3) = 0.0d0
      Xp(csp,3,3) = 4.0d0*pryld(ic+3)**2
      Xp(csp,2,1) = Xp(csp,1,2)
      Xp(csp,3,1) = Xp(csp,1,3)
      Xp(csp,3,2) = Xp(csp,2,3)
c
      return
      end subroutine jancae_bbc2008_setup_MN_tensors
c
c
c
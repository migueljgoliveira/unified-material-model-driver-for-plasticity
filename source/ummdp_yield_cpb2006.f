************************************************************************
c     CPB2006 YIELD FUNCTION AND DERIVATIVES
c
c       doi: 
c
c     !!! CAUTION !!!
c     Plane stress condition is NOT implemented in this code.
c
      subroutine jancae_cpb2006 ( s,se,dseds,d2seds2,nreq,
     1                            pryld,ndyld )
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
      integer i,j,k,m,n,l,iq,ip,ir
      real*8 pi,eps,a,ck,ai,H1,H2,H3,p,q,theta,F,D,DseDF,denom,D2seDF2,
     1       del,sea,seb,abc1,abc2,seaa,seba,seab,sebb,cpb2006_seND
      real*8 s0(6),sigma(6),psigma(3),phi(3),psi(3),omega(3),DFDH(3),
     1       DFDpsigma(3),DFDs(6)
      real*8 c(6,6),ct(6,6),DpsigmaDH(3,3),DHdsigma(3,6),DsigmaDs(6,6),
     1       D2FDpsigma2(3,3),D2FDH2(3,3),D2FDs2(6,6),dummat(3,6)
      real*8 D2psigmaDH2(3,3,3),D2HDsigma2(3,6,6)     
c-----------------------------------------------------------------------
c                                                         ---- variables
c
c     sigma(6)  : Linear-transformed stress
c     psigma(3) : Principal values of sigma
c     ct(6,6)   : Transformation matrix from Cauchy stress to sigma(6)
c
c                                         --- For 1st order differential
c
c     dseds(6)       : D(se)/D(s) -> 1st differential
c     DFDH(3)        : D(F)/D(H)
c     DFDpsigma(3)   : D(F)/D(psigma)
c     DpsigmaDH(3,3) : D(psigma)/D(H)
c     DHDsigma(3,6)  : D(H)/D(sigma)
c     DsigmaDs(6,6)  : D(sigma)/D(s)
c     DFDs(6)        : D(F)/D(s)
c
c                                        ---- For 2nd order differential
c
c     d2seds2(6,6)       : D2(se)/D(s)2 -> 2nd differential
c     D2FDpsigma(3,3)    : D2(F)/D(psigma)2
c     D2psigmaDH2(3,3,3) : D2(psigma)/D(H)2
c     D2FDH2(3,3)        : D2(F)/D(H)2
c     D2HDsigma2(3,6,6)  : D2(H)/D(sigma)2
c     D2FDs2(6,6)        : D2(F)/D(s)2
c-----------------------------------------------------------------------
c
      pi = acos(-1.0d0)
      eps = 1.0d-5
c                                        ---- set anisotropic parameters
c
      call jancae_clear2 ( c,6,6 )
c
      c(1,1) = pryld(1+1)            ! C11
      c(1,2) = pryld(1+2)            ! C12
      c(1,3) = pryld(1+3)            ! C13
      c(2,1) = pryld(1+4)            ! C21
      c(2,2) = pryld(1+5)            ! C22
      c(2,3) = pryld(1+6)            ! C23
      c(3,1) = pryld(1+7)            ! C31
      c(3,2) = pryld(1+8)            ! C32
      c(3,3) = pryld(1+9)            ! C33
      c(4,4) = pryld(1+10)           ! C44     ! tau_xy
      c(5,5) = pryld(1+11)           ! C55     ! tau_xz
      c(6,6) = pryld(1+12)           ! C66     ! tau_yz
      a      = pryld(1+13)           ! a
      ck     = pryld(1+14)           ! k
c
      ai = 1.0d0 / a
c
c                                         ---- Calculate phi, psi, omega
      phi(1) = (2.0d0*c(1,1) - c(1,2) - c(1,3)) / 3.0d0
      phi(2) = (2.0d0*c(1,2) - c(2,2) - c(2,3)) / 3.0d0
      phi(3) = (2.0d0*c(1,3) - c(2,3) - c(3,3)) / 3.0d0
c
      psi(1) = (-c(1,1) + 2.0d0*c(1,2) - c(1,3)) / 3.0d0
      psi(2) = (-c(1,2) + 2.0d0*c(2,2) - c(2,3)) / 3.0d0
      psi(3) = (-c(1,3) + 2.0d0*c(2,3) - c(3,3)) / 3.0d0
c
      omega(1) = (c(1,1) + c(1,2) - 2.0d0*c(1,3)) / 3.0d0
      omega(2) = (c(1,2) + c(2,2) - 2.0d0*c(2,3)) / 3.0d0
      omega(3) = (c(1,3) + c(2,3) - 2.0d0*c(3,3)) / 3.0d0
c
c      ---- Calculate 4th order orthotropic tensor "L" ( named ct here )
      call jancae_clear2 ( ct,6,6 )
c
      ct(1,1) = phi(1)
      ct(1,2) = psi(1)
      ct(1,3) = -omega(1)
      ct(2,1) = phi(2)
      ct(2,2) = psi(2)
      ct(2,3) = -omega(2)
      ct(3,1) = phi(3)
      ct(3,2) = psi(3)
      ct(3,3) = -omega(3)
      ct(4,4) = c(4,4)
      ct(5,5) = c(5,5)
      ct(6,6) = c(6,6)
c
c                               ---- Calculate linear transformed stress
      call jancae_mv ( sigma,ct,s,6,6 )
c
c            ---- Calculate principal values of transformed stress sigma
c                                                     by Cardan's method
c
c                                       ---- 1st, 2nd and 3rd invariants
      H1 = (sigma(1) + sigma(2) + sigma(3)) / 3.0d0
      H2 = (sigma(5)**2.0d0 + sigma(6)**2.0d0 + sigma(4)**2.0d0 -
     1      sigma(2)*sigma(3) - sigma(3)*sigma(1) - 
     2      sigma(1)*sigma(2)) / 3.0d0
      H3 = (2.0d0*sigma(5)*sigma(6)*sigma(4) + 
     1            sigma(1)*sigma(2)*sigma(3) - 
     2            sigma(1)*sigma(6)**2.0d0 - 
     3            sigma(2)*sigma(5)**2.0d0 -
     4      sigma(3)*sigma(4)**2.0d0) / 2.0d0 ! sigma(5) <-> sigma(6)
c
      p = H1**2.0d0 + H2
      q = (2.0d0*H1**3.0d0 + 3.0d0*H1*H2 + 2.0d0*H3) / 2.0d0
      if ( abs(p) >= 1.0d-16 ) then
        theta = q / (p**1.5d0)
        if ( theta > 1.0d0 ) theta = 1.0d0
        if ( theta < -1.0d0 ) theta =-1.0d0
        theta = acos(theta)
      else
        theta = 0.0
      end if
c
c                               ---- calculate principal values of sigma
      psigma(1) = 2.0d0*sqrt(p)*cos(theta/3.0d0) + H1
      psigma(2) = 2.0d0*sqrt(p)*cos((theta+4.0d0*pi)/3.0d0) + H1
      psigma(3) = 2.0d0*sqrt(p)*cos((theta+2.0d0*pi)/3.0d0) + H1
c
c                                                 ---- equivalent stress
c
c                                          ---- calculate yield function
      F = (abs(psigma(1)) - ck*psigma(1))**a + 
     1    (abs(psigma(2)) - ck*psigma(2))**a + 
     2    (abs(psigma(3)) - ck*psigma(3))**a
c                                           ---- denominator coefficient
      D = (abs(phi(1)) - ck*phi(1))**a + (abs(phi(2)) - ck*phi(2))**a +
     1    (abs(phi(3)) - ck*phi(3))**a
c
      se = (F/D) ** ai
c
c                                            ---- 1st order differential
      if ( nreq >= 1 ) then
c                                              ---- D(se)/D(F) -> Scalar
        DseDF = (1.0d0/D)**ai * ai * F**(ai-1.0d0)
c                                    ---- D(F)/D(psigma) -> 1 x 3 Vector
        do i = 1,3
          DFDpsigma(i) = a * (psigma(i)/abs(psigma(i))-ck) *
     1                  (abs(psigma(i))-ck*psigma(i))**(a-1.0d0)
        end do
c
c                                 ---- D(F)/D(H) by using D(psigma)/D(H)
c                                         D(F)/D(H)      -> 1 x 3 Vector
c                                         D(psigma)/D(H) -> 3 x 3 Matrix
        call jancae_clear1 ( DFDH,3 )
        call jancae_clear2 ( DpsigmaDH,3,3 )
c
        if ( abs(psigma(2)-psigma(3)) / se > eps .and.
     1       abs(psigma(2)-psigma(1)) / se > eps ) then
c                                                 ---- not Singular case
          do i = 1,3
            denom = psigma(i)**2.0d0 - 2.0d0*H1*psigma(i) - H2
            DpsigmaDH(i,1) = psigma(i)**2.0d0/denom
            DpsigmaDH(i,2) = psigma(i)/denom
            DpsigmaDH(i,3) = 2.0d0/3.0d0/denom
          end do
c
          do i = 1,3
            do j = 1,3
              DFDH(i) = DFDH(i) + DFDpsigma(j)*DpsigmaDH(j,i)
            end do
          end do
        else
c                                                     ---- singular case
          if ( abs(psigma(2)-psigma(3)) / se <= eps) then
c                                           ---- Case1 S2=S3 ( theta=0 )
            denom = psigma(1)**2.0d0 - 2.0d0*H1*psigma(1) - H2
            DpsigmaDH(1,1) = psigma(1)**2.0d0/denom
            DpsigmaDH(1,2) = psigma(1) / denom
            DpsigmaDH(1,3) = 2.0d0/3.0d0/denom
c
            DFDH(1) = DpsigmaDH(1,1)*(DFDpsigma(1)-DFDpsigma(2)) +
     1                3.0d0*DFDpsigma(2)
            DFDH(2) = DpsigmaDH(1,2) * (DFDpsigma(1)-DFDpsigma(2))
            DFDH(3) = DpsigmaDH(1,3) * (DFDpsigma(1)-DFDpsigma(2))
          else if ( abs(psigma(2)-psigma(1)) / se <= eps ) then
c
c                                          ---- Case2 S2=S1 ( theta=pi )
            denom = psigma(3)**2.0d0 - 2.0d0*H1*psigma(3) - H2
            DpsigmaDH(3,1) = psigma(3)**2.0d0/denom
            DpsigmaDH(3,2) = psigma(3) / denom
            DpsigmaDH(3,3) = 2.0d0/3.0d0/denom
c
            DFDH(1) = DpsigmaDH(3,1)*(DFDpsigma(3)-DFDpsigma(2)) +
     1                3.0d0*DFDpsigma(2)
            DFDH(2) = DpsigmaDH(3,2) * (DFDpsigma(3)-DFDpsigma(2))
            DFDH(3) = DpsigmaDH(3,3) * (DFDpsigma(3)-DFDpsigma(2))
          end if
        end if
c
c                                     ---- D(H)/D(sigma) -> 3 x 6 Matrix
        call jancae_clear2 ( DHDsigma,3,6 )
c
        DHDsigma(1,1) = 1.0d0 / 3.0d0
        DHDsigma(1,2) = 1.0d0 / 3.0d0
        DHDsigma(1,3) = 1.0d0 / 3.0d0
c
        DHDsigma(2,1) = -1.0d0/3.0d0*(sigma(2)+sigma(3))
        DHDsigma(2,2) = -1.0d0/3.0d0*(sigma(3)+sigma(1))
        DHDsigma(2,3) = -1.0d0/3.0d0*(sigma(1)+sigma(2))
        DHDsigma(2,4) = 2.0d0/3.0d0*sigma(4)
        DHDsigma(2,5) = 2.0d0/3.0d0*sigma(5)
        DHDsigma(2,6) = 2.0d0/3.0d0*sigma(6)
c
c                                            !!-sigma(5)**2.0d0)
        DHDsigma(3,1) = 0.5d0 * (sigma(2)*sigma(3)-sigma(6)**2.0d0)
c                                            !!-sigma(6)**2.0d0)
        DHDsigma(3,2) = 0.5d0 * (sigma(3)*sigma(1)-sigma(5)**2.0d0)
        DHDsigma(3,3) = 0.5d0 * (sigma(1)*sigma(2)-sigma(4)**2.0d0)
        DHDsigma(3,4) = sigma(5)*sigma(6) - sigma(3)*sigma(4)
        DHDsigma(3,6) = sigma(6)*sigma(4) - sigma(1)*sigma(5) !!...(3,5)=
        DHDsigma(3,5) = sigma(4)*sigma(5) - sigma(2)*sigma(6) !!...(3,6)=
c
c                                     ---- D(sigma)/D(s) -> 6 x 6 Matrix
        do i = 1,6
          do j = 1,6
            DsigmaDs(i,j) = ct(i,j)
          end do
        end do
c
c                                        ---- D(se)/D(s) -> 1 x 3 Vector
        call jancae_clear1 ( DFDs,6 )
        call jancae_clear2 ( dummat,3,6 )
c
        do i = 1,6
          do j = 1,3
            do k = 1,6
              dummat(j,i) = dummat(j,i) + DHDsigma(j,k)*DsigmaDs(k,i)
            end do
            DFDs(i) = DFDs(i) + DFDH(j)*dummat(j,i)
          end do
        end do
      end if
c
      do i = 1,6
        dseds(i) = DseDF*DFDs(i)
      end do
c
c                                            ---- 2nd order differential
      if ( nreq >= 2 ) then
c                                              ---- D(se)/D(F) -> Scalar
        D2seDF2 = (1.0d0/D)**ai * ai * (ai-1.0d0) * F**(ai-2.0d0)
c
c                                  ---- D2(F)/D(psigma)2 -> 3 x 3 Matrix
        call jancae_clear2 ( D2FDpsigma2,3,3 )
c
        do i = 1,3
          D2FDpsigma2(i,i) = a*(psigma(i)/abs(psigma(i))-ck)**2.0d0 *
     1                       (abs(psigma(i))-ck*psigma(i))**(a-2.0d0)
        end do
c
c                              ---- D2(psigma)/D(H)2 -> 3 x 3 x 3 Matrix
        call jancae_clear3 ( D2psigmaDH2,3,3,3 )
c
        if ( abs(psigma(2)-psigma(3)) / se > eps .and.
     1       abs(psigma(2)-psigma(1)) / se > eps ) then
c                                                 ---- Not Singular case
c
          do i = 1,3
            denom = (psigma(i)**2.0d0-2.0d0*H1*psigma(i)-H2) ** 3.0d0
            D2psigmaDH2(i,1,1) = 2.0d0 * psigma(i)**3.0d0 * 
     1                         (psigma(i)**2.0d0-3.0d0*H1*psigma(i)
     2                         -2.0d0*H2) / denom
            D2psigmaDH2(i,2,2) = -2.0d0 * psigma(i) * 
     1                           (H1*psigma(i)+H2) / denom
            D2psigmaDH2(i,3,3) = -8.0d0/9.0d0 * (psigma(i)-H1) / denom
            D2psigmaDH2(i,1,2) = psigma(i)**2.0d0 * (psigma(i)**2.0d0 -
     1                           4.0d0*H1*psigma(i)-3.0d0*H2) / denom
            D2psigmaDH2(i,2,3) = (-2.0d0/3.0d0) * 
     1                           (psigma(i)**2.0d0+H2) / denom
            D2psigmaDH2(i,3,1) = -4.0d0/3.0d0 * psigma(i) * 
     1                           (H1*psigma(i)+H2) / denom
c
            D2psigmaDH2(i,2,1) = D2psigmaDH2(i,1,2)
            D2psigmaDH2(i,3,2) = D2psigmaDH2(i,2,3)
            D2psigmaDH2(i,1,3) = D2psigmaDH2(i,3,1)
          end do
c
c                                       ---- D2(F)/D(H)2 -> 3 x 3 Matrix
          call jancae_clear2 ( D2FDH2,3,3 )
c
          do iq = 1,3
            do m = 1,3
              do ip = 1,3
                do l = 1,3
                  D2FDH2(iq,m) = D2FDH2(iq,m) + D2FDpsigma2(ip,l)*
     1                           DpsigmaDH(l,m)*DpsigmaDH(ip,iq)
                end do
                D2FDH2(iq,m) = D2FDH2(iq,m) + DFDpsigma(ip)*
     1                         D2psigmaDH2(ip,iq,m)
              end do
            end do
          end do
c
c                               ---- D2(H)/D(sigma)2 -> 3 x 6 x 6 Matrix
          call jancae_clear3 ( D2HDsigma2,3,6,6 )
c
          D2HDsigma2(2,1,2) = -1.0d0 / 3.0d0
          D2HDsigma2(2,2,3) = -1.0d0 / 3.0d0
          D2HDsigma2(2,3,1) = -1.0d0 / 3.0d0
c
          D2HDsigma2(2,2,1) = D2HDsigma2(2,1,2)
          D2HDsigma2(2,3,2) = D2HDsigma2(2,2,3)
          D2HDsigma2(2,1,3) = D2HDsigma2(2,3,1)
c
          D2HDsigma2(2,4,4) = 2.0d0 / 3.0d0
          D2HDsigma2(2,5,5) = 2.0d0 / 3.0d0
          D2HDsigma2(2,6,6) = 2.0d0 / 3.0d0
c
          D2HDsigma2(3,1,2) = sigma(3) / 2.0d0
          D2HDsigma2(3,2,3) = sigma(1) / 2.0d0
          D2HDsigma2(3,3,1) = sigma(2) / 2.0d0
c
          D2HDsigma2(3,2,1) = D2HDsigma2(3,1,2)
          D2HDsigma2(3,3,2) = D2HDsigma2(3,2,3)
          D2HDsigma2(3,1,3) = D2HDsigma2(3,3,1)
c
          D2HDsigma2(3,4,4) = -sigma(3)
          D2HDsigma2(3,6,6) = -sigma(1)         !!...(3,5,5)
          D2HDsigma2(3,5,5) = -sigma(2)         !!...(3,6,6)
c
          D2HDsigma2(3,4,5) = sigma(6)
          D2HDsigma2(3,5,6) = sigma(4)
          D2HDsigma2(3,6,4) = sigma(5)
c
          D2HDsigma2(3,5,4) = D2HDsigma2(3,4,5)
          D2HDsigma2(3,6,5) = D2HDsigma2(3,5,6)
          D2HDsigma2(3,4,6) = D2HDsigma2(3,6,4)
c
          D2HDsigma2(3,1,6) = -sigma(6)         !!...(3,1,5)=-sigma(5)
          D2HDsigma2(3,6,1) = D2HDsigma2(3,1,6) !!...(3,5,1)=...(3,1,5)
c
          D2HDsigma2(3,2,5) = -sigma(5)         !!...(3,2,6)=-sigma(6)
          D2HDsigma2(3,5,2) = D2HDsigma2(3,2,5) !!...(3,6,2)=...(3,2,6)
c
          D2HDsigma2(3,3,4) = -sigma(4)
          D2HDsigma2(3,4,3) = D2HDsigma2(3,3,4)
c
c                                       ---- D2(F)/D(s)2 -> 6 x 6 Matrix
          call jancae_clear2 ( D2FDs2,6,6 )
          call jancae_clear2 ( dummat,3,6 )
c
          do i = 1,3
            do j = 1,6
              do ip = 1,6
                dummat(i,j) = dummat(i,j) + 
     1                        DHDsigma(i,ip)*DsigmaDs(ip,j)
              end do
            end do
          end do
c
          do i = 1,6
            do j = 1,6
              do iq = 1,3
                do m = 1,3
                  D2FDs2(i,j) = D2FDs2(i,j) +
     1                          D2FDH2(iq,m)*dummat(iq,i)*dummat(m,j)
                end do
c
                do n = 1,6
                  do ir = 1,6
                    D2FDs2(i,j) = D2FDs2(i,j) + D2HDsigma2(iq,ir,n) *
     1                            DFDH(iq)*DsigmaDs(ir,i)*DsigmaDs(n,j)
                  end do
                end do
              end do
            end do
          end do
c                                      ---- D2(se)/D(s)2 -> 6 x 6 Matrix
          do i = 1,6
            do j = 1,6
              d2seds2(i,j) = D2seDF2*DfDs(i)*DfDs(j) + DseDF*D2FDs2(i,j)
            end do
          end do
        else
c                                                     ---- Singular case
          del = eps
c
          do i = 1,6
            s0(i) = s(i)
          end do
c
          do i = 1,6
            do j = 1,6
              if ( i == j ) then
                s0(i) = s(i) - del
                sea = cpb2006_seND(s0,ct,phi,ck,a,ai)
                s0(i) = s(i) + del
                seb = cpb2006_seND(s0,ct,phi,ck,a,ai)
c
                s0(i) = s(i)
                abc1 = (se-sea) / del
                abc2 = (seb-se) / del
                d2seds2(i,j) = (abc2-abc1) / del
              else
                s0(i) = s(i) - del
                s0(j) = s(j) - del
                seaa = cpb2006_seND(s0,ct,phi,ck,a,ai)
c
                s0(i) = s(i) + del
                s0(j) = s(j) - del
                seba = cpb2006_seND(s0,ct,phi,ck,a,ai)
c
                s0(i) = s(i) - del
                s0(j) = s(j) + del
                seab = cpb2006_seND(s0,ct,phi,ck,a,ai)
c
                s0(i) = s(i) + del
                s0(j) = s(j) + del
                sebb = cpb2006_seND(s0,ct,phi,ck,a,ai)
c
                s0(i) = s(i)
                s0(j) = s(j)
                abc1 = (seba-seaa) / (2.0d0*del)
                abc2 = (sebb-seab) / (2.0d0*del)
                d2seds2(i,j) = (abc2-abc1) / (2.0d0*del)
              end if
            end do
          end do
        end if
      end if
c
      return
      end subroutine jancae_cpb2006
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     NUMERICAL DIFFERENTIATION FOR EQUIVALENT STRESS
c
      real*8 function cpb2006_seND ( s,ct,phi,ck,a,ai )
c-----------------------------------------------------------------------
      implicit none
c 
      real*8 ck,a,ai
      real*8 s(6),phi(3)
      real*8 ct(6,6)
c
      real*8 pi,H1,H2,H3,p,q,theta,F,D
      real*8 sigma(6),psigma(3)
c-----------------------------------------------------------------------
c
      pi = acos(-1.0d0)
c                               ---- calculate linear transformed stress
      call jancae_mv ( sigma,ct,s,6,6 )
c                                       ---- 1st, 2nd and 3rd invariants
      H1 = (sigma(1)+sigma(2)+sigma(3)) / 3.0d0
      H2 = (sigma(5)**2.0d0+sigma(6)**2.0d0+sigma(4)**2.0d0 - 
     1     sigma(2)*sigma(3)-sigma(3)*sigma(1)-sigma(1)*sigma(2))/3.0d0
      H3 = (2.0d0*sigma(5)*sigma(6)*sigma(4) + 
     1            sigma(1)*sigma(2)*sigma(3) -
     2            sigma(1)*sigma(6)**2.0d0 - 
     3            sigma(2)*sigma(5)**2.0d0 -
     4            sigma(3)*sigma(4)**2.0d0) / 2.0d0 ! sigma(5) <-> sigma(6)
c
      p = H1**2.0d0+H2
      q = (2.0d0*H1**3.0d0+3.0d0*H1*H2+2.0d0*H3) / 2.0d0
      theta = q/(p**1.5d0)
      if ( theta > 1.0d0 ) theta = 1.0d0
      if ( theta < -1.0d0 ) theta = -1.0d0
      theta = acos(theta)
c                               ---- calculate principal values of sigma
      psigma(1) = 2.0d0*sqrt(p)*cos(theta/3.0d0) + H1
      psigma(2) = 2.0d0*sqrt(p)*cos((theta+4.0d0*pi)/3.0d0) + H1
      psigma(3) = 2.0d0*sqrt(p)*cos((theta+2.0d0*pi)/3.0d0) + H1
c
c                                                 ---- equivalent stress
c                                          ---- calculate yield function
      F = (abs(psigma(1))-ck*psigma(1))**a +
     &    (abs(psigma(2))-ck*psigma(2))**a +
     &    (abs(psigma(3))-ck*psigma(3))**a
c                                           ---- denominator coefficient
      D = (abs(phi(1))-ck*phi(1))**a +
     &    (abs(phi(2))-ck*phi(2))**a +
     &    (abs(phi(3))-ck*phi(3))**a
c
      cpb2006_seND = (F/D) ** ai
c
      return
      end function cpb2006_seND
c
c
c

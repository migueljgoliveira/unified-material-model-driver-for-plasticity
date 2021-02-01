c------------------------------------------------------(karafillisboyce)
c     Karafillis-Boyce Yield Function
c
      subroutine jancae_KarafillisBoyce ( s,se,dseds,d2seds2,nreq,
     &                                    pryld,ndyld )
c-----------------------------------------------------------------------
      implicit none
c
      real*8, parameter :: TOL=1e-7
c                                                    ---- input & output
      integer, intent(in) :: nreq, ndyld
      real*8, intent(in) :: s(6), pryld(ndyld)
      real*8, intent(out) :: se, dseds(6), d2seds2(6,6)
c                                                   ---- local variables
      real*8 c_p, L(6,6), smallS(6), largeS(3), Jinvar(3), phi, phiN(2),
     &   DseDphi, DjDss(3,6), DphiDs(2,6), DphiDj(2,3),
     &   DphiDls(2,3), DlsDj(3,3),
     &   DDphiDDs(2,6,6), DDjDDss(3,6,6), DDphiDDj(2,3,3),
     &   DDphiDDls(2,3,3), X12, X13, X23, DDlsDDj(3,3,3),
     &   dum, workmat(3,6), coef(2), workmat1(2,6,6)
      integer k_p
      integer i, j, k, m, n, eqFlag
      real*8  beta(3),alpha(2)
c
c     local variables  : symbols in Karafillis-Boyce Yield Function's 
c                        document
c
c     c_p       : c parameter
c     k_p       : k parameter
c     L(6,6)    : linear transformation matrix
c     smallS(i) : IPE deviatoric Cauchy stress tensor (s_tilde_i, i=1~6)
c     largeS(i) : principal values of smallS (S_tilde_i, i=1,2,3)
c     Jinvar(i) : invariants (J_i, i=1,2,3)
c     phi       : Phi
c     phiN(i)   : Phi1 and Phi2
    
c     DseDphi          : D(se)/D(phi)
c     DjDss(i,j)       : D(J_i)/D(s_tilde_j)               (i=1~3, j=1~6)
c     DphiDs(i,j)      : D(Phi_i)/D(sigma_j)               (i=1~2, j=1~6)
c     DphiDj(i,j)      : D(Phi_i)/D(J_j)                   (i=1~2, j=1~3)
c     DphiDls(i,j)     : D(Phi_i)/D(S_tilde_j)             (i=1~2, j=1~3)
c     DlsDj(i,j)       : D(S_tilde_i)/D(J_j)               (i=1~3, j=1~3)
c     DDphiDDs(i,j,k)  : D2(Phi_i)/D(sigma_j)D(sigma_k)     
c                                                  (i=1~2, j=1~6, k=1~6)
c
c     DDjDDss(i,j,k)   : D2(J_i)/D(s_tilde_j)D(s_tilde_k)   
c                                                  (i=1~3, j=1~6, k=1~6)
c
c     DDphiDDj(i,j,k)  : D2(Phi_i)/D(J_j)D(J_k)             
c                                                  (i=1~2, j=1~3, k=1^3)
c
c     DDphiDDls(i,j,k) : D2(Phi_i)/D(S_tilde_j)D(S_tilde_k) 
c                                                  (i=1~2, j=1~3, k=1^3)
c
c     DDlsDDj(i,j,k)   : D2(S_tilde_i)/D(J_j)D(J_k) 
c                                                  (i=1~3, j=1~3, k=1^3)
c
c     X12      : X12
c     X13      : X13
c     X23      : X23
c     coef(2)  : coefficient at "Phi = coef(1) * Phi1 + coef(2) * Phi2"
c      
c                                                        ---- parameters
      do i=1,6
        do j=1,6
          L(i,j) = 0.0d0
        end do
      end do
      alpha(1)=pryld(1+2)
      alpha(2)=pryld(1+3)
      beta(1)=(alpha(2)-1.0d0   -alpha(1))*0.5d0
      beta(2)=(alpha(1)-alpha(2)-1.0d0   )*0.5d0
      beta(3)=(1.0d0   -alpha(1)-alpha(2))*0.5d0
      L(1,1)=pryld(1+1)* 1.0d0
      L(1,2)=pryld(1+1)* beta(1)
      L(1,3)=pryld(1+1)* beta(2)
      L(2,2)=pryld(1+1)* alpha(1)
      L(2,3)=pryld(1+1)* beta(3)
      L(3,3)=pryld(1+1)* alpha(2)
      L(4,4)=pryld(1+1)* pryld(1+4)
      L(5,5)=pryld(1+1)* pryld(1+5)
      L(6,6)=pryld(1+1)* pryld(1+6)
      L(2,1)=L(1,2)
      L(3,1)=L(1,3)
      L(3,2)=L(2,3)
      k_p = pryld(1+7)
      c_p = pryld(1+8)
c                                                 ---- equivalent stress
      call jancae_mv (smallS, L, s, 6, 6)
      call jancae_KarafillisBoyce_principalStress(smallS, Jinvar,largeS)
c
      phiN(1) = (largeS(1)-largeS(2))**(2*k_p)
     &        + (largeS(2)-largeS(3))**(2*k_p)
     &        + (largeS(3)-largeS(1))**(2*k_p)
      phiN(2) = largeS(1)**(2*k_p) 
     &        + largeS(2)**(2*k_p)
     &        + largeS(3)**(2*k_p)
      coef(1) = (1d0-c_p)
      coef(2) = c_p * 3d0**(2*k_p)/(2d0**(2*k_p-1)+1d0)
      phi = coef(1)*phiN(1) + coef(2)*phiN(2)
      se=(0.5d0 * phi)**(0.5d0/k_p)
c
c                                            ---- 1st order differential
      if ( nreq.ge.1 ) then
c                 ---- check if there are two components of largeS equal 
c                                                          to each other
c               ---- if so, rearrange largeS so that largeS(1)=largeS(2)
        eqFlag = 0
        if (abs(largeS(1)-largeS(2)).le.TOL) then
          eqFlag = 1
        else if (abs(largeS(2)-largeS(3)).le.TOL) then
          eqFlag = 2
          dum = largeS(3)
          largeS(3) = largeS(1)
          largeS(1) = dum
        else if (abs(largeS(1)-largeS(3)).le.TOL) then
          eqFlag = 3
          dum = largeS(3)
          largeS(3) = largeS(2)
          largeS(2) = dum
        end if

        if (eqFlag.eq.0) then
          DphiDls(1,1) = 2d0 * k_p * ((largeS(1)-largeS(2))**(2*k_p-1)
     &       + (largeS(1)-largeS(3))**(2*k_p-1))
          DphiDls(1,2) = 2d0 * k_p * ((largeS(2)-largeS(1))**(2*k_p-1)
     &       + (largeS(2)-largeS(3))**(2*k_p-1))
          DphiDls(1,3) = 2d0 * k_p * ((largeS(3)-largeS(1))**(2*k_p-1)
     &       + (largeS(3)-largeS(2))**(2*k_p-1))
          do i=1,3
            DphiDls(2,i) = 2d0 * k_p * largeS(i)**(2*k_p-1)
          end do
c 
          do i=1,3
            dum = 1d0 /(3d0*largeS(i)**2 - 2d0*Jinvar(1)*largeS(i)
     &         + Jinvar(2))
            DlsDj(i,1) = largeS(i)**2 * dum
            DlsDj(i,2) = -largeS(i) * dum
            DlsDj(i,3) = dum
          end do
c
          do i=1,2
            do j=1,3
              DphiDj(i,j) = 0d0
              do k=1,3
                DphiDj(i,j) = DphiDj(i,j) + DphiDls(i,k) * DlsDj(k,j)
              end do
            end do
          end do
        else
          if (k_p.eq.1) then
            DphiDj(1,1) = 4d0*(2d0*largeS(1)+largeS(3))
            DphiDj(1,2) = -6d0
            DphiDj(1,3) = 0d0
          else
            dum = (largeS(1)-largeS(3))**(2*k_p-3)
            DphiDj(1,1) = 4d0*k_p*dum
     &         *(k_p*largeS(1)**2-largeS(1)*largeS(3)-largeS(3)**2)
            DphiDj(1,2) = -2d0*k_p*dum
     &         *((2*k_p-1)*largeS(1)-3d0*largeS(3))
            DphiDj(1,3) = 4d0*k_p*(k_p-2)*dum
          end if
c
          if (abs(largeS(1)-largeS(3)).le.TOL) then
            DphiDj(2,1) = 2d0*k_p**2*(2*k_p+1)*largeS(1)**(2*k_p-1)
            DphiDj(2,2) = -2d0*k_p**2*(2*k_p-1)*largeS(1)**(2*k_p-2)
            DphiDj(2,3) = 2d0*k_p*(2*k_p-1)*(k_p-1)
     &         *largeS(1)**(2*k_p-3)
          else
            dum = 2d0*k_p/(largeS(1)-largeS(3))**2
            DphiDj(2,1) = dum * ((2*k_p+1)*largeS(1)**(2*k_p)*
     &         (largeS(1)-largeS(3))
     &         - largeS(1)**(2*k_p+1) + largeS(3)**(2*k_p+1))
            DphiDj(2,2) = -dum * (2*k_p*largeS(1)**(2*k_p-1)*
     &         (largeS(1)-largeS(3))
     &         - largeS(1)**(2*k_p) + largeS(3)**(2*k_p))
            DphiDj(2,3) = dum * ((2*k_p-1)*largeS(1)**(2*k_p-2)*
     &         (largeS(1)-largeS(3))
     &         - largeS(1)**(2*k_p-1) + largeS(3)**(2*k_p-1))
          end if
        end if
c
        DjDss(1,1) = 1d0
        DjDss(1,2) = 1d0
        DjDss(1,3) = 1d0
        DjDss(1,4) = 0d0
        DjDss(1,5) = 0d0
        DjDss(1,6) = 0d0
        DjDss(2,1) = smallS(2) + smallS(3)
        DjDss(2,2) = smallS(1) + smallS(3)
        DjDss(2,3) = smallS(1) + smallS(2)
        DjDss(2,4) = -2d0*smallS(4)
        DjDss(2,5) = -2d0*smallS(5)
        DjDss(2,6) = -2d0*smallS(6)
        DjDss(3,1) = smallS(2)*smallS(3) - smallS(5)**2
        DjDss(3,2) = smallS(1)*smallS(3) - smallS(6)**2
        DjDss(3,3) = smallS(1)*smallS(2) - smallS(4)**2
        DjDss(3,4) = 2d0*(smallS(5)*smallS(6) - smallS(3)*smallS(4))
        DjDss(3,5) = 2d0*(smallS(4)*smallS(6) - smallS(1)*smallS(5))
        DjDss(3,6) = 2d0*(smallS(4)*smallS(5) - smallS(2)*smallS(6))
c
        do i=1,3
          do j=1,6
            workmat(i,j) = 0d0
            do k=1,6
              workmat(i,j) = workmat(i,j) + DjDss(i,k) * L(k,j)
            end do
          end do
        end do
        do i=1,2
          do j=1,6
            DphiDs(i,j) = 0d0
            do k=1,3
              DphiDs(i,j) = DphiDs(i,j) + DphiDj(i,k) * workmat(k,j)
            end do
          end do
        end do
c
        DseDphi = se / (2d0*k_p*phi)
c
        do i=1,6
          DphiDs(1,i) = coef(1)*DphiDs(1,i)+coef(2)*DphiDs(2,i)
          dseds(i) = DseDphi * DphiDs(1,i)
        end do
      end if
c                                            ---- 2nd order differential
      if ( nreq.ge.2 ) then
        if (eqFlag.eq.0) then
          X12 = (largeS(1)-largeS(2))**(2*k_p-2)
          X23 = (largeS(2)-largeS(3))**(2*k_p-2)
          X13 = (largeS(1)-largeS(3))**(2*k_p-2)
c
          DDphiDDls(1,1,1) =  2*k_p*(2*k_p-1)*(X12+X13)
          DDphiDDls(1,2,2) =  2*k_p*(2*k_p-1)*(X12+X23)
          DDphiDDls(1,3,3) =  2*k_p*(2*k_p-1)*(X13+X23)
          DDphiDDls(1,1,2) = -2*k_p*(2*k_p-1)*X12
          DDphiDDls(1,2,3) = -2*k_p*(2*k_p-1)*X23
          DDphiDDls(1,1,3) = -2*k_p*(2*k_p-1)*X13
          do i=1,3
            do j=i+1,3
              DDphiDDls(1,j,i) = DDphiDDls(1,i,j)
            end do
          end do
c
          do i=1,3
            do j=1,3
              DDphiDDls(2,i,j) = 0d0
            end do
            DDphiDDls(2,i,i) = 2*k_p*(2*k_p-1)*largeS(i)**(2*k_p-2)
          end do
c
          do i=1,3
            dum = 1d0 /(3d0*largeS(i)**2 - 2d0*Jinvar(1)*largeS(i)
     &         + Jinvar(2))**3
            DDlsDDj(i,1,1) = dum*largeS(i)**3*
     &         (6d0*largeS(i)**2-6d0*Jinvar(1)*largeS(i)+4d0*Jinvar(2))
            DDlsDDj(i,1,2) = dum*largeS(i)**2*
     &         (-3d0*largeS(i)**2+4d0*Jinvar(1)*largeS(i)-3d0*Jinvar(2))
            DDlsDDj(i,1,3) = 2d0*dum*(-Jinvar(1)*largeS(i)**2
     &         +Jinvar(2)*largeS(i))
            DDlsDDj(i,2,2) = DDlsDDj(i,1,3)
            DDlsDDj(i,2,3) = dum*(3d0*largeS(i)**2-Jinvar(2))
            DDlsDDj(i,3,3) = -dum*(6d0*largeS(i)-2d0*Jinvar(1))
            DDlsDDj(i,2,1) = DDlsDDj(i,1,2)
            DDlsDDj(i,3,1) = DDlsDDj(i,1,3)
            DDlsDDj(i,3,2) = DDlsDDj(i,2,3)
          end do
c
          do i=1,2
            do j=1,3
              do k=j,3
                DDphiDDj(i,j,k) = 0d0
                do m=1,3
                  do n=1,3
                    DDphiDDj(i,j,k) = DDphiDDj(i,j,k)
     &                 + DDphiDDls(i,m,n)*DlsDj(m,j)*DlsDj(n,k)
                  end do
                  DDphiDDj(i,j,k) = DDphiDDj(i,j,k)
     &               + DphiDls(i,m)*DDlsDDj(m,j,k)
                end do
              end do
            end do
            do j=1,3
              do k=j+1,3
                DDphiDDj(i,k,j) = DDphiDDj(i,j,k)
              end do
            end do
          end do
        else
          if (k_p.eq.1) then
            do i=1,3
              do j=1,3
                DDphiDDj(1,i,j) = 0d0
              end do
            end do
            DDphiDDj(1,1,1) = 4d0
          else if (k_p.eq.2) then
            do i=1,3
              do j=1,3
                DDphiDDj(1,i,j) = 0d0
              end do
            end do
            DDphiDDj(1,1,1) = 24d0 * (Jinvar(1)**2 - Jinvar(2))
            DDphiDDj(1,1,2) = -24d0 * Jinvar(1)
            DDphiDDj(1,2,1) = DDphiDDj(1,1,2)
            DDphiDDj(1,2,2) = 36d0
          else
            dum = 2d0 * k_p * (largeS(1)-largeS(3))**(2*(k_p-3))
            DDphiDDj(1,1,1) = (k_p*(2*k_p-1)*(2*k_p+1)*largeS(1)**4/3d0
     &         -4*k_p*(2*k_p-1)*largeS(1)**3*largeS(3)
     &         -4*(2*k_p-1)*(k_p-2)*largeS(1)**2*largeS(3)**2
     &         +8*(k_p-2)*largeS(1)*largeS(3)**3
     &         +2*(2*k_p-1)*largeS(3)**4) * dum
            DDphiDDj(1,1,2) = (-2*k_p*(2*k_p-1)*(k_p-1)*largeS(1)**3/3d0
     &         +(2*k_p-1)*(5*k_p-4)*largeS(1)**2*largeS(3)
     &         +4*(k_p-2)**2*largeS(1)**1*largeS(3)**2
     &         -6*(k_p-1)*largeS(3)**3) * dum
            DDphiDDj(1,1,3) = dum *
     &         ((2*k_p-1)*(2*k_p**2-11*k_p+6)*largeS(1)**2/3d0
     &         -2*(k_p-2)*(2*k_p-3)*largeS(3)*(largeS(1)+largeS(3)))
            DDphiDDj(1,2,2) = dum *
     &         ((k_p-1)*(2*k_p-1)*(2*k_p-3)*largeS(1)**2/3d0
     &         -(12*k_p**2-22*k_p+14)*largeS(1)*largeS(3)
     &         +(10*k_p-11)*largeS(3)**2)
            DDphiDDj(1,2,3) = dum *
     &         ((-4*k_p**3+30*k_p**2-44*k_p+24)*largeS(1)/3d0
     &         +3*(k_p-2)*(2*k_p-3)*largeS(3))
            DDphiDDj(1,3,3) = dum * (4*k_p**3-48*k_p**2+107*k_p-78)/3d0
            do i=1,3
              do j=i+1,3
                DDphiDDj(1,j,i) = DDphiDDj(1,i,j)
              end do
            end do
          end if
c
          do i=1,3
            do j=i,3
              DDphiDDj(2,i,j) = 0d0
            end do
          end do
          do i=0,2*k_p-2
            DDphiDDj(2,1,1) = DDphiDDj(2,1,1) +
     &         k_p*(i+1)*(i+2)*(i+3)*(2*k_p-i-1)/3*
     &         largeS(1)**i * largeS(3)**(2*k_p-2-i)
          end do
          do i=0,2*k_p-3
            DDphiDDj(2,1,2) = DDphiDDj(2,1,2) -
     &         k_p*(i+1)*(i+2)*(i+3)*(2*k_p-i-2)/3*
     &         largeS(1)**i * largeS(3)**(2*k_p-3-i)
          end do
          do i=0,2*k_p-4
            DDphiDDj(2,1,3) = DDphiDDj(2,1,3) +
     &         k_p*(i+1)*(i+2)*(i+3)*(2*k_p-i-3)/3*
     &         largeS(1)**i * largeS(3)**(2*k_p-4-i)
          end do
          DDphiDDj(2,2,2) = DDphiDDj(2,1,3)
          do i=0,2*k_p-5
            DDphiDDj(2,2,3) = DDphiDDj(2,2,3) -
     &         k_p*(i+1)*(i+2)*(i+3)*(2*k_p-i-4)/3*
     &         largeS(1)**i * largeS(3)**(2*k_p-5-i)
          end do
          do i=0,2*k_p-6
            DDphiDDj(2,3,3) = DDphiDDj(2,3,3) +
     &         k_p*(i+1)*(i+2)*(i+3)*(2*k_p-i-5)/3*
     &         largeS(1)**i * largeS(3)**(2*k_p-6-i)
          end do
          do i=1,3
            do j=i+1,3
              DDphiDDj(2,j,i) = DDphiDDj(2,i,j)
            end do
          end do
        end if
c
        do i=1,3
          do j=1,6
            do k=1,6
              DDjDDss(i,j,k) = 0d0
            end do
          end do
        end do
        DDjDDss(2,1,2) = 1d0
        DDjDDss(2,1,3) = 1d0
        DDjDDss(2,2,3) = 1d0
        DDjDDss(2,4,4) = -2d0
        DDjDDss(2,5,5) = -2d0
        DDjDDss(2,6,6) = -2d0
        DDjDDss(3,1,2) = smallS(3)
        DDjDDss(3,1,3) = smallS(2)
        DDjDDss(3,2,3) = smallS(1)
        DDjDDss(3,1,5) = -2d0 * smallS(5)
        DDjDDss(3,2,6) = -2d0 * smallS(6)
        DDjDDss(3,3,4) = -2d0 * smallS(4)
        DDjDDss(3,4,4) = -2d0 * smallS(3)
        DDjDDss(3,5,5) = -2d0 * smallS(1)
        DDjDDss(3,6,6) = -2d0 * smallS(2)
        DDjDDss(3,4,5) = 2d0 * smallS(6)
        DDjDDss(3,4,6) = 2d0 * smallS(5)
        DDjDDss(3,5,6) = 2d0 * smallS(4)
        do i=1,6
          do j=i+1,6
            DDjDDss(2,j,i) = DDjDDss(2,i,j)
            DDjDDss(3,j,i) = DDjDDss(3,i,j)
          end do
        end do
c
        do i=1,2
          do j=1,6
            do k=1,6
              workmat1(i,j,k) = 0d0
              do m=1,3
                do n=1,3
                  workmat1(i,j,k) = workmat1(i,j,k)
     &               + DDphiDDj(i,m,n)*DjDss(m,j)*DjDss(n,k)
                end do
                workmat1(i,j,k) = workmat1(i,j,k)
     &             + DphiDj(i,m)*DDjDDss(m,j,k)
              end do
            end do
          end do
        end do
c
        do i=1,2
          do j=1,6
            do k=1,6
              DDphiDDs(i,j,k) = 0d0
              do m=1,6
                do n=1,6
                  DDphiDDs(i,j,k) = DDphiDDs(i,j,k)
     &               + workmat1(i,m,n)*L(m,j)*L(n,k)
                end do
              end do
            end do
          end do
        end do
c
        do j=1,6
          do k=1,6
            DDphiDDs(1,j,k) = coef(1)*DDphiDDs(1,j,k)
     &         + coef(2)*DDphiDDs(2,j,k)
          end do
        end do
c
        do i=1,6
          do j=i,6
            d2seds2(i,j) = (1-2*k_p)*se/(4*k_p**2*phi**2)
     &         *DphiDs(1,i)*DphiDs(1,j) + se/(2*k_p*phi)*DDphiDDs(1,i,j)
          end do
        end do
        do i=1,6
          do j=i+1,6
            d2seds2(j,i) = d2seds2(i,j)
          end do
        end do
      end if
c
      return
      end
c
c
c
c------------------------------------------------------(karafillisboyce)
c     Find principal stress and invariants.
c     Solving cubic equation by Francois Viete method.
c
      subroutine jancae_KarafillisBoyce_principalStress
     &                                 (stress,invar, pStress)
c-----------------------------------------------------------------------
      implicit none
      real*8, parameter :: PI=3.141592653589793d0, TOL=1e-5
      real*8, intent(inout) :: stress(6)
      real*8, intent(out) :: invar(3), pStress(3)
      real*8 p, q, alpha, c, dum
      invar(1) = stress(1) + stress(2) + stress(3)
      invar(2) = stress(1)*stress(2) + stress(2)*stress(3)
     &   + stress(1)*stress(3) - stress(4)**2 - stress(5)**2
     &   - stress(6)**2
      invar(3) = stress(1)*stress(2)*stress(3)
     &   + 2d0*stress(4)*stress(5)*stress(6) - stress(1)*stress(5)**2
     &   - stress(2)*stress(6)**2 - stress(3)*stress(4)**2
      p = invar(1)**2/9d0 - invar(2)/3d0
      q = invar(1)**3/27d0 + 0.5d0*invar(3) - invar(1)*invar(2)/6d0
      if (p.le.TOL*abs(q)) then
        pStress(1) = (2d0*q)**(1d0/3d0) + invar(1)/3d0
        pStress(2) = pStress(1)
        pStress(3) = pStress(1)
      else
        dum = q/sqrt(p)**3
        if (abs(dum).gt.1d0) then
          if (abs(abs(dum)-1d0).le.TOL) then
            dum = dum/abs(dum)
          else
            call jancae_exit(1000)
          end if
        end if
        alpha = acos(dum)/3d0
        c = 2d0*sqrt(p)
        pStress(1) = c*cos(alpha) + invar(1)/3d0
        pStress(2) = c*cos(alpha+2d0/3d0*PI) + invar(1)/3d0
        pStress(3) = c*cos(alpha+4d0/3d0*PI) + invar(1)/3d0
      end if
      end subroutine jancae_KarafillisBoyce_principalStress
c
c
c
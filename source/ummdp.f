c
c                               OVER THIS LINE IS CODE DEPENDENT
c<<-->><<-->><<-->><<->><<-->><<-->><<-->><<-->><<-->><<-->><<-->><<-->>
c                            UNDER THIS LINE IS CODE INDEPENDENT
c
c     UMMDp: UNIFIED MATERIAL MODEL DRIVER FOR PLASTICITY
c
c
************************************************************************
c
c     PLASTICITY DUMMY
c
      subroutine ummdp_plasticity ( s1,s2,de,p,dp,dpe,de33,x1,x2,mxpbs,
     1                              ddsdde,nnrm,nshr,nttl,nvbs,mjac,
     2                              prop,nprop,propdim )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: mxpbs,nnrm,nshr,nttl,nvbs,nprop,propdim
      real*8 ,intent(in) :: p
      real*8 ,intent(in) :: s1(nttl),de(nttl)
c
      real*8,intent(out) :: dp,de33
      real*8,intent(out) :: s2(nttl),dpe(nttl)
      real*8,intent(out) :: ddsdde(nttl,nttl)
c
      integer,intent(inout) :: mjac
      real*8 ,intent(inout) :: prop(nprop)
      real*8 ,intent(inout) :: x1(mxpbs,nttl),x2(mxpbs,nttl)
c
      integer i,n,ndela,ndyld,ndihd,ndkin,npbs,ndrup,nnn
      character*100 text
c-----------------------------------------------------------------------
c
      if ( prop(1) >= 1000.0d+0 ) then
        mjac = -1
        prop(1) = prop(1) - 1000.0d0
      end if
c
      call ummdp_prop_dim ( prop,nprop,propdim,ndela,ndyld,ndihd,ndkin,
     1                      npbs,ndrup )
c
      n = ndela + ndyld + ndihd + ndkin + ndrup
      if ( n > nprop ) then
        write (6,*) 'nprop error in ummdp_plasticity'
        write (6,*) 'nprop=',nprop
        write (6,*) 'n    =',n
        do i = 1,5
          write (6,*) 'props(',i,')=',prop(i)
        end do
        call ummdp_exit ( 303 )
      end if
      if ( nvbs >= 4 ) then
        write(6,'(//8xA/)') '>> Properties'
        do i = 1,n
          write (6,'(12xA8,I2,A3,E20.12)') '. props(',i,') =',prop(i)
        end do
      end if
c
      nnn = (npbs+1) * nttl
c
      call ummdp_plasticity_core ( s1,s2,de,p,dp,dpe,de33,x1,x2,mxpbs,
     1                             ddsdde,nnrm,nshr,nttl,nvbs,mjac,
     2                             prop,nprop,npbs,ndela,ndyld,ndihd,
     3                             ndkin,ndrup,nnn )
c
      return
      end subroutine ummdp_plasticity
c
c
c
************************************************************************
c
c     PLASTICITY CORE
c
      subroutine ummdp_plasticity_core ( s1,s2,de,p,dp,dpe,de33,x1,x2,
     1                                   mxpbs,ddsdde,nnrm,nshr,nttl,
     2                                   nvbs,mjac,prop,nprop,npbs,
     3                                   ndela,ndyld,ndihd,ndkin,ndrup,
     4                                   nnn )
c
c-----------------------------------------------------------------------
      implicit none
c
      common /ummdp1/ne,ip,lay
      common /ummdpa/n1234
c
      integer,intent(in) :: mxpbs,nnrm,nshr,nttl,nvbs,mjac,nprop,npbs,
     1                      ndela,ndyld,ndihd,ndkin,ndrup,nnn
      real*8,intent(in) :: p
      real*8,intent(in) :: s1(nttl),de(nttl),prop(nprop)
c
      real*8,intent(out) :: de33,dp
      real*8,intent(out) :: s2(nttl),dpe(nttl),ddsdde(nttl,nttl)
     1
      real*8,intent(inout) :: x1(mxpbs,nttl),x2(mxpbs,nttl)
c
      integer ne,ip,lay,n1234,i,j,k,n,m,maxnr,ndiv,maxnest,nout,i1,i2,
     1        j1,j2,k1,k2,nest,newmstg,nite,nstg,mstg,knr,ip1,ip2,nsym
      real*8 tol,xe,se,sy,dsydp,d2sydp2,dpconv,sgapi,sgapb,dsgap,sgap,
     1       pt,g1,g2n,g3nn,top,top0,bot,bot0,ddp,d,a,dd,aa,aaa,sc1,det
      real*8 prela(ndela),pryld(ndyld),prihd(ndihd),prkin(ndkin),
     1       prrup(ndrup),dseds(nttl),stry(nttl),g2(nttl),d33d(nttl),
     2       eta(nttl),xt1(nttl),xt2(nttl),g3n(npbs),s2conv(nttl),
     3       vv(nttl),uv(nttl),v1(nttl),gv(nnn),wv(nnn)
      real*8 delast(nttl,nttl),d2seds2(nttl,nttl),vk(npbs,nttl),
     1       dvkdp(npbs,nttl),g3(npbs,nttl),x2conv(mxpbs,nttl),
     2       em(nttl,nttl),em3(nttl,nttl)
      real*8 dvkds(npbs,nttl,nttl),dvkdxt(npbs,nttl,nttl),am(nnn,nnn),
     1       ami(nnn,nnn),um(nttl,nnn),cm(nttl,nnn),em1(nttl,nnn),
     2       em2(nttl,nnn),bm(nnn,nttl)
      real*8 dvkdx(npbs,npbs,nttl,nttl)
      character*100 text,tmp
c
c-----------------------------------------------------------------------
c     >>> Arguments List
c
c       ne      | index number of element                           (in)
c       ip      | index number of integration point                 (in)
c       lay     | index number of layer (shell and membrane)        (in)
c
c       nnrm    | number of normal components                       (in)
c       nshr    | number of shear components                        (in)
c       nttl    | number of total components                        (in)
c
c             problem type  | nnrm | nshr |  stress components
c            ---------------+----+----+---------------------------
c             plane stress  |  2   |  1   |  sx,sy,   txy
c             thin shell    |  2   |  1   |  sx,sy,   txy
c             plane strain  |  3   |  1   |  sx,sy,sz,txy
c             axi-symmetric |  3   |  1   |  sr,sq,sz,trz
c             thick shell   |  2   |  3   |  sx,sy,   txy,tyz,tzx
c             3D solid      |  3   |  3   |  sx,sy,sz,txy,txz,tyx
c
c       npbs    | number of terms for partial back stresses         (in)
c       mxpbs   | array size of terms for partial back stresses     (in)
c
c       s1      | stress before update                              (in)
c       de      | strain increment                                  (in)
c       p       | equivalent plastic strain                         (in)
c       x1      | partial back stress before update                 (in)
c
c       s2      | stress after update                              (out)
c       dp      | equivalent plastic strain increment              (out)
c       dpe     | plastic strain increment components              (out)
c       de33    | total strain increment in thickness direction    (out)
c       ddsdde  | material Jacobian Dds/Dde                        (out)
c       x2      | partial back stress after update                 (out)
c
c       nvbs    | verbose mode                                      (in)
c
c                  0 | error messages only
c                  1 | summary of MsRM
c                  2 | details of MsRM and summary of NR
c                  3 | details of NR
c                  4 | input/output
c                  5 | all status for debug
c
c                  MsRM | Multistage Return Mapping
c                  NR   | Newton-Raphson
c
c       mjac    | flag for material jacobian                        (in)
c
c                  0 | only stress update
c                  1 | update stress and calculate material jacobian
c                 -1 | use elastic matrix (emergency mode)
c
c       nprop   | dimensions of material properties                 (in)
c       prop    | material properties                               (in)
c
c
c     >>> Local Variables List
c
c       delast  | elastic material Jacobian
c       se      | equivalent stress
c       dseds   | 1st derivative of yield function wrt stress
c       d2seds2 | 2nd derivative of yield function wrt stress
c       stry    | trial stress predicted elastically
c       sy      | flow stress (function of equivalent plastic strain)
c       dsydp   | 1st derivative of flow stress wrt
c                  equivalent plastic strain
c       g1      | residual of stress point to yield surface
c       g2      | residual of direction of plastic strain increment to
c                  normal of yield surface (vector)
c       g2n     | residual of direction of plastic strain increment to
c                  normal of yield surface (norm)
c       eta     | stress for yield function {s}-{xt}
c       xt1     | total back stress before update
c       xt2     | total back stress after update
c       vk      | evolution for kinematic hardening dx=dp*vk
c       dvkdp   | derivative of v wrt equivalent plastic strain
c       dvkds   | derivative of v wrt stress
c       dvkdx   | derivative of v wrt partial back stress
c       dvkdxt  | derivative of v wrt total back stress
c       g3      | residual of evolution for back stress (vector)
c       g3n     | residual of evolution for back stress (norm)
c       sgap    | stress gap to be eliminated in multistage steps
c       tol     | covergence tolerance
c       maxnr   | maximum iterations of Newton-Raphson
c       maxnest | maximum trial times of multistage gap reduction
c       ndiv    | division number of multistage
c
c-----------------------------------------------------------------------
c
      tol = 1.0d-8
      maxnr = 25
      ndiv =  5
      maxnest = 10
c
      nout = 0
      if ( n1234 /= 1234 ) then
        n1234 = 1234
        nout = 1
      end if
c
c                                             ---- print ummdp separator
      if ( (nvbs >= 1) .or. (nout /= 0) ) then
        call ummdp_print_ummdp ( )
      end if
c
c                                          ---- copy material properties
      n = 0
      do i = 1,ndela
        n = n + 1
        prela(i) = prop(n)
      end do
      do i = 1,ndyld
        n = n + 1
        pryld(i) = prop(n)
      end do
      do i = 1,ndihd
        n = n + 1
        prihd(i) = prop(n)
      end do
      do i = 1,ndkin
        n = n + 1
        prkin(i) = prop(n)
      end do
      do i = 1,ndrup
        n = n + 1
        prrup(i) = prop(n)
      enddo
c
      if ( nout /= 0 ) then
        write (6,'(/8xA)') '>> Material Parameters'
        call ummdp_print_elastic   ( prela,ndela )
        call ummdp_print_yield     ( pryld,ndyld )
        call ummdp_print_isotropic ( prihd,ndihd )
        call ummdp_print_kinematic ( prkin,ndkin,npbs )
        call ummdp_print_rupture   ( prrup,ndrup )
      end if
c
c                                                           ---- set [U]
      um = 0.0d0
      i1 = 1
      do i2 = 1,npbs+1
        do j = 1,nttl
          k1 = (i1-1)*nttl + j
          k2 = (i2-1)*nttl + j
          if ( i2 == 1 ) then
            um(k1,k2) = 1.0
          else
            um(k1,k2) = -1.0
          end if
        end do
      end do
c
c                                                     ---- default value
      if ( npbs == 0 ) then
        x2 = 0.0d0
        x1 = 0.0d0
      end if
c
      dp = 0.0d0
      dpe = 0.0d0
      de33 = 0.0d0
      do n = 1,npbs
        do i = 1,nttl
          x2(n,i) = x1(n,i)
        end do
      end do  
      call ummdp_backsum ( npbs,xt1,x1,nttl,mxpbs )
      call ummdp_backsum ( npbs,xt2,x2,nttl,mxpbs )
c
c                                                  ---- print out arrays
      if ( nvbs >= 4 ) then
        write(6,'(//8xA)') '>> Input'
        text = 'Stress'
        call ummdp_utility_print1 ( text,s1,nttl,0 )
        text = 'Strain Increment'
        call ummdp_utility_print1 ( text,de,nttl,0 )
        if ( npbs /= 0 ) then
          text = 'Partial Back Stress'
          call ummdp_backprint ( text,npbs,x1,nttl,mxpbs )
          text = 'Total Back Stress'
          call ummdp_utility_print1 ( text,xt1,nttl,0 )
        end if
      end if
c
c                                                ---- elastic prediction
      if ( nvbs >= 5 ) then
        write(6,'(//8xA)') '>> Elastic Prediction'
      end if
c                                                ---- set elastic matrix
      call ummdp_setdelast ( delast,prela,ndela,nttl,nnrm,nshr,d33d )
c
c                                             ---- copy delast to ddsdde
      do i = 1,nttl
        do j = 1,nttl
          ddsdde(i,j) = delast(i,j)
        end do
      end do
      if ( nvbs >= 5 ) then
        text = 'Elastic Matrix'
        call ummdp_utility_print2 ( text,ddsdde,nttl,nttl,0 )
      end if
c
      call ummdp_utility_mv ( vv,ddsdde,de,nttl,nttl )
      do i = 1,nttl
        s2(i) = s1(i) + vv(i)
      end do
      if ( nvbs >= 5 ) then
        text = 'Elastic Predicted Stress'
        call ummdp_utility_print1 ( text,s2,nttl,0 )
      end if
c                                                       ---- back stress
      do i = 1,nttl
        eta(i) = s2(i) - xt2(i)
      end do
c                                                       ---- check yield
      call ummdp_yield ( se,dseds,d2seds2,0,eta,nttl,nnrm,nshr,pryld,
     1                   ndyld )
      call ummdp_isotropic ( sy,dsydp,d2sydp2,0,p,prihd,ndihd )
c
      if ( nvbs >= 3 ) then
        text = 'Plastic Strain'
        call ummdp_utility_print3 ( text,p,0 )
        text = 'Flow Stress'
        call ummdp_utility_print3 ( text,sy,0 )
        text = 'Equivalent Stress'
        call ummdp_utility_print3 ( text,se,0 )
        if ( npbs /= 0 ) then
          call ummdp_yield  ( xe,dseds,d2seds2,0,xt1,nttl,nnrm,nshr,
     1                        pryld,ndyld )
          text = 'Equivalent Back Stress'
          call ummdp_utility_print3 ( text,xe,0 )
        end if
      end if
      if ( se <= sy ) then
        if ( nvbs >= 3 ) write (6,'(/12xA)') 'JUDGE : ELASTIC'
        if ( (nttl == 3) .or. (nttl == 5) ) then
          do i = 1,nttl
            de33 = de33 + d33d(i)*de(i)
          end do
          if ( nvbs >= 4 ) then
            text = 'Thickness Strain'
            call ummdp_utility_print3 ( text,de33,0 )
          end if
        end if
        goto 500
      else
        if ( nvbs >= 3 ) write (6,'(/12xA)') 'JUDGE : PLASTIC'
      end if
c
c                                                   ---- initialize loop
      do i = 1,nttl
        stry(i) = s2(i)
        s2conv(i) = s2(i)
        do j = 1,npbs
          x2(j,i) = x1(j,i)
          x2conv(j,i) = x1(j,i)
        end do
      end do
      dp = 0.0d0
      dpconv = 0.0d0
      nest = 0
      newmstg = 1
      sgapi = se - sy
      sgapb = sgapi
      nite = 0
      nstg = 0
c
  300 continue
      if ( nest > 0 ) then
        if ( nvbs >= 2 ) then
          write (6,*) '********** Nest of Multistage :',nest
        end if
      end if
      mstg = newmstg
      dsgap = sgapb / float(mstg)
      sgap = sgapb
      dp = dpconv
      do i = 1,nttl
        s2(i) = s2conv(i)
        do n = 1,npbs
          x2(n,i) = x2conv(n,i)
        end do
      end do
      call ummdp_backsum ( npbs,xt2,x2,nttl,mxpbs )
c
c                                          ---- start of multistage loop
      do m = 1,mstg
        nstg = nstg+1
        sgapb = sgap
        sgap = sgapb-dsgap
        if ( m == mstg ) sgap = 0.0
        if ( mstg > 1 ) then
          if ( nvbs >= 2 ) then
            write (6,*) '******** Multistage :',m,'/',mstg
            write (6,*) 'gap% in stress =',sgap/sgapi*100.0d0
          end if
        end if
c
c                                      ---- start of Newton-Raphson loop
        knr = 0
        if ( nvbs >= 3 ) then
          write (6,'(//8xA)') ' >> Newton-Raphson Loop'
        end if
c
  100   continue
        knr = knr + 1
        nite = nite + 1
        if ( nvbs >= 3 ) then
          write(tmp,'(I)') knr
          write (6,'(//12xA,A)') 'Iteration : ',adjustl(tmp)
          text = 'Equivalent Plastic Strain Increment'
          call ummdp_utility_print3 ( text,dp,4 )
        end if
c
        pt = p + dp
c                                 ---- equivalent stress and derivatives
        do i = 1,nttl
          eta(i) = s2(i) - xt2(i)
        end do
        call ummdp_yield ( se,dseds,d2seds2,2,eta,nttl,nnrm,nshr,pryld,
     1                     ndyld )
c
        if ( nvbs >= 5 ) then
          text = 'Updated Stress'
          call ummdp_utility_print1 ( text,s2,nttl,4 )
          if ( npbs /= 0 ) then
            text = 'Updated Total Back Stress'
            call ummdp_utility_print1 ( text,xt2,nttl,4 )
            text = 'eta'
            call ummdp_utility_print1 ( text,eta,nttl,4 )
          end if
          text = '1st Yield Function Derivative'
          call ummdp_utility_print1 ( text,dseds,nttl,4 )
          text = '2nd Yield Function Derivative'
          call ummdp_utility_print2 ( text,d2seds2,nttl,nttl,4 )
        end if
c
c                                       ---- flow stress and derivatives
        call ummdp_isotropic ( sy,dsydp,d2sydp2,1,pt,prihd,ndihd )
c
        if ( nvbs >= 5 ) then
          text = 'Plastic Strain'
          call ummdp_utility_print3 ( text,pt,4 )
          text = 'Flow Stress'
          call ummdp_utility_print3 ( text,sy,4 )
          text = 'Hardening'
          call ummdp_utility_print3 ( text,dsydp,4 )
        end if
c                                                          ---- calc. g1
        g1 = se - sy - sgap
c                                                          ---- calc. g2
        call ummdp_utility_mv ( vv,delast,dseds,nttl,nttl )
        do i = 1,nttl
          g2(i) = s2(i) - stry(i) + dp*vv(i)
        end do
        call ummdp_utility_vvs ( g2n,g2,g2,nttl )
        g2n = sqrt(g2n)
c                                                          ---- calc. g3
        if ( npbs /= 0 ) then
          call ummdp_kinematic ( vk,dvkdp,dvkds,dvkdx,dvkdxt,pt,s2,x2,
     1                           xt2,nttl,nnrm,nshr,mxpbs,npbs,prkin,
     2                           ndkin, pryld,ndyld )
          do n = 1,npbs
            do i = 1,nttl
              g3(n,i) = x2(n,i) - x1(n,i) - dp*vk(n,i)
            end do
          end do
          g3nn = 0.0d0
          do n = 1,npbs
            g3n(n) = 0.0d0
            do i = 1,nttl
              g3n(n) = g3n(n) + g3(n,i)*g3(n,i)
            end do
            g3n(n) = sqrt(g3n(n))
            g3nn = g3nn + g3n(n)*g3n(n)
          end do
          g3nn = sqrt(g3nn)
        else
          g3nn = 0.0d0
        end if
c
        if ( nvbs >= 3 ) then
          write(6,'(//16xA)') '> Residuals'
          text = 'Yield Surface'
          call ummdp_utility_print3 ( text,g1,8 )
          text = 'Normality (Norm)'
          call ummdp_utility_print3 ( text,g2n,8 )
          if ( nvbs >= 5 ) then
            text = 'Normality (Vector)'
            call ummdp_utility_print1 ( text,g2,nttl,8 )
          end if
          if ( npbs /= 0 ) then
            if ( nvbs >= 4 ) then
              do n = 1,npbs
                write(text,'(A,I1)') 'Back Stress Evolution ',n
                call ummdp_utility_print3 ( text,g3n(n),8 )
                if ( nvbs >= 5 ) then
                  do i = 1,nttl
                    uv(i) = g3(n,i)
                  end do
                  write(text,'(A,I1)') 'Back Stress Evolution ',n
                  call ummdp_utility_print1 ( text,uv,nttl,8 )
                end if
              end do
            end if
          end if
        end if
c                  ---- calculate dependencies common for NR and Dds/Dde
c                                                              * set [A]
        call ummdp_utility_setunitm ( am,nnn )
        call ummdp_utility_mm ( em,delast,d2seds2,nttl,nttl,nttl )
        do i1 = 1,npbs+1
          do i2 = 1,npbs+1
            do j1 = 1,nttl
              do j2 = 1,nttl
                k1 = (i1-1)*nttl + j1
                k2 = (i2-1)*nttl + j2
                if ( i1 == 1 ) then
                  if ( i2 == 1 ) then
                    am(k1,k2) = am(k1,k2) + dp*em(j1,j2)
                  else
                    am(k1,k2) = am(k1,k2) - dp*em(j1,j2)
                  end if
                else
                  ip1 = i1 - 1
                  ip2 = i2 - 1
                  if ( i2 == 1 ) then
                    am(k1,k2) = am(k1,k2) - dp*dvkds(ip1,j1,j2)
                  else
                    am(k1,k2) = am(k1,k2) - dp*dvkdx(ip1,ip2,j1,j2)
     1                                    - dp*dvkdxt(ip1,j1,j2)
                  end if
                end if
              end do
            end do
          end do
        end do
c                                                           ---- set {W}
        wv = 0.0d0
        do i1 = 1, npbs+1
          do j1 = 1,nttl
            k1 = (i1-1)*nttl + j1
            if ( i1 == 1 ) then
              do k2 = 1,nttl
                wv(k1) = wv(k1) + delast(j1,k2)*dseds(k2)
              end do
            else
              ip1 = i1-1
              wv(k1) = -vk(ip1,j1) - dp*dvkdp(ip1,j1)
            end if
          end do
        end do
c                                                  ---- calculate [A]^-1
        call ummdp_utility_minv ( ami,am,nnn,det )
c                                                     ---- [C]=[U][A]^-1
        call ummdp_utility_mm ( cm,um,ami,nttl,nnn,nnn )
c
c                                                 ---- check convergence
        if ( (abs(g1  /sy) <= tol) .and.
     1       (abs(g2n /sy) <= tol) .and.
     2       (abs(g3nn/sy) <= tol) ) then
c
          if ( nvbs >= 2 ) then
            write (6,'(//16xA)') 'JUDGE: CONVERGENCE'
          end if
          dpconv = dp
          do i = 1,nttl
            s2conv(i) = s2(i)
            do j = 1,npbs
              x2conv(j,i) = x2(j,i)
            end do
          end do
          goto 200
        else
          if ( nvbs >= 2 ) then
            write (6,'(//16xA)') 'JUDGE: NO CONVERGENCE'
          end if
        end if
c                                                         ---- solve ddp
c                                                           ---- set {G}
        do i = 1,nttl
          gv(i) = g2(i)
        end do
        do n = 1,npbs
          do i = 1,nttl
            gv(n*nttl+i) = g3(n,i)
          end do
        end do
c                              ---- ddp=(g1-{m}^T[C]{G})/(H+{m}^T[C]{W})
        call ummdp_utility_mv ( vv,cm,gv,nttl,nnn )
        call ummdp_utility_vvs ( top0,dseds,vv,nttl )
        top = g1 - top0
        call ummdp_utility_mv ( vv,cm,wv,nttl,nnn )
        call ummdp_utility_vvs ( bot0,dseds,vv,nttl )
        bot = dsydp + bot0
        ddp = top / bot
c
c                        ---- update equivalent plastic strain increment
        dp = dp + ddp
        if ( nvbs >= 3 ) then
          write(6,'(//16xA)') '> Update'
          text = 'Equivalent Plastic Strain Increment'
          call ummdp_utility_print3 ( text,ddp,8 )
          text = 'Updated Equivalent Plastic Strain'
          call ummdp_utility_print3 ( text,dp,8 )
        end if
        if ( dp <= 0.0 ) then
          if ( nvbs >= 3 ) then
            write (6,*) 'negative dp is detected.'
            write (6,*) 'multistage is subdivided.'
          end if
          goto 400
        end if
c
c                                     ---- update stress and back stress
        do i1 = 1,npbs+1
          vv = 0.0d0
          do j1 = 1,nttl
            k1 = (i1-1)*nttl + j1
            do k2 = 1,nnn
              vv(j1) = vv(j1) - ami(k1,k2)*(gv(k2)+ddp*wv(k2))
            end do
          end do
          do j1 = 1,nttl
            if ( i1 == 1 ) then
              s2(j1) = s2(j1) + vv(j1)
            else
              x2(i1-1,j1) = x2(i1-1,j1) + vv(j1)
            end if
          end do
        end do
        call ummdp_backsum ( npbs,xt2,x2,nttl,mxpbs )
c
        if ( knr <= maxnr ) goto 100
c                                        ---- end of Newton-Raphson loop
c
  400   continue
        if ( nvbs >= 2 ) then
          write (6,*) 'Newton Raphson loop is over.',knr
          write (6,*) 'convergence is failed.'
        end if
        if ( nest < maxnest ) then
          nest = nest + 1
          newmstg = (mstg-m+1) * ndiv
          goto 300
        else
          write (6,*) 'Nest of multistage is over.',nest
          text = 'Current stress (input)'
          call ummdp_utility_print1 ( text,s1,nttl,0 )
          text = 'Strain Increment (input)'
          call ummdp_utility_print1 ( text,de,nttl,0 )
          write (6,*) 'Equivalent Plastic Strain (input)'
          write (6,*) p
          write (6,*) 'the proposals to fix this error'
          write (6,*) ' reduce the amount of strain per inc.'
          write (6,*) ' increase maxnest in program',maxnest
          write (6,*) ' increase ndiv    in program',ndiv
          write (6,*) ' increase maxnr   in program',maxnr
          write (6,*) ' increase tol     in program',tol
          call ummdp_exit ( 403 )
        end if
c
  200   continue
c
      end do
c                                            ---- end of multistage loop
c
c                                          ---- plastic strain increment
      do i = 1,nttl
        dpe(i) = dp * dseds(i)
      end do
c
c                                               ---- print out converged
      if ( nvbs >= 4 ) then
        write(6,'(//8xA)') '>> Convergence'
        text = 'Stress'
        call ummdp_utility_print1 ( text,s2,nttl,0 )
        text = 'Plastic Strain Increment'
        call ummdp_utility_print1 ( text,dpe,nttl,0 )
        if ( npbs /= 0 ) then
          text = 'Partial Back Stress'
          call ummdp_backprint ( text,npbs,x2,nttl,mxpbs )
          text = 'Total Back Stress'
          call ummdp_utility_print1 ( text,xt2,nttl,0 )
        end if
      end if
c
c                                        ---- thickness strain increment
      if ( (nttl == 3) .or. (nttl == 5) ) then
        de33 = -dpe(1) - dpe(2)
        do i = 1,nttl
          de33 = de33 + d33d(i)*(de(i) - dpe(i))
        end do
        if ( nvbs >= 4 ) then
          text = 'Thickness Strain'
          call ummdp_utility_print3 ( text,de33,0 )
        end if
      end if
c
c
      if ( nvbs >= 1 ) then
        if ( nest /= 0 ) then
          write (6,*) 'Nest        :',nest
          write (6,*) 'Total Stages        :',nstg
          write (6,*) 'Total NR Ierations  :',nite
          write (6,*) 'Initial Stress Gap         :',sgapi
          write (6,*) 'Increment of Equivalent Plastic Strain :',dp
          write (6,*) 'Updated Equivalent Plastic Strain :',p+dp
        end if
      end if
c
      if ( mjac == 0 ) then
        ddsdde = 0.0d0
        goto 500
      end if
c
      if ( mjac == -1 ) then
        do i = 1,nttl
          do j = 1,nttl
            ddsdde(i,j) = delast(i,j)
          end do
        end do
        goto 500
      end if
c
c                                      ---- consistent material jacobian
      if ( nvbs >= 4 ) then
        write(6,'(//8xA)') '>> Material Jacobian'
      end if
c                                                           ---- set [B]
      bm = 0.0d0
      i1 = 1
      i2 = 1
      do j1 = 1,nttl
        do j2 = 1,nttl
          k1 = (i1-1)*nttl + j1
          k2 = (i2-1)*nttl + j2
          bm(k1,k2) = delast(j1,j2)
        end do
      end do
c                                                       ---- [M1]=[N][C]
      call ummdp_utility_mm ( em1,d2seds2,cm,nttl,nnn,nttl )
c                                               ---- {V1}={m}-dp*[M1]{W}
      call ummdp_utility_mv ( vv ,em1,wv,nttl,nnn )
      do i = 1,nttl
        v1(i) = dseds(i) - dp*vv(i)
      end do
c                                                 ---- [M2]={V1}{m}^T[C]
      em2 = 0.0d0
      do i = 1,nttl
        do j = 1,nnn
          do k = 1,nttl
            em2(i,j) = em2(i,j) + v1(i)*dseds(k)*cm(k,j)
          end do
        end do
      end do
c                                                  ---- S1=H+{m}^T[C]{W}
      sc1 = dsydp
      do i = 1,nttl
        do j = 1,nnn
          sc1 = sc1 + dseds(i)*cm(i,j)*wv(j)
        end do
      end do
c                                     ---- [M3]=[I]-[dp*[M1]-[M2]/S1][B]
      call ummdp_utility_setunitm ( em3,nttl )
      do i = 1,nttl
        do j = 1,nttl
          do k = 1,nnn
            em3(i,j) = em3(i,j) - (dp*em1(i,k)+em2(i,k)/sc1)*bm(k,j)
          end do
        end do
      end do
c                                                     ---- [Dc]=[De][M3]
      call ummdp_utility_mm ( ddsdde,delast,em3,nttl,nttl,nttl )
c
c                                                    ---- check symmetry
      nsym = 1
      d = 0.0d0
      a = 0.0d0
      do i = 1,nttl
        do j = i,nttl
          dd = ddsdde(i,j) - ddsdde(j,i)
          aa = 0.5*(ddsdde(i,j)+ddsdde(i,j))
          d = d + dd*dd
          a = a + aa*aa
        end do
      end do
      a = sqrt(d/a)
      if ( a > 1.0d-8 ) then
        if ( nvbs >= 4 ) then
          write (6,'(/12xA,F)') 'Material Jacobian is not symmetric :',a
          text = 'Material Jacobian | No Symmetry'
          call ummdp_utility_print2 ( text,ddsdde,nttl,nttl,0 )
        end if
c                                                    ---- symmetrization
        if ( nsym == 1 ) then
          do i = 1,nttl
            do j = i+1,nttl
              aaa = 0.5d0 * (ddsdde(i,j)+ddsdde(j,i))
              ddsdde(i,j) = aaa
              ddsdde(j,i) = aaa
            end do
          end do
        end if
      end if
c
      if ( nvbs >= 4 ) then
        text = 'Material Jacobian'
        call ummdp_utility_print2 ( text,ddsdde,nttl,nttl,0 )
      end if
c
  500 if ( (nvbs >= 1) .or. (nout /= 0) ) then
        call ummdp_print_ummdp ( )
      end if
c
      return
      end subroutine ummdp_plasticity_core
c
c
c
************************************************************************
c
c     SET DEBUG AND PRINT VERBOSE MODE
c
      subroutine ummdp_debugmode ( nvbs,nvbs0 )
c
c-----------------------------------------------------------------------
      implicit none
c
      common /ummdp1/ne,ip,lay
c
      integer,intent(in) :: nvbs0
c
      integer,intent(out) :: nvbs
c
      integer ne,ip,lay,nechk,ipchk,laychk,nchk
c-----------------------------------------------------------------------
c
      nechk = 1     ! element number to be checked
      ipchk = 1     ! integration point number to checked
      laychk = 1    ! layer number to be checked
c
      nvbs = 0
      nchk = nechk * ipchk * laychk
      if ( nchk > 0 ) then
        if ( (ne == nechk) .and.
     1       (ip == ipchk) .and.
     2       (lay == laychk) ) then
          nvbs = nvbs0
        end if
      end if
c
      return
      end subroutine ummdp_debugmode
c
c
c
************************************************************************
c
c     SET ELASTIC MATERIAL JACOBIAN MATRIX
c
      subroutine ummdp_setdelast ( delast,prela,ndela,nttl,nnrm,nshr,
     1                             d33d )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: ndela,nttl,nnrm,nshr
      real*8 ,intent(in) :: prela(ndela)
c
      real*8,intent(out) :: delast(nttl,nttl),d33d(nttl)
c
      integer i,j,k,ib,ni,jb,nj,is,i3,js,j3,id,ntela
      real*8 eyoung,epoas,erigid,ek,eg,coef,d33
      real*8 delast3d(6,6)
c-----------------------------------------------------------------------
c
      ntela = nint(prela(1))
      select case ( ntela )
c
      case ( 0 ) ! Young Modulus and Poisson Ratio
        eyoung = prela(2)                           ! Young modulus
        epoas = prela(3)                            ! Poisson ratio
        erigid = eyoung / 2.0d0 / (1.0d0+epoas)     ! Rigidity
c
      case ( 1 ) ! Bulk Modulus and Modulus of Rigidity
        ek = prela(2)                               ! Bulk modulus
        eg = prela(3)                               ! Rigidity
        eyoung = 9.0d0*ek*eg / (3.0d0*ek+eg)        ! Young modulus
        epoas = (eyoung-2.0d0*eg) / 2.0d0 / eg      ! Poisson ratio
        erigid = eg
      case default
        write (6,*) 'error in ummdp_setdelast'
        write (6,*) 'ntela error :',ntela
        call ummdp_exit ( 201 )
      end select
c                                       ---- set 6*6 matrix for 3d solid
      delast3d = 0.0d0
      do i = 1,3
        do j = 1,3
          if ( i == j ) then
            delast3d(i,j) = 1.0d0 - epoas
          else
            delast3d(i,j) = epoas
          end if
        end do
      end do
      do i = 4,6
        delast3d(i,i) = 0.5d0 - epoas
      end do
      coef = erigid / (0.5d0-epoas)
      do i = 1,6
        do j = 1,6
          delast3d(i,j) = coef * delast3d(i,j)
        end do
      end do
c                                       ---- condensation for 2D problem
      do ib = 1,2
        if ( ib == 1 ) then
          ni = nnrm
        else
          ni = nshr
        end if
        do jb = 1,2
          if ( jb == 1 ) then
            nj = nnrm
          else
            nj = nshr
          end if
          do is = 1,ni
            i = (ib-1)*nnrm + is
            i3 = (ib-1)*3 + is
            do js = 1,nj
              j = (jb-1)*nnrm + js
              j3 = (jb-1)*3 + js
              delast(i,j) = delast3d(i3,j3)
            end do
          end do
        end do
      end do
c                                     ---- plane stress or shell element
      if ( nnrm == 2 ) then
        d33 = delast3d(3,3)
        do ib = 1,2
          if ( ib == 1 ) then
            ni = nnrm
          else
            ni = nshr
          end if
          do jb = 1,2
            if ( jb == 1 ) then
              nj = nnrm
            else
              nj = nshr
            end if
            do is = 1,ni
              i = (ib-1)*nnrm + is
              i3 = (ib-1)*3 + is
              do js = 1,nj
                j = (jb-1)*nnrm + js
                j3 = (jb-1)*3 + js
                delast(i,j) = delast(i,j)
     1                        - delast3d(i3,3)*delast3d(3,j3)/d33
              end do
            end do
          end do
        end do
c                 ---- elastic strain coefficient in thickness direction
        do i = 1,nttl
          if ( i <= nnrm ) then
            id = i
          else
            id = i - nnrm + 3
          end if
          d33d(i) = -delast3d(3,id) / d33
        end do
      end if
c
      return
      end subroutine ummdp_setdelast
c
c
c
************************************************************************
c
c     CHECK DIMENSIONS OF INTERNAL STATE VARIABLES
c
      subroutine ummdp_check_nisv ( nisv,nttl,npbs )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: nisv,nttl,npbs
c
      integer isvrsvd,isvsclr,isvtnsr,isvttl
c-----------------------------------------------------------------------
c
      call ummdp_isvprof ( isvrsvd,isvsclr )
c
      if ( npbs == 0 ) then
        isvtnsr = nttl
      else
        isvtnsr = nttl * (1+npbs)
      end if
      isvttl = isvrsvd + isvsclr + isvtnsr
      if ( nisv < isvttl ) then
        write (6,*) 'check number of internal state variables (isv)'
        write (6,*) 'nisv must be larger than',isvttl
        write (6,*) 'nisv=',nisv
        write (6,*) 'isv : required       ',isvttl
        write (6,*) 'isv : system reserved',isvrsvd
        write (6,*) 'isv : for scaler     ',isvsclr
        write (6,*) 'isv : for tensor     ',isvtnsr
        call ummdp_exit ( 302 )
      end if
c
      return
      end subroutine ummdp_check_nisv
c
c
c
************************************************************************
c
c     SET VARIABLES FROM STATE VARIABLES
c
      subroutine ummdp_isv2pex ( isvrsvd,isvsclr,stv,nstv,p,pe,x,nttl,
     1                           mxpbs,npbs )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: nstv,nttl,mxpbs,npbs,isvrsvd,isvsclr
      real*8 ,intent(in) :: stv(nstv)
c
      real*8,intent(out) :: p
      real*8,intent(out) :: pe(nttl)
      real*8,intent(out) :: x(mxpbs,nttl)
c
      integer i,nb,it
c-----------------------------------------------------------------------
c
c                                         ---- equivalent plastic strain
      p = stv(isvrsvd+1)
c
c                                         ---- plastic strain components
      do i = 1,nttl
        pe(i) = stv(isvrsvd + isvsclr + i)
      end do
c
c                                    ---- partial back stress components
      if ( npbs /= 0 ) then
        do nb = 1,npbs
          do i = 1,nttl
            it = isvrsvd + isvsclr + nttl*nb + i
            x(nb,i) = stv(it)
          end do
        end do
      end if
c
      return
      end subroutine ummdp_isv2pex
c
c
c
************************************************************************
c
c     SUM PARTIAL BACK STRESS FOR TOTAL BACK STREESS
c
      subroutine ummdp_backsum ( npbs,xt,x,nttl,mxpbs )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: npbs,nttl,mxpbs
      real*8 ,intent(in) :: x(mxpbs,nttl)
c
      real*8,intent(out) :: xt(nttl)
c
      integer i,j
c-----------------------------------------------------------------------
c
      xt = 0.0d0
      if ( npbs == 0 ) return
c
      do i = 1,nttl
        do j = 1,npbs
          xt(i) = xt(i) + x(j,i)
        end do
      end do
c
      return
      end subroutine ummdp_backsum
c
c
c
************************************************************************
c
c     PRINT BACK STRESS
c
      subroutine ummdp_backprint ( text,npbs,x,nttl,mxpbs )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer      ,intent(in) :: npbs,nttl,mxpbs
      real*8       ,intent(in) :: x(mxpbs,nttl)
      character*100,intent(in) :: text
c
      integer i,j
      real*8 xx(npbs,nttl)
c-----------------------------------------------------------------------
c
      do i = 1,nttl
        do j = 1,npbs
          xx(j,i) = x(j,i)
        end do
      end do
      call ummdp_utility_print2 ( text,xx,npbs,nttl,0 )
c
      return
      end subroutine ummdp_backprint
c
c
c
************************************************************************
c
c     SET DIMENSIONS OF MATERIAL PROPERTIES
c
      subroutine ummdp_prop_dim ( prop,mxprop,propdim,ndela,ndyld,ndihd,
     1                            ndkin,npbs,ndrup )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: mxprop,propdim
      real*8 ,intent(in) :: prop(mxprop)
c
      integer,intent(out) :: ndela,ndyld,ndihd,ndkin,npbs,ndrup
c
      integer n,nd,nela,nyld,nihd,nkin,nrup
      real*8 p
c-----------------------------------------------------------------------
c
      n = 0
      p = prop(n+1)
      if ( p >= 1000.0d0 ) p = p - 1000.d0
c
      nela = nint(p)
      select case ( nela )
        case (0)
          nd = 2
        case (1)
          nd = 2
        case default
          write (6,*) 'Elasticity ID :',nela
          call ummdp_exit ( 201 )
      end select
      ndela = nd + 1
c
      n = ndela
      nyld = nint(prop(n+1))
      select case (nyld)
        case ( 0 ) ! von Mises
          nd = 0
        case ( 1 ) ! Hill 1948
          nd = 6
        case ( 2 ) ! Yld2004-18p
          nd = 19
        case ( 3 ) ! CPB 2005
          nd = 14
        case ( 4 ) ! Karafillis-Boyce
          nd = 8
        case ( 5 ) ! Hu 2005
          nd = 10
        case ( 6 ) ! Yoshida 2011
          nd = 16
c
        case ( -1 ) ! Gotoh
          nd = 9
        case ( -2 ) ! Yld2000-2d
          nd = 9
        case ( -3 ) ! Vegter
          nd = 3 + 4*nint(prop(n+2))
        case ( -4 ) ! BBC 2005
          nd = 9
        case ( -5 ) ! Yld89
          nd = 4
        case ( -6 ) ! BBC 2008
          nd = 2 + 8*nint(prop(n+2))
        case ( -7 ) ! Hill 1990
          nd = 5
        case default
          write (6,*) 'Yield Function ID :',nyld
          call ummdp_exit ( 202 )
      end select
      ndyld = nd + 1
c
      n = ndela + ndyld
      nihd = nint(prop(n+1))
      select case ( nihd )
        case ( 0 ) ! Perfecty Plastic
          nd = 1
        case ( 1 ) ! Linear
          nd = 2
        case ( 2 ) ! Swift
          nd = 3
        case ( 3 ) ! Ludwik
          nd = 3
        case ( 4 ) ! Voce
          nd = 3
        case ( 5 ) ! Voce & Linear
          nd = 4
        case ( 6 ) ! Voce & Swift
          nd = 7
        case ( 7 ) ! p-Model
          nd = 5
        case default
          write (6,*) 'Isotropic Hardening Law ID :',nihd
          call ummdp_exit ( 203 )
      end select
      ndihd = nd + 1
c
      n = ndela + ndyld + ndihd
      nkin = nint(prop(n+1))
      select case ( nkin )
        case ( 0 ) ! None
          nd = 0
          npbs = 0
        case ( 1 ) ! Prager
          nd = 1
          npbs = 1
        case ( 2 ) ! Ziegler
          nd = 1
          npbs = 1
        case ( 3 ) ! Armstrong & Frederick
          nd = 2
          npbs = 1
        case ( 4 ) ! Chaboche I
          npbs = nint(prop(n+2))
          nd = 2*npbs + 1
        case ( 5 ) ! Chaboche II
          npbs = nint(prop(n+2))
          nd = 2*npbs + 1
        case ( 6 ) ! Yoshida-Uemori
          nd = 5
          npbs = 2
        case default
          write (6,*) 'Kinematic Hardening Law ID :',nkin
          call ummdp_exit ( 204 )
c
      end select
      ndkin = nd + 1
c
      n = ndela + ndyld + ndihd + ndkin
      nrup = nint(prop(n+1))
      select case (nrup)
        case ( 0 ) ! None
          nd = 0
        case ( 1 ) ! Equivalent Plastic Strain
          nd = 1 + 1
        case ( 2 ) ! Cockroft & Latham
          nd = 1 + 1
        case ( 3 ) ! Rice & Tracey
          nd = 1 + 1
        case ( 4 ) ! Ayada
          nd = 1 + 1
        case ( 5 ) ! Brozzo
          nd = 1 + 1
        case ( 6 ) ! Forming Limit Diagram
          nd = 1 + 1
        case default
          write (6,*) 'Uncoupled Rupture Criterion ID :',nrup
          call ummdp_exit ( 205 )
      end select
      ndrup = nd + 1
c
      return
      end subroutine ummdp_prop_dim
c
c
c

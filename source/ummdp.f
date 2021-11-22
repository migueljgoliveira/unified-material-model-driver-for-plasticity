c 171230 : start check process by venders' supporters
c 180110 : nela=1 for LS-Dyna like input (bulk & shear modulus)
c
c
c
c
c                               OVER THIS LINE DEPENDS ON CODE
c***********************************************************************
c                            UNDER THIS LINE INDEPENDS ON CODE
c
c     UMMDp  : Unified Material Model Driver for Plasticity
c
c     JANCAE : Japan Association for Nonlinear CAE
c     MMSM   : Material Modeling Sub Meeting
c     MPWG   : Metal Plasticity Working Group
c
c
c
c-----------------------------------------------------------------------
c     this is dummy routine
c
      subroutine jancae_plasticity ( s1,s2,de,
     &                               p,dp,dpe,de33,
     &                               x1,x2,mxpbs,
     &                               ddsdde,
     &                               nnrm,nshr,nttl,
     &                               nvbs,mjac,
     &                               prop,nprop,propdim )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension s1(nttl),s2(nttl),
     &          de(nttl),
     &          dpe(nttl),
     &          x1(mxpbs,nttl),x2(mxpbs,nttl),
     &          ddsdde(nttl,nttl),
     &          prop(nprop)
c
      character text*32
c
      if ( prop(1) >= 1000.0d+0 ) then
        mjac = -1
        prop(1) = prop(1) - 1000.0d+0
      end if
c
      call jancae_prop_dim ( prop,nprop,propdim,
     &                       ndela,ndyld,ndihd,ndkin,
     &                       npbs,ndrup )
c
      n = ndela + ndyld + ndihd + ndkin + ndrup
      if ( n > nprop ) then
        write (6,*) 'nprop error in jancae_plasticity'
        write (6,*) 'nprop=',nprop
        write (6,*) 'n    =',n
        do i = 1,5
          write (6,*) 'prop(',i,')=',prop(i)
        end do
        call jancae_exit ( 9000 )
      end if
      if ( nvbs >= 4 ) then
        do i = 1,n
          write (6,*) 'prop(',i,')=',prop(i)
        end do
      end if
c
      nnn = (npbs+1) * nttl
c
      call jancae_plasticity_core
     &                  ( s1,s2,de,
     &                    p,dp,dpe,de33,
     &                    x1,x2,mxpbs,
     &                    ddsdde,
     &                    nnrm,nshr,nttl,
     &                    nvbs,mjac,
     &                    prop,nprop,
     &                    npbs,ndela,ndyld,ndihd,ndkin,ndrup,nnn )
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     core routine
c
      subroutine jancae_plasticity_core
     &                  ( s1,s2,de,
     &                    p,dp,dpe,de33,
     &                    x1,x2,mxpbs,
     &                    ddsdde,
     &                    nnrm,nshr,nttl,
     &                    nvbs,mjac,
     &                    prop,nprop,
     &                    npbs,ndela,ndyld,ndihd,ndkin,ndrup,nnn )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /jancae1/ne,ip,lay
      common /jancae2/n1234
c
      dimension s1(nttl),s2(nttl),
     &          de(nttl),
     &          dpe(nttl),
     &          x1(mxpbs,nttl),x2(mxpbs,nttl),
     &          ddsdde(nttl,nttl),
     &          prop(nprop)
c
c     arguments list
c
c     ne     : index no. of element
c     ip     : index no. of integration point
c     lay    : index no. of layer (shell and menbrane)
c
c     nnrm   : no. of normal components
c     nshr   : no. of shear  components
c     nttl   : total number of components =nnrm+nshr
c       nnrm and nshr indicate the type of problem
c
c       problem type  |nnrm|nshr | stress comp.
c       --------------+----+-----+-----------------------
c       plane stress  | 2  |  1  | sx,sy,   txy
c       thin shell    | 2  |  1  | sx,sy,   txy
c       plane strain  | 3  |  1  | sx,sy,sz,txy
c       axi-symmetric | 3  |  1  | sr,sq,sz,trz
c       thick shell   | 2  |  3  | sx,sy,   txy,tyz,tzx
c       3D solid      | 3  |  3  | sx,sy,sz,txy,txz,tyx
c
c     npbs   : number terms for partial back stresses
c     mxpbs  : array size of terms for partial back stresses
c
c     s1     : stress before update                   (input)
c     s2     : stress after  update                  (output)
c     de     : strain increment                       (input)
c     p      : equivalent plastic strain              (input)
c              (energetic conjugate to equivalent stress)
c     dp     : equivalent plastic strain inc.        (output)
c     dpe    : plastic strain inc. component         (output)
c     de33   : strain inc. in thickness direction    (output)
c              (for plane stress thin/thick shells)
c     x1     : partial back stress before update      (input)
c     x2     : partial back stress after  update     (output)
c
c     ddsdde : material Jacobian Dds/Dde             (output)
c
c     nvbs   : verbose mode                           (input)
c              0  error message only
c              1  summary of MsRM
c              2  detail of MsRM and summary of NR
c              3  detail of NR
c              4  input/output
c              5  all status for debug
c                 MsRM : Multistage Return Mapping
c                 NR   : Newton-Raphson
c     mjac   : flag for material jacobian             (input)
c              0  only stress update
c              1  update stress and calc. material jacobian
c             -1  use elastic matrix (emergency mode)
c
c     nprop  : dimensions of prop                     (input)
c     prop   : material parameters                    (input)
c
c     these variables are written in Voigt notation
c
c
      dimension prela(ndela),pryld(ndyld),prihd(ndihd),prkin(ndkin),
     &          prrup(ndrup)
c
      dimension delast(nttl,nttl),
     &          dseds(nttl),d2seds2(nttl,nttl),
     &          stry(nttl),g2(nttl),
     &          d33d(nttl)
c
      dimension eta(nttl),xt1(nttl),xt2(nttl),
     &          vk(npbs,nttl),dvkdp(npbs,nttl),dvkds(npbs,nttl,nttl),
     &          dvkdx(npbs,npbs,nttl,nttl),dvkdxt(npbs,nttl,nttl),
     &          g3(npbs,nttl),g3n(npbs)
c
c     local variables list
c
c     delast  : elastic material Jacobian
c     se      : equivalent stress
c     dseds   : dse/ds 1st order differential of eq.stress
c                             with respect to stress components
c     d2seds2 : d2se/ds2 2nd order differential of eq.stress
c                             with respect to stress components
c     stry    : trial stress predicted elastically
c     sy      : flow stress (function of eq.plast.strain)
c     dsydp   : dsy/dp 1st order differential of flow stress
c                             with respect to eq.plast.strain
c     g1      : error of stress point to yield surface
c     g2      : error of direction of plastic strain inc. to
c                       normal of yield surface (error vector)
c     g2n     : norm of g2 vector
c
c     eta     : stress for yield function {s}-{xt}
c     xt1     : total back stress before update
c     xt2     : total back stress after update
c     vk      : eq. of evolution for kinematic hardening dx=dp*vk
c     dvkdp   : dvk/dp differential of v w.r.t eq.plast.strain
c     dvkds   : dvk/ds differential of v w.r.t stress
c     dvkdx   : dvk/dx differential of v w.r.t partial back stress
c     dvkdxt  : dvk/dX differential of v w.r.t total back stress
c     g3      : error of eq. of evolution for back stress
c               (error vector)
c     g3n     : norm of g3 vectors
c
c     sgap    : stress gap to be eliminated in multistage steps
c     tol     : covergence tolerance (default : 1.0d-5 )
c     maxnr   : max iteration of Newton-Raphson
c     maxnest : max trial times of multistage gap reduction
c     ndiv    : division number of multistage
c               max no. of multistage is ndiv^maxnest
c
      dimension s2conv(nttl),x2conv(mxpbs,nttl)
      dimension vv(nttl),uv(nttl),v1(nttl),
     &          gv(nnn),wv(nnn),
     &          em(nttl,nttl),em3(nttl,nttl),
     &          am(nnn,nnn),ami(nnn,nnn),
     &          um(nttl,nnn),cm(nttl,nnn),
     &          em1(nttl,nnn),em2(nttl,nnn),
     &          bm(nnn,nttl)
      character text*32
      logical   debug
c
c
      debug = .true.
      debug = .false.
      tol = 1.0d-5    ! tolerrance of convergence
      maxnr = 25      ! max iterations of Newton-Raphson loop
      ndiv =  5       ! division of multistage loop
      maxnest = 10    ! max re-division of multistage loop
c
      nout = 0
      if ( n1234 /= 1234 ) then
        n1234 = 1234
        nout = 1
      end if
c
      if ( ( nvbs >= 1 ) .or. ( nout /= 0 ) ) then
        write (6,*)
        write (6,*) '**************************************'
        write (6,*) '******* START OF JANCAE/UMMDp ********'
        write (6,*) '**************************************'
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
        write (6,*)
        write (6,*) 'MATERIAL DATA LIST --------------------'
        call jancae_elast_print     ( prela,ndela )
        call jancae_yfunc_print     ( pryld,ndyld )
        call jancae_harden_print    ( prihd,ndihd )
        call jancae_kinematic_print ( prkin,ndkin,npbs )
        call jancae_rupture_print   ( prrup,ndrup )
      end if
c                                                           ---- set [U]
      call jancae_clear2( um,nttl,nnn )
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
c                                                     ---- default value
      if ( npbs == 0 ) then
        do n = 1,mxpbs
          do i = 1,nttl
            x2(n,i) = 0.0
            x1(n,i) = 0.0
          end do
        end do
      end if
c
      de33 = 0.0
      dp = 0.0
      do i = 1,nttl
        dpe(i) = 0.0
      end do
      do n = 1,npbs
        do i = 1,nttl
          x2(n,i) = x1(n,i)
        end do
      end do
      call jancae_backsum ( npbs,xt1,x1,nttl,mxpbs )
      call jancae_backsum ( npbs,xt2,x2,nttl,mxpbs )
c
c                                                  ---- print out arrays
      if ( nvbs >= 4 ) then
        text = 'current stress (input)'
        call jancae_print1 ( text,s1,nttl )
        text = 'strain inc. (input)'
        call jancae_print1 ( text,de,nttl )
        if ( npbs /= 0 ) then
          text = 'part. back stess (input)'
          call jancae_backprint ( text,npbs,x1,nttl,mxpbs )
          text = 'total  back stess (input)'
          call jancae_print1 ( text,xt1,nttl )
        end if
      end if
c                                            ---- set elastic [D] matrix
      call jancae_setdelast ( delast,prela,ndela,
     &                        nttl,nnrm,nshr,d33d )
c                                             ---- copy delast to ddsdde
      do i = 1,nttl
        do j = 1,nttl
          ddsdde(i,j) = delast(i,j)
        end do
      end do
      if ( nvbs >= 5 ) then
        text = 'elastic matrix'
        call jancae_print2 ( text,ddsdde,nttl,nttl )
      end if
c                                                ---- elastic prediction
      call jancae_mv ( vv,ddsdde,de,nttl,nttl )
      do i = 1,nttl
        s2(i) = s1(i) + vv(i)
      end do
      if ( nvbs >= 5 ) then
        text = 'elastic predicted stress'
        call jancae_print1 ( text,s2,nttl )
      end if
c                                                       ---- back stress
      do i = 1,nttl
        eta(i) = s2(i) - xt2(i)
      end do
c                                                       ---- check yield
      call jancae_yfunc  ( se,dseds,d2seds2,0,
     &                     eta,nttl,nnrm,nshr,
     &                     pryld,ndyld )
      call jancae_hardencurve ( sy,dsydp,d2sydp2,
     &                          0,p,prihd,ndihd )
c
      if ( nvbs >= 3 ) then
        write (6,*) 'plastic strain p=',p
        write (6,*) 'flow stress   sy=',sy
        write (6,*) 'equiv.stress  se=',se
        if ( npbs /= 0 ) then
          call jancae_yfunc  ( xe,dseds,d2seds2,0,
     &                         xt1,nttl,nnrm,nshr,
     &                         pryld,ndyld )
          write (6,*) 'equiv.back.s  xe=',xe
        end if
      end if
      if ( se <= sy ) then
        if ( nvbs >= 3 ) write (6,*) 'judge : elastic'
        if ( ( nttl == 3 ) .or. ( nttl == 5 ) ) then
          de33 = 0.0
          do i = 1,nttl
            de33 = de33 + d33d(i)*de(i)
          end do
          if ( nvbs >= 4 ) write (6,*) 'de33=',de33
        end if
        return
      else
        if ( nvbs >= 3 ) write (6,*) 'judge : plastic'
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
      dp = 0.0
      dpconv = 0.0
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
      call jancae_backsum ( npbs,xt2,x2,nttl,mxpbs )
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
        knr = 0
c                                      ---- start of Newton-Raphson loop
        if ( nvbs >= 3 ) then
          write (6,*)
          write (6,*) '**** start of Newton-Raphson loop'
        end if
c
  100   continue
        knr = knr + 1
        nite = nite + 1
        if ( nvbs >= 3 ) then
          write (6,*) '----- NR iteration',knr
          write (6,*) 'inc of p : dp   =',dp
        end if
c
        pt = p + dp
c                                        ---- calc. se and differentials
        do i = 1,nttl
          eta(i) = s2(i) - xt2(i)
        end do
        call jancae_yfunc ( se,dseds,d2seds2,2,
     &                      eta,nttl,nnrm,nshr,
     &                      pryld,ndyld )
c
        if ( nvbs >= 5 ) then
          text = 's2'
          call jancae_print1 ( text,s2,nttl )
          if ( npbs /= 0 ) then
            text = 'xt2'
            call jancae_print1 ( text,xt2,nttl )
            text = 'eta'
            call jancae_print1 ( text,eta,nttl )
          end if
          text = 'dse/ds'
          call jancae_print1 ( text,dseds,nttl )
          text = 'd2se/ds2'
          call jancae_print2 ( text,d2seds2,nttl,nttl )
        end if
c                                        ---- calc. sy and differentials
        call jancae_hardencurve ( sy,dsydp,d2sydp2,
     &                            1,pt,prihd,ndihd )

        if ( nvbs >= 5 ) then
          write (6,*) 'plastic strain p=',pt
          write (6,*) 'flow stress   sy=',sy
          write (6,*) 'hardening dsy/dp=',dsydp
        end if
c                                                          ---- calc. g1
        g1 = se - sy - sgap
c                                                          ---- calc. g2
        call jancae_mv ( vv,delast,dseds,nttl,nttl )
        do i = 1,nttl
          g2(i) = s2(i) - stry(i) + dp*vv(i)
        end do
        call jancae_vvs ( g2n,g2,g2,nttl )
        g2n = sqrt(g2n)
c                                                          ---- calc. g3
        if ( npbs /= 0 ) then
          call jancae_kinematic ( vk,dvkdp,
     &                            dvkds,dvkdx,dvkdxt,
     &                            pt,s2,x2,xt2,
     &                            nttl,nnrm,nshr,
     &                            mxpbs,npbs,
     &                            prkin,ndkin,
     &                            pryld,ndyld )
          do n = 1,npbs
            do i = 1,nttl
              g3(n,i) = x2(n,i) - x1(n,i) - dp*vk(n,i)
            end do
          end do
          g3nn = 0.0
          do n = 1,npbs
            g3n(n) = 0.0
            do i = 1,nttl
              g3n(n) = g3n(n) + g3(n,i)*g3(n,i)
            end do
            g3n(n) = sqrt(g3n(n))
            g3nn = g3nn + g3n(n)*g3n(n)
          end do
          g3nn = sqrt(g3nn)
        else
          g3nn = 0.0
        end if
c
        if ( nvbs >= 3 ) then
          write (6,*) 'g1 (yield surf) =',g1
          write (6,*) 'g2n (normality) =',g2n
          if ( nvbs >= 5 ) then
            text = 'g2 vector'
            call jancae_print1 ( text,g2,nttl )
          end if
          if ( npbs /= 0 ) then
            if ( nvbs >= 4 ) then
              do n = 1,npbs
                write (6,*) 'g3n(',n,')=',g3n(n)
                if ( nvbs >= 5 ) then
                  do i = 1,nttl
                    uv(i) = g3(n,i)
                  end do
                  text = 'g3 vector'
                  call jancae_print1 ( text,uv,nttl )
                end if
              end do
            end if
          end if
        end if
c                      ---- calc. dependencies common for NR and Dds/Dde
c                                                              * set [A]
        call jancae_setunitm ( am,nnn )
        call jancae_mm ( em,delast,d2seds2,nttl,nttl,nttl )
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
     &                                    - dp*dvkdxt(ip1,j1,j2)
                  end if
                end if
              end do
            end do
          end do
        end do
c                                                           ---- set {W}
        call jancae_clear1( wv,nnn )
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
c                                                      ---- calc. [A]^-1
        call jancae_minv ( ami,am,nnn,det )
c                                                     ---- [C]=[U][A]^-1
        call jancae_mm ( cm,um,ami,nttl,nnn,nnn )
c
c
c                                                 ---- check convergence
        if ( ( abs(g1  /sy) <= tol ) .and.
     &       ( abs(g2n /sy) <= tol ) .and.
     &       ( abs(g3nn/sy) <= tol )      ) then
c
          if ( nvbs >= 2 ) then
            write (6,*) '**** Newton-Raphson converged.',knr
          end if
          dpconv = dp
          do i = 1,nttl
            s2conv(i) = s2(i)
            do j = 1,npbs
              x2conv(j,i) = x2(j,i)
            end do
          end do
          goto 200
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
        call jancae_mv ( vv,cm,gv,nttl,nnn )
        call jancae_vvs ( top0,dseds,vv,nttl )
        top = g1 - top0
        call jancae_mv ( vv,cm,wv,nttl,nnn )
        call jancae_vvs ( bot0,dseds,vv,nttl )
        bot = dsydp + bot0
        ddp = top / bot
c                                                         ---- update dp
        dp = dp + ddp
        if ( nvbs >= 3 ) then
          write (6,*) 'modification of dp:ddp=',ddp
          write (6,*) 'updated             dp=',dp
        end if
        if ( dp <= 0.0 ) then
          if ( nvbs >= 3 ) then
            write (6,*) 'negative dp is detected.'
            write (6,*) 'multistage is subdivided.'
          end if
          goto 400
        end if
c                                                  ---- update s2 and x2
        do i1 = 1,npbs+1
          call jancae_clear1( vv,nttl )
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
        call jancae_backsum ( npbs,xt2,x2,nttl,mxpbs )
c
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
          text = 'current stress (input)'
          call jancae_print1 ( text,s1,nttl )
          text = 'strain inc. (input)'
          call jancae_print1 ( text,de,nttl )
          write (6,*) 'eq.plast.strain (input)'
          write (6,*) p
          write (6,*) 'the proposals to fix this error'
          write (6,*) ' reduce the amount of strain per inc.'
          write (6,*) ' increase maxnest in program',maxnest
          write (6,*) ' increase ndiv    in program',ndiv
          write (6,*) ' increase maxnr   in program',maxnr
          write (6,*) ' increase tol     in program',tol
          call jancae_exit ( 9000 )
        end if
c
  200   continue
c
      end do
c                                            ---- end of multistage loop
c
c
c                                                 ---- plast.strain inc.
      do i = 1,nttl
        dpe(i) = dp * dseds(i)
      end do
c                                               ---- print out converged
      if ( nvbs >= 4 ) then
        text = 'updated stress'
        call jancae_print1 ( text,s2,nttl )
        text = 'plastic strain inc'
        call jancae_print1 ( text,dpe,nttl )
        if ( npbs /= 0 ) then
          text = 'updated part. back stess'
          call jancae_backprint ( text,npbs,x2,nttl,mxpbs )
          text = 'updated total back stess'
          call jancae_print1 ( text,xt2,nttl )
        end if
      end if
c                                    ---- calc. strain inc. in thickness
      if ( ( nttl == 3 ) .or. ( nttl == 5 ) ) then
        de33 = -dpe(1) - dpe(2)
        do i = 1,nttl
          de33 = de33 + d33d(i)*(de(i)-dpe(i))
        end do
        if ( nvbs >= 4 ) then
          write (6,*) 'de33=',de33
        end if
      end if
c
      if ( nvbs >= 1 ) then
        if ( nest /= 0 ) then
          write (6,*) 'nest of MsRM               :',nest
          write (6,*) 'total no. of stages        :',nstg
          write (6,*) 'total no. of NR iteration  :',nite
          write (6,*) 'initial stress gap         :',sgapi
          write (6,*) 'inc. of equiv.plast.strain :',dp
          write (6,*) 'equiv.plast.strain updated :',p+dp
          write (6,*) 'location ne,ip,lay         :',ne,ip,lay
        end if
      end if
c
      if ( mjac == 0 ) then
        do i = 1,nttl
          do j = 1,nttl
            ddsdde(i,j) = 0.0
          end do
        end do
        return
      end if
c
      if ( mjac == -1 ) then
        do i = 1,nttl
          do j = 1,nttl
            ddsdde(i,j) = delast(i,j)
          end do
        end do
        return
      end if
c
c                                      ---- consistent material jacobian
c                                                           ---- set [B]
      call jancae_clear2( bm,nnn,nttl )
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
      call jancae_mm ( em1,d2seds2,cm,nttl,nnn,nttl )
c                                               ---- {V1}={m}-dp*[M1]{W}
      call jancae_mv ( vv ,em1,wv,nttl,nnn )
      do i = 1,nttl
        v1(i) = dseds(i) - dp*vv(i)
      end do
c                                                 ---- [M2]={V1}{m}^T[C]
      call jancae_clear2 ( em2,nttl,nnn )
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
      call jancae_setunitm ( em3,nttl )
      do i = 1,nttl
        do j = 1,nttl
          do k = 1,nnn
            em3(i,j) = em3(i,j) - (dp*em1(i,k)+em2(i,k)/sc1)*bm(k,j)
          end do
        end do
      end do
c                                                     ---- [Dc]=[De][M3]
      call jancae_mm ( ddsdde,delast,em3,nttl,nttl,nttl )
c
c                                                    ---- check symmetry
      nsym = 0
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
          write (6,*) 'ddsdde is not symmetric.',a
          text = 'material jacobian (nonsym)'
          call jancae_print2 ( text,ddsdde,nttl,nttl )
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
        text = 'material jacobian (output)'
        call jancae_print2 ( text,ddsdde,nttl,nttl )
      end if
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     set debug and print mode
c
      subroutine jancae_debugmode ( nvbs,nvbs0 )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /jancae1/ne,ip,lay
c                             specify verbose level and point
c      nvbs0 = 0   ! verbose mode
c
c           0  error message only
c           1  summary of MsRM
c           2  detail of MsRM and summary of NR
c           3  detail of NR
c           4  input/output
c           5  all status for debug
c       MsRM : Multistage Return Mapping
c       NR   : Newton-Raphson
c
      nechk = 1     ! element no. to be checked
      ipchk = 1     ! integration point no. to checked
      laychk = 1    ! layer no. to be checked
c
      nvbs = 0
      nchk = nechk * ipchk * laychk
      if ( nchk > 0 ) then
        if ( (ne == nechk ) .and.
     &       (ip == ipchk ) .and.
     &       (lay == laychk ) ) then
          nvbs = nvbs0
        end if
      end if
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     set elastic material jacobian marix
c
      subroutine jancae_setdelast ( delast,prela,ndela,
     &                              nttl,nnrm,nshr,d33d )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension delast(nttl,nttl),prela(ndela),d33d(nttl)
c
      dimension delast3d(6,6)
c
      ntela = nint(prela(1))
      select case ( ntela )
c
      case ( 0:1 )    !  isotropic linear elasticity (Hooke)
c
        if ( ntela == 0 ) then
          eyoung = prela(2)                           ! Young's modulus
          epoas = prela(3)                            ! Poisson's ratio
          erigid = eyoung / 2.0d0 / (1.0d0+epoas)     ! Rigidity
        else
          ek = prela(2)                               ! Bulk modulus
          eg = prela(3)                               ! Rigidity
          eyoung = 9.0d0*ek*eg / (3.0d0*ek+eg)        ! Young's modulus
          epoas = (eyoung-2.0d0*eg) / 2.0d0 / eg      ! Poisson's ratio
          erigid = eg
        end if
c                                       ---- set 6*6 matrix for 3d solid
        call jancae_clear2( delast3d,6,6 )
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
c
      case default                                              ! Error
        write (6,*) 'elasticity code error in jancae_setelast'
        write (6,*) 'ntela=',ntela
        call jancae_exit ( 9000 )
c
      end select
c
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
              j3 = (jb-1)*3    + js
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
     &                      - delast3d(i3,3)*delast3d(3,j3)/d33
              end do
            end do
          end do
        end do
c                         ---- elastic strain in thickness direction e_t
c                                             ---- e_t=SUM(d33d(i)*e(i))
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
      end
c
c
c
c-----------------------------------------------------------------------
c     check dimensions of internal state variables
c
      subroutine jancae_check_nisv ( nisv,nttl,npbs )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      call jancae_isvprof ( isvrsvd,isvsclr )
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
        call jancae_exit ( 9000 )
      end if
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     set variables from state variables
c
      subroutine jancae_isv2pex ( isvrsvd,isvsclr,
     &                            stv,nstv,
     &                            p,pe,x,nttl,mxpbs,npbs )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension stv(nstv),pe(nttl),x(mxpbs,nttl)
c
c                                                 ---- eq.plastic strain
      p = stv(isvrsvd+1)
c                                              ---- plastic strain comp.
      do i = 1,nttl
        pe(i) = stv(isvrsvd+isvsclr+i)
      end do
c                                         ---- partial back stress comp.
      if ( npbs /= 0 ) then
        do nb = 1,npbs
          do i = 1,nttl
            it = isvrsvd + isvsclr + nttl*nb+i
            x(nb,i) = stv(it)
          end do
        end do
      end if
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     sum partial back stress for total back stress
c
      subroutine jancae_backsum ( npbs,xt,x,nttl,mxpbs )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension xt(nttl),x(mxpbs,nttl)
c
      do i = 1,nttl
        xt(i) = 0.0
      end do
      if ( npbs == 0 ) return
c
      do i = 1,nttl
        do j = 1,npbs
          xt(i) = xt(i) + x(j,i)
        end do
      end do
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     print back stress
c
      subroutine jancae_backprint ( text,npbs,x,nttl,mxpbs )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension x(mxpbs,nttl)
      character text*32
      dimension xx(npbs,nttl)
c
      if ( npbs == 0 ) return
c
      do i = 1,nttl
        do j = 1,npbs
          xx(j,i) = x(j,i)
        end do
      end do
      call jancae_print2 ( text,xx,npbs,nttl )
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     set dimensions of material properties
c
      subroutine jancae_prop_dim ( prop,mxprop,propdim,
     &                             ndela,ndyld,ndihd,ndkin,
     &                             npbs,ndrup )
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      dimension prop(mxprop)
c
      n = 0
      p = prop(n+1)
      if ( p >= 1000.0d0 ) p = p - 1000.d0
c
      nela = nint(p)
      select case ( nela )
c
      case (0) ; nd = 2
        case (1) ; nd = 2
c
        case default
          write (6,*) 'error elastic property id :',nela
          call jancae_exit ( 9000 )
c
      end select
      ndela = nd + 1
c
      n = ndela
      nyld = nint(prop(n+1))
      select case (nyld)
c
        case ( 0 ) ; nd = 0                   ! von Mises
        case ( 1 ) ; nd = 6                   ! Hill 1948
        case ( 2 ) ; nd = 19                  ! Yld2004-18p
        case ( 3 ) ; nd = 14                  ! CPB 2005
        case ( 4 ) ; nd = 8                   ! Karafillis-Boyce
        case ( 5 ) ; nd = 10                  ! Hu 2005
        case ( 6 ) ; nd = 16                  ! Yoshida 2011
c
        case ( -1 ) ; nd = 9                  ! Gotoh biquadratic (1978)
        case ( -2 ) ; nd = 9                  ! Yld2000-2d
        case ( -3 )                           ! Vegter
          nd = 3 + 4*nint(prop(n+2))
        case ( -4 ) ; nd = 9                  ! BBC 2005
        case ( -5 ) ; nd = 4                  ! Yld89
        case ( -6 )                           ! BBC 2008
          nd = 2 + 8*nint(prop(n+2))
        case ( -7 ) ; nd = 0.5d0              ! Hill 1990
c
        case default
          write (6,*) 'error yield function id :',nyld
          call jancae_exit ( 9000 )
c
      end select
      ndyld = nd + 1
c
      n = ndela + ndyld
      nihd = nint(prop(n+1))
      select case ( nihd )
c
        case ( 0 ) ; nd = 1                           ! Perfecty Plastic
        case ( 1 ) ; nd = 2                           ! Linear
        case ( 2 ) ; nd = 3                           ! Swift
        case ( 3 ) ; nd = 3                           ! Ludwick
        case ( 4 ) ; nd = 3                           ! Voce
        case ( 5 ) ; nd = 4                           ! Voce + Linear
        case ( 6 ) ; nd = 7                           ! Voce + Swift
c
        case default
          write (6,*) 'error work hardening curve id :',nihd
          call jancae_exit ( 9000 )
c
      end select
      ndihd = nd + 1
c
      n = ndela + ndyld + ndihd
      nkin = nint(prop(n+1))
      select case ( nkin )
c
        case ( 0 ) ; nd = 0 ; npbs = 0          ! No Kinematic Hardening
        case ( 1 ) ; nd = 1 ; npbs = 1          ! Prager
        case ( 2 ) ; nd = 1 ; npbs = 1          ! Ziegler
        case ( 3 ) ; nd = 2 ; npbs = 1          ! Armstrong & Frederick
        case ( 4 )                              ! Chaboche
          npbs = nint(prop(n+2))
          nd = 2*npbs + 1
        case ( 5 )                              ! Chaboche - Ziegler
          npbs = nint(prop(n+2))
          nd = 2*npbs + 1
        case ( 6 ) ; nd = 5 ; npbs = 2          ! Yoshida-Uemori
c
        case default
          write (6,*) 'error kinematic hardening id :',nkin
          call jancae_exit ( 9000 )
c
      end select
      ndkin = nd + 1
c
      n = ndela + ndyld + ndihd + ndkin
      nrup = nint(prop(n+1))
      select case (nrup)
c
        case ( 0 ) ; nd = 0             ! No Uncoupled Rupture Criterion
        case ( 1 ) ; nd = 2             ! Equivalent Plastic Strain
        case ( 2 ) ; nd = 2             ! Cockroft and Latham
        case ( 3 ) ; nd = 2             ! Rice and Tracey
        case ( 4 ) ; nd = 2             ! Ayada
        case ( 5 ) ; nd = 2             ! Brozzo
        case ( 6 ) ; nd = 2             ! Forming Limit Diagram
c
        case default
          write (6,*) 'error uncoupled rupture criterion id :',nrup
          call jancae_exit ( 9000 )
c
      end select
      ndrup = nd + 1
c
      return
      end
c
c
c

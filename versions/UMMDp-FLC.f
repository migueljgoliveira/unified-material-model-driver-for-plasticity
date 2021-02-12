c
c     JANCAE.UMMDp/Simulia.Abaqus
c
c     150908 H.Takizawa mxpbs,prop(nprop)
c     171229 H.Takizawa material data from PROPS(NPROPS)
c            jancae_prop_dim
c            common /jancae3/prop (from UMAT to UVARM)
c
c
c
c
c
c-----------------------------------------------------------------------
c
c
c
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1    RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2    TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3    NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4    DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
c
      INCLUDE 'ABA_PARAM.INC'
c
      CHARACTER*80 MATERL
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),
     4 DFGRD0(3,3),DFGRD1(3,3)
c
c-----------------------------------------------------------------------
      common /jancae1/ne,ip,lay
      common /jancae3/prop
      common /jancaea/nsdv
      common /jancaeb/propdim
      parameter (mxpbs=10)
c
      dimension s2(ntens),dpe(ntens),x1(mxpbs,ntens),x2(mxpbs,ntens),
     &          pe(ntens)
c
      dimension ustatev(6)
c
      parameter (mxprop=100)
      dimension prop(mxprop)
c
c                        ne  : element no.
c                        ip  : integration point no.
c                        lay : layer no. of shell
      ne=noel
      ip=npt
      lay=kspt
      if ( lay.eq.0 ) lay=1
      nsdv=nstatv
      nprop=mxprop
      propdim=nprops-1
c
c                                        ---- set debug and verbose mode
      nvbs0 = props(1)
      call jancae_debugmode ( nvbs,nvbs0 )
c                                       ---- output detailed information
      if ( nvbs.ge.4 ) then
        call jancae_printinfo  ( kinc,ndi,nshr )
        call jancae_printinout ( 0,stress,dstran,ddsdde,ntens,
     &                           statev,nstatv )
      endif
c
c                                           ---- set material properties
      do i=2,nprops
        prop(i-1)=props(i)
      enddo
c
      call jancae_prop_dim ( prop,nprop,propdim,
     &                       ndela,ndyld,ndihd,ndkin,
     &                       npbs )
      if ( npbs.gt.mxpbs ) then
        write (6,*) 'npbs > mxpbs error in umat'
        write (6,*) 'npbs =',npbs
        write (6,*) 'mxpbs=',mxpbs
        call jancae_exit ( 9000 )
      endif
c                                                      ---- check nstatv
      call jancae_check_nisv ( nstatv,ntens,npbs )
c                             ---- copy current internal state variables
      call jancae_isvprof ( isvrsvd,isvsclr )
      call jancae_isv2pex ( isvrsvd,isvsclr,
     &                      statev,nstatv,
     &                      p,pe,x1,ntens,mxpbs,npbs )
c
c                             ---- update stress and set tangent modulus
      mjac=1
      call jancae_plasticity  ( stress,s2,dstran,
     &                          p,dp,dpe,de33,
     &                          x1,x2,mxpbs,
     &                          ddsdde,
     &                          ndi,nshr,ntens,
     &                          nvbs,mjac,
     &                          prop,nprop,propdim )
c                                                     ---- update stress
      do i=1,ntens
        stress(i)=s2(i)
      enddo
c                                            ---- update eq.plast,strain
      statev(isvrsvd+1)=p+dp
c                                         ---- update plast.strain comp.
      call rotsig(statev(isvrsvd+2), drot, ustatev, 2, ndi, nshr)
c
      do i=1,ntens
        is=isvrsvd+isvsclr+i
        statev(is)=ustatev(i)+dpe(i)
      enddo
c                                       ---- update of back stress comp.
      if ( npbs.ne.0 ) then
        do n=1,npbs
          do i=1,ntens
            is=isvrsvd+isvsclr+ntens*n+i
            statev(is)=x2(n,i)
          enddo
        enddo
      endif
c                           ----  if debug mode, output return arguments
      if ( nvbs.ge.4 ) then
        call jancae_printinout ( 1,stress,dstran,ddsdde,ntens,
     &                           statev,nstatv )
      endif
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,
     1                  LAYER,KSPT)
c
      INCLUDE 'ABA_PARAM.INC'
c
      DIMENSION STATEV(NSTATV),COORDS(NCRDS)
c-----------------------------------------------------------------------
c
      ne =noel
      ip =npt
      lay=kspt
      if ( lay.eq.0 ) lay=1
c
      if ( ne*ip*lay.eq.1 ) then
        write (6,*) 'SDVINI is called. '
      endif
c
      do n=1,nstatv
        statev(n)=0.0
      enddo
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c
      SUBROUTINE UVARM(UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1             NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2             JMAC,JMATYP,MATLAYO,LACCFLA)
c
      INCLUDE 'ABA_PARAM.INC'
c
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
c-----------------------------------------------------------------------
C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.
c
      parameter   (maxsdv=50)
      parameter   (mxpbs=10)
      parameter   (mxprop=100)
      common /jancae3/prop
      common /jancaea/nsdv
      common /jancaeb/propdim
c
      dimension s(ndi+nshr),xsum(ndi+nshr),x(mxpbs,ndi+nshr),
     &          pe(ndi+nshr),eta(ndi+nshr),
     &          dseds(ndi+nshr),d2seds2(ndi+nshr,ndi+nshr)
c
      dimension   ARRAY2(maxsdv),JARRAY2(maxsdv)
      character*3 FLGRAY2(maxsdv)
      dimension   sdv(maxsdv)
c
      dimension prop(mxprop)
      real*8,allocatable,dimension(:) :: prela,pryld,prihd,prkin
c
c--------------------------------------------------------- flc variables
      dimension   ARRAY3(3),JARRAY3(3)
      character*3 FLGRAY3(3)
c
c-----------------------------------------------------------------------
      nprop=mxprop
c
c     variables list :
c        uvar codes and state variables arrays
c        nt:ntens
c
c     statev(1                     ) : equivalent plastic strain
c     statev(2        ~  1+nt      ) : plastic strain comp.
c     statev(1+nt*i+1 ~  1+nt*(i+1)) : partial back stress component xi
c
c     uvar(1     ) : equivalent stress
c     uvar(2     ) : flow stress
c     uvar(3~2+nt) : back stress
c     uvar(3~2+nt) : rupture criterion
c
      ne =noel
      ip =npt
      lay=kspt
      if ( lay.eq.0 ) lay=1
      ntens=ndi+nshr
c                                                        ---- get stress
      CALL GETVRM('S',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1            MATLAYO,LACCFLA)
      if ( JRCD.ne.0 ) then
        write (6,*) 'request error in uvarm for s'
        write (6,*) 'stop in uvrm.'
        call jancae_exit ( 9000 )
      endif
c
      do i=1,ndi
        s(i)=array(i)
      enddo
      do i=1,nshr
        i1=ndi+i
        i2=3+i
        s(i1)=array(i2)
      enddo
c                                               ---- get state variables
      if ( nsdv.gt.maxsdv ) then
        write (6,*) 'increase dimension of ARRAY2 and JARRAY2'
        write (6,*) 'stop in uvrm.'
        call jancae_exit ( 9000 )
      endif
      CALL GETVRM('SDV',ARRAY2,JARRAY2,FLGRAY2,JRCD,JMAC,JMATYP,
     1            MATLAYO,LACCFLA)
      if ( JRCD.ne.0 ) then
        write (6,*) 'request error in uvarm for sdv'
        write (6,*) 'stop in uvrm.'
        call jancae_exit ( 9000 )
      endif
      do i=1,nsdv
        sdv(i)=array2(i)
      enddo
c                                           ---- set material properties
      call jancae_prop_dim ( prop,nprop,propdim,
     &                       ndela,ndyld,ndihd,ndkin,
     &                       npbs )
      allocate( prela(ndela) )
      allocate( pryld(ndyld) )
      allocate( prihd(ndihd) )
      allocate( prkin(ndkin) )
      k=0
      do i=1,ndela
        k=k+1
        prela(i)=prop(k)
      enddo
      do i=1,ndyld
        k=k+1
        pryld(i)=prop(k)
      enddo
      do i=1,ndihd
        k=k+1
        prihd(i)=prop(k)
      enddo
      do i=1,ndkin
        k=k+1
        prkin(i)=prop(k)
      enddo
c
c                                                  ---- calc back stress
      call jancae_isvprof ( isvrsvd,isvsclr )
      call jancae_isv2pex ( isvrsvd,isvsclr,
     &                      sdv,maxsdv,
     &                      p,pe,x,ntens,mxpbs,npbs )
      do i=1,ntens
         xsum(i)=0.0
      enddo
      if ( npbs.ne.0 ) then
        do i=1,ntens
          do nb=1,npbs
            xsum(i)=xsum(i)+x(nb,i)
          enddo
        enddo
      endif
c                                                 ---- equivalent stress
      if ( nuvarm.ge.1 ) then
        do i=1,ntens
          eta(i)=s(i)-xsum(i)
        enddo
        call jancae_yfunc ( se,dseds,d2seds2,0,
     &                      eta,ntens,ndi,nshr,
     &                      pryld,ndyld )
        uvar(1)=se
      endif
c                                                       ---- flow stress
      if ( nuvarm.ge.2 ) then
        call jancae_hardencurve ( sy,dsydp,d2sydp2,
     &                            0,p,prihd,ndihd )
        uvar(2)=sy
      endif
c                                                       ---- back stress
      if ( npbs.ne.0 ) then
        if ( nuvarm.ge.3 ) then
          do i=1,ntens
            uvar(2+i)=xsum(i)
          enddo
        endif
      endif
c
c------------------------------------------------- forming limit diagram
c
      CALL GETVRM('LEP',ARRAY3,JARRAY3,FLGRAY3,JRCD,JMAC,JMATYP,
     1            MATLAYO,LACCFLA)
      if ( JRCD.ne.0 ) then
        write (6,*) 'request error in uvarm for ep'
        write (6,*) 'stop in uvrm.'
        call jancae_exit ( 9000 )
      endif
c
      E2 = ARRAY3(1)
      E1 = ARRAY3(2)
c
      D1 = -1.0d0
      D2 =  0.5d0
c
      if ( prela(2).eq.104.0d3 ) then
        mat = 1.0d0
      elseif ( prela(2).eq.210.0d3 ) then
        mat = 2.0d0
      elseif ( prela(2).eq.69.0d3 ) then
        mat = 3.0d0
      endif
c
      select case ( int(mat) )
      case ( 1 )                      ! Cu
        CU_B = 0.129d0
        if ( E2.le.0.0d0 ) then
          E1_FLC = D1*E2 + CU_B
        else
          E1_FLC = D2*E2 + CU_B
        endif
      case ( 2 )                      ! DP600
        DP_B = 0.194d0
        if ( E2.le.0.0d0 ) then
          E1_FLC = D1*E2 + DP_B
        else
          E1_FLC = D2*E2 + DP_B
        endif
      case ( 3 )                      ! AA2090-T3
        AA_B = 0.227d0
        if ( E2.le.0.0d0 ) then
          E1_FLC = D1*E2 + AA_B
        else
          E1_FLC = D2*E2 + AA_B
        endif
      end select

      if (E1.eq.0.d0) THEN
        W_FLD = 0.d0
      else
        W_FLD = E1/E1_FLC
      endif
      uvar(3) = W_FLD
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     set internal state variables profile
c
      subroutine jancae_isvprof ( isvrsvd,isvsclr )
c-----------------------------------------------------------------------
      INCLUDE 'ABA_PARAM.INC'
c
c               no reserved variables
      isvrsvd=0
c               statev(1) is for eq.plast.strain
      isvsclr=1
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     exit program by error
c
      subroutine jancae_exit (nexit)
c-----------------------------------------------------------------------
      INCLUDE 'ABA_PARAM.INC'
      common /jancae1/ne,ip,lay
c                                  nexit : exit code
      write (6,*) 'error code :',nexit
      write (6,*) 'element no.           :',ne
      write (6,*) 'integration point no. :',ip
      write (6,*) 'layer no.             :',lay
c
      call xit
c
      return
      end
c
c
c
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
      if ( prop(1).ge.1000.0d+0 ) then
        mjac=-1
        prop(1)=prop(1)-1000.0d+0
      endif
c
      call jancae_prop_dim ( prop,nprop,propdim,
     &                       ndela,ndyld,ndihd,ndkin,
     &                       npbs )
c
      n=ndela+ndyld+ndihd+ndkin
      if ( n.gt.nprop ) then
        write (6,*) 'nprop error in jancae_plasticity'
        write (6,*) 'nprop=',nprop
        write (6,*) 'n    =',n
        do i=1,5
          write (6,*) 'prop(',i,')=',prop(i)
        enddo
        call jancae_exit ( 9000 )
      endif
      if ( nvbs.ge.4 ) then
        do i=1,n
          write (6,*) 'prop(',i,')=',prop(i)
        enddo
      endif
c
      nnn=(npbs+1)*nttl
c
      call jancae_plasticity_core
     &                  ( s1,s2,de,
     &                    p,dp,dpe,de33,
     &                    x1,x2,mxpbs,
     &                    ddsdde,
     &                    nnrm,nshr,nttl,
     &                    nvbs,mjac,
     &                    prop,nprop,
     &                    npbs,ndela,ndyld,ndihd,ndkin,nnn )
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
     &                    npbs,ndela,ndyld,ndihd,ndkin,nnn )
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
      dimension prela(ndela),pryld(ndyld),prihd(ndihd),prkin(ndkin)
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
      debug  = .true.
      debug  = .false.
      tol    = 1.0d-5  ! tolerrance of convergence
      maxnr  = 25      ! max iterations of Newton-Raphson loop
      ndiv   =  5      ! division of multistage loop
      maxnest= 10      ! max re-division of multistage loop
c
      nout=0
      if ( n1234.ne.1234 ) then
        n1234=1234
        nout=1
      endif
c
      if ( (nvbs.ge.1).or.(nout.ne.0) ) then
        write (6,*)
        write (6,*) '**************************************'
        write (6,*) '******* START OF JANCAE/UMMDp ********'
        write (6,*) '**************************************'
      endif
c
c                                          ---- copy material properties
      n=0
      do i=1,ndela
        n=n+1
        prela(i)=prop(n)
      enddo
      do i=1,ndyld
        n=n+1
        pryld(i)=prop(n)
      enddo
      do i=1,ndihd
        n=n+1
        prihd(i)=prop(n)
      enddo
      do i=1,ndkin
        n=n+1
        prkin(i)=prop(n)
      enddo
c
      if ( nout.ne.0 ) then
        write (6,*)
        write (6,*) 'MATERIAL DATA LIST --------------------'
        call jancae_elast_print     ( prela,ndela )
        call jancae_yfunc_print     ( pryld,ndyld )
        call jancae_harden_print    ( prihd,ndihd )
        call jancae_kinematic_print ( prkin,ndkin,npbs )
      endif
c                                                           ---- set [U]
      call jancae_clear2( um,nttl,nnn )
      i1=1
      do i2=1,npbs+1
        do j=1,nttl
          k1=(i1-1)*nttl+j
          k2=(i2-1)*nttl+j
          if ( i2.eq.1 ) then
            um(k1,k2)= 1.0
          else
            um(k1,k2)=-1.0
          endif
        enddo
      enddo
c                                                     ---- default value
      if ( npbs.eq.0 ) then
        do n=1,mxpbs
          do i=1,nttl
            x2(n,i)=0.0
            x1(n,i)=0.0
          enddo
        enddo
      endif
c
      de33=0.0
      dp=0.0
      do i=1,nttl
        dpe(i)=0.0
      enddo
      do n=1,npbs
        do i=1,nttl
          x2(n,i)=x1(n,i)
        enddo
      enddo
      call jancae_backsum ( npbs,xt1,x1,nttl,mxpbs )
      call jancae_backsum ( npbs,xt2,x2,nttl,mxpbs )
c
c                                                  ---- print out arrays
      if ( nvbs.ge.4 ) then
        text='current stress (input)'
        call jancae_print1 ( text,s1,nttl )
        text='strain inc. (input)'
        call jancae_print1 ( text,de,nttl )
        if ( npbs.ne.0 ) then
          text='part. back stess (input)'
          call jancae_backprint ( text,npbs,x1,nttl,mxpbs )
          text='total  back stess (input)'
          call jancae_print1 ( text,xt1,nttl )
        endif
      endif
c                                            ---- set elastic [D] matrix
      call jancae_setdelast ( delast,prela,ndela,
     &                        nttl,nnrm,nshr,d33d )
c                                             ---- copy delast to ddsdde
      do i=1,nttl
        do j=1,nttl
          ddsdde(i,j)=delast(i,j)
        enddo
      enddo
      if ( nvbs.ge.5 ) then
        text='elastic matrix'
        call jancae_print2 ( text,ddsdde,nttl,nttl )
      endif
c                                                ---- elastic prediction
      call jancae_mv ( vv,ddsdde,de,nttl,nttl )
      do i=1,nttl
        s2(i)=s1(i)+vv(i)
      enddo
      if ( nvbs.ge.5 ) then
        text='elastic predicted stress'
        call jancae_print1 ( text,s2,nttl )
      endif
c                                                       ---- back stress
      do i=1,nttl
        eta(i)=s2(i)-xt2(i)
      enddo
c                                                       ---- check yield
      call jancae_yfunc  ( se,dseds,d2seds2,0,
     &                     eta,nttl,nnrm,nshr,
     &                     pryld,ndyld )
      call jancae_hardencurve ( sy,dsydp,d2sydp2,
     &                          0,p,prihd,ndihd )
c
      if ( nvbs.ge.3 ) then
        write (6,*) 'plastic strain p=',p
        write (6,*) 'flow stress   sy=',sy
        write (6,*) 'equiv.stress  se=',se
        if ( npbs.ne.0 ) then
          call jancae_yfunc  ( xe,dseds,d2seds2,0,
     &                         xt1,nttl,nnrm,nshr,
     &                         pryld,ndyld )
          write (6,*) 'equiv.back.s  xe=',xe
        endif
      endif
      if ( se.le.sy ) then
        if ( nvbs.ge.3 ) write (6,*) 'judge : elastic'
        if ( (nttl.eq.3).or.(nttl.eq.5) ) then
          de33=0.0
          do i=1,nttl
            de33=de33+d33d(i)*de(i)
          enddo
          if ( nvbs.ge.4 ) write (6,*) 'de33=',de33
        endif
        return
      else
        if ( nvbs.ge.3 ) write (6,*) 'judge : plastic'
      endif
c
c                                                   ---- initialize loop
      do i=1,nttl
        stry(i)  =s2(i)
        s2conv(i)=s2(i)
        do j=1,npbs
          x2(j,i)    =x1(j,i)
          x2conv(j,i)=x1(j,i)
        enddo
      enddo
      dp=0.0
      dpconv=0.0
      nest=0
      newmstg=1
      sgapi=se-sy
      sgapb=sgapi
      nite=0
      nstg=0
c
  300 continue
      if ( nest.gt.0 ) then
        if ( nvbs.ge.2 ) then
          write (6,*) '********** Nest of Multistage :',nest
        endif
      endif
      mstg=newmstg
      dsgap=sgapb/float(mstg)
      sgap =sgapb
      dp=dpconv
      do i=1,nttl
        s2(i)=s2conv(i)
        do n=1,npbs
          x2(n,i)=x2conv(n,i)
        enddo
      enddo
      call jancae_backsum ( npbs,xt2,x2,nttl,mxpbs )
c
c                                          ---- start of multistage loop
      do m=1,mstg
        nstg=nstg+1
        sgapb=sgap
        sgap =sgapb-dsgap
        if ( m.eq.mstg ) sgap=0.0
        if ( mstg.gt.1 ) then
          if ( nvbs.ge.2 ) then
            write (6,*) '******** Multistage :',m,'/',mstg
            write (6,*) 'gap% in stress =',sgap/sgapi*100.0d0
          endif
        endif
c
        knr=0
c                                      ---- start of Newton-Raphson loop
        if ( nvbs.ge.3 ) then
          write (6,*)
          write (6,*) '**** start of Newton-Raphson loop'
        endif
c
  100   continue
        knr=knr+1
        nite=nite+1
        if ( nvbs.ge.3 ) then
          write (6,*) '----- NR iteration',knr
          write (6,*) 'inc of p : dp   =',dp
        endif
c
        pt=p+dp
c                                        ---- calc. se and differentials
        do i=1,nttl
          eta(i)=s2(i)-xt2(i)
        enddo
        call jancae_yfunc ( se,dseds,d2seds2,2,
     &                      eta,nttl,nnrm,nshr,
     &                      pryld,ndyld )
c
        if ( nvbs.ge.5 ) then
          text='s2'
          call jancae_print1 ( text,s2,nttl )
          if ( npbs.ne.0 ) then
            text='xt2'
            call jancae_print1 ( text,xt2,nttl )
            text='eta'
            call jancae_print1 ( text,eta,nttl )
          endif
          text='dse/ds'
          call jancae_print1 ( text,dseds,nttl )
          text='d2se/ds2'
          call jancae_print2 ( text,d2seds2,nttl,nttl )
        endif
c                                        ---- calc. sy and differentials
        call jancae_hardencurve ( sy,dsydp,d2sydp2,
     &                            1,pt,prihd,ndihd )

        if ( nvbs.ge.5 ) then
          write (6,*) 'plastic strain p=',pt
          write (6,*) 'flow stress   sy=',sy
          write (6,*) 'hardening dsy/dp=',dsydp
        endif
c                                                          ---- calc. g1
        g1=se-sy-sgap
c                                                          ---- calc. g2
        call jancae_mv ( vv,delast,dseds,nttl,nttl )
        do i=1,nttl
          g2(i)=s2(i)-stry(i)+dp*vv(i)
        enddo
        call jancae_vvs ( g2n,g2,g2,nttl )
        g2n=sqrt(g2n)
c                                                          ---- calc. g3
        if ( npbs.ne.0 ) then
          call jancae_kinematic ( vk,dvkdp,
     &                            dvkds,dvkdx,dvkdxt,
     &                            pt,s2,x2,xt2,
     &                            nttl,nnrm,nshr,
     &                            mxpbs,npbs,
     &                            prkin,ndkin,
     &                            pryld,ndyld )
          do n=1,npbs
            do i=1,nttl
              g3(n,i)=x2(n,i)-x1(n,i)-dp*vk(n,i)
            enddo
          enddo
          g3nn=0.0
          do n=1,npbs
            g3n(n)=0.0
            do i=1,nttl
              g3n(n)=g3n(n)+g3(n,i)*g3(n,i)
            enddo
            g3n(n)=sqrt(g3n(n))
            g3nn=g3nn+g3n(n)*g3n(n)
          enddo
          g3nn=sqrt(g3nn)
        else
          g3nn=0.0
        endif
c
        if ( nvbs.ge.3 ) then
          write (6,*) 'g1 (yield surf) =',g1
          write (6,*) 'g2n (normality) =',g2n
          if ( nvbs.ge.5 ) then
            text='g2 vector'
            call jancae_print1 ( text,g2,nttl )
          endif
          if ( npbs.ne.0 ) then
            if ( nvbs.ge.4 ) then
              do n=1,npbs
                write (6,*) 'g3n(',n,')=',g3n(n)
                if ( nvbs.ge.5 ) then
                  do i=1,nttl
                    uv(i)=g3(n,i)
                  enddo
                  text='g3 vector'
                  call jancae_print1 ( text,uv,nttl )
                endif
              enddo
            endif
          endif
        endif
c                      ---- calc. dependencies common for NR and Dds/Dde
c                                                              * set [A]
        call jancae_setunitm ( am,nnn )
        call jancae_mm ( em,delast,d2seds2,nttl,nttl,nttl )
        do i1=1,npbs+1
          do i2=1,npbs+1
            do j1=1,nttl
              do j2=1,nttl
                k1=(i1-1)*nttl+j1
                k2=(i2-1)*nttl+j2
                if ( i1.eq.1 ) then
                  if ( i2.eq.1 ) then
                    am(k1,k2)=am(k1,k2)+dp*em(j1,j2)
                  else
                    am(k1,k2)=am(k1,k2)-dp*em(j1,j2)
                  endif
                else
                  ip1=i1-1
                  ip2=i2-1
                  if ( i2.eq.1 ) then
                    am(k1,k2)=am(k1,k2)
     &                              -dp*dvkds( ip1,    j1,j2)
                  else
                    am(k1,k2)=am(k1,k2)
     &                              -dp*dvkdx( ip1,ip2,j1,j2)
     &                              -dp*dvkdxt(ip1,    j1,j2)
                  endif
                endif
              enddo
            enddo
          enddo
        enddo
c                                                           ---- set {W}
        call jancae_clear1( wv,nnn )
        do i1=1,npbs+1
          do j1=1,nttl
            k1=(i1-1)*nttl+j1
            if ( i1.eq.1 ) then
              do k2=1,nttl
                wv(k1)=wv(k1)+delast(j1,k2)*dseds(k2)
              enddo
            else
              ip1=i1-1
              wv(k1)=-vk(ip1,j1)-dp*dvkdp(ip1,j1)
            endif
          enddo
        enddo
c                                                      ---- calc. [A]^-1
        call jancae_minv ( ami,am,nnn,det )
c                                                     ---- [C]=[U][A]^-1
        call jancae_mm   ( cm,um,ami,nttl,nnn,nnn )
c
c
c                                                 ---- check convergence
        if ( (abs(g1  /sy).le.tol ).and.
     &       (abs(g2n /sy).le.tol ).and.
     &       (abs(g3nn/sy).le.tol )      ) then
c
          if ( nvbs.ge.2 ) then
            write (6,*) '**** Newton-Raphson converged.',knr
          endif
          dpconv=dp
          do i=1,nttl
            s2conv(i)=s2(i)
            do j=1,npbs
              x2conv(j,i)=x2(j,i)
            enddo
          enddo
          goto 200
        endif
c                                                         ---- solve ddp
c                                                           ---- set {G}
        do i=1,nttl
          gv(i)=g2(i)
        enddo
        do n=1,npbs
          do i=1,nttl
            gv(n*nttl+i)=g3(n,i)
          enddo
        enddo
c                              ---- ddp=(g1-{m}^T[C]{G})/(H+{m}^T[C]{W})
        call jancae_mv   ( vv,cm,gv,nttl,nnn )
        call jancae_vvs  ( top0,dseds,vv,nttl )
        top=g1-top0
        call jancae_mv   ( vv,cm,wv,nttl,nnn )
        call jancae_vvs  ( bot0,dseds,vv,nttl )
        bot=dsydp+bot0
        ddp=top/bot
c                                                         ---- update dp
        dp=dp+ddp
        if ( nvbs.ge.3 ) then
          write (6,*) 'modification of dp:ddp=',ddp
          write (6,*) 'updated             dp=',dp
        endif
        if ( dp.le.0.0 ) then
          if ( nvbs.ge.3 ) then
            write (6,*) 'negative dp is detected.'
            write (6,*) 'multistage is subdivided.'
          endif
          goto 400
        endif
c                                                  ---- update s2 and x2
        do i1=1,npbs+1
          call jancae_clear1( vv,nttl )
          do j1=1,nttl
            k1=(i1-1)*nttl+j1
            do k2=1,nnn
              vv(j1)=vv(j1)-ami(k1,k2)*(gv(k2)+ddp*wv(k2))
            enddo
          enddo
          do j1=1,nttl
            if ( i1.eq.1 ) then
              s2(     j1)=s2(     j1)+vv(j1)
            else
              x2(i1-1,j1)=x2(i1-1,j1)+vv(j1)
            endif
          enddo
        enddo
        call jancae_backsum ( npbs,xt2,x2,nttl,mxpbs )
c
c
        if ( knr.le.maxnr ) goto 100
c                                        ---- end of Newton-Raphson loop
c
  400   continue
        if ( nvbs.ge.2 ) then
          write (6,*) 'Newton Raphson loop is over.',knr
          write (6,*) 'convergence is failed.'
        endif
        if ( nest.lt.maxnest ) then
          nest=nest+1
          newmstg=(mstg-m+1)*ndiv
          goto 300
        else
          write (6,*) 'Nest of multistage is over.',nest
          text='current stress (input)'
          call jancae_print1 ( text,s1,nttl )
          text='strain inc. (input)'
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
        endif
c
  200   continue
c
      enddo
c                                            ---- end of multistage loop
c
c
c                                                 ---- plast.strain inc.
      do i=1,nttl
        dpe(i)=dp*dseds(i)
      enddo
c                                               ---- print out converged
      if ( nvbs.ge.4 ) then
        text='updated stress'
        call jancae_print1 ( text,s2,nttl )
        text='plastic strain inc'
        call jancae_print1 ( text,dpe,nttl )
        if ( npbs.ne.0 ) then
          text='updated part. back stess'
          call jancae_backprint ( text,npbs,x2,nttl,mxpbs )
          text='updated total back stess'
          call jancae_print1 ( text,xt2,nttl )
        endif
      endif
c                                    ---- calc. strain inc. in thickness
      if ( (nttl.eq.3).or.(nttl.eq.5) ) then
        de33=-dpe(1)-dpe(2)
        do i=1,nttl
          de33=de33+d33d(i)*(de(i)-dpe(i))
        enddo
        if ( nvbs.ge.4 ) then
          write (6,*) 'de33=',de33
        endif
      endif
c
      if ( nvbs.ge.1 ) then
        if ( nest.ne.0 ) then
          write (6,*) 'nest of MsRM               :',nest
          write (6,*) 'total no. of stages        :',nstg
          write (6,*) 'total no. of NR iteration  :',nite
          write (6,*) 'initial stress gap         :',sgapi
          write (6,*) 'inc. of equiv.plast.strain :',dp
          write (6,*) 'equiv.plast.strain updated :',p+dp
          write (6,*) 'location ne,ip,lay         :',ne,ip,lay
        endif
      endif
c
      if ( mjac.eq.0 ) then
        do i=1,nttl
          do j=1,nttl
            ddsdde(i,j)=0.0
          enddo
       enddo
       return
      endif
c
      if ( mjac.eq.-1 ) then
        do i=1,nttl
          do j=1,nttl
            ddsdde(i,j)=delast(i,j)
          enddo
       enddo
       return
      endif
c
c                                      ---- consistent material jacobian
c                                                           ---- set [B]
      call jancae_clear2( bm,nnn,nttl )
      i1=1
      i2=1
      do j1=1,nttl
        do j2=1,nttl
          k1=(i1-1)*nttl+j1
          k2=(i2-1)*nttl+j2
          bm(k1,k2)=delast(j1,j2)
        enddo
      enddo
c                                                       ---- [M1]=[N][C]
      call jancae_mm ( em1,d2seds2,cm,nttl,nnn,nttl )
c                                               ---- {V1}={m}-dp*[M1]{W}
      call jancae_mv ( vv ,em1,wv,nttl,nnn )
      do i=1,nttl
        v1(i)=dseds(i)-dp*vv(i)
      enddo
c                                                 ---- [M2]={V1}{m}^T[C]
      call jancae_clear2 ( em2,nttl,nnn )
      do i=1,nttl
        do j=1,nnn
          do k=1,nttl
            em2(i,j)=em2(i,j)+v1(i)*dseds(k)*cm(k,j)
          enddo
        enddo
      enddo
c                                                  ---- S1=H+{m}^T[C]{W}
      sc1=dsydp
      do i=1,nttl
        do j=1,nnn
          sc1=sc1+dseds(i)*cm(i,j)*wv(j)
        enddo
      enddo
c                                     ---- [M3]=[I]-[dp*[M1]-[M2]/S1][B]
      call jancae_setunitm ( em3,nttl )
      do i=1,nttl
        do j=1,nttl
          do k=1,nnn
            em3(i,j)=em3(i,j)
     &                    -(dp*em1(i,k)+em2(i,k)/sc1)*bm(k,j)
          enddo
        enddo
      enddo
c                                                     ---- [Dc]=[De][M3]
      call jancae_mm ( ddsdde,delast,em3,nttl,nttl,nttl )
c
c                                                    ---- check symmetry
      nsym=0
      d=0.0d0
      a=0.0d0
      do i=1,nttl
        do j=i,nttl
          dd=     ddsdde(i,j)-ddsdde(j,i)
          aa=0.5*(ddsdde(i,j)+ddsdde(i,j))
          d=d+dd*dd
          a=a+aa*aa
        enddo
      enddo
      a=sqrt(d/a)
      if ( a.gt.1.0d-8 ) then
        if ( nvbs.ge.4 ) then
          write (6,*) 'ddsdde is not symmetric.',a
          text='material jacobian (nonsym)'
          call jancae_print2 ( text,ddsdde,nttl,nttl )
        endif
c                                                    ---- symmetrization
        if ( nsym.eq.1 ) then
          do i=1,nttl
            do j=i+1,nttl
              aaa=0.5d0*(ddsdde(i,j)+ddsdde(j,i))
              ddsdde(i,j)=aaa
              ddsdde(j,i)=aaa
            enddo
          enddo
        endif
      endif
c
      if ( nvbs.ge.4 ) then
        text='material jacobian (output)'
        call jancae_print2 ( text,ddsdde,nttl,nttl )
      endif
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
c      nvbs0=0   ! verbose mode
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
      nechk=1    ! element no. to be checked
      ipchk=1    ! integration point no. to checked
      laychk=1   ! layer no. to be checked
c
      nvbs=0
      nchk=nechk*ipchk*laychk
      if ( nchk.gt.0 ) then
        if ( (ne .eq.nechk ).and.
     &       (ip .eq.ipchk ).and.
     &       (lay.eq.laychk)      ) then
          nvbs=nvbs0
        endif
      endif
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
      ntela=nint(prela(1))
      select case ( ntela )
c
      case ( 0:1 )    !  isotropic linear elasticity (Hooke)
c
        if ( ntela.eq.0 ) then
          eyoung=prela(2)                   ! Young's modulus
          epoas =prela(3)                   ! Poisson's ratio
          erigid=eyoung/2.0d0/(1.0d0+epoas) ! Rigidity
        else
          ek=prela(2)                       ! Bulk modulus
          eg=prela(3)                       ! Rigidity
          eyoung=9.0d0*ek*eg/(3.0d0*ek+eg)  ! Young's modulus
          epoas =(eyoung-2.0d0*eg)/2.0d0/eg ! Poisson's ratio
          erigid=eg
        endif
c                                       ---- set 6*6 matrix for 3d solid
        call jancae_clear2( delast3d,6,6 )
        do i=1,3
          do j=1,3
            if ( i.eq.j ) then
              delast3d(i,j)=1.0d0-epoas
            else
              delast3d(i,j)=epoas
            endif
          enddo
        enddo
        do i=4,6
          delast3d(i,i)=0.5d0-epoas
        enddo
        coef=erigid/(0.5d0-epoas)
        do i=1,6
          do j=1,6
            delast3d(i,j)=coef*delast3d(i,j)
          enddo
        enddo
c
      case default  !  error
        write (6,*) 'elasticity code error in jancae_setelast'
        write (6,*) 'ntela=',ntela
        call jancae_exit ( 9000 )
c
      end select
c
c                                       ---- condensation for 2D problem
      do ib=1,2
        if (ib.eq.1) then
          ni=nnrm
        else
          ni=nshr
        endif
        do jb=1,2
          if (jb.eq.1) then
            nj=nnrm
          else
            nj=nshr
          endif
          do is=1,ni
            i =(ib-1)*nnrm+is
            i3=(ib-1)*3   +is
            do js=1,nj
              j =(jb-1)*nnrm+js
              j3=(jb-1)*3   +js
              delast(i,j)=delast3d(i3,j3)
            enddo
          enddo
        enddo
      enddo
c                                     ---- plane stress or shell element
      if ( nnrm.eq.2 ) then
        d33=delast3d(3,3)
        do ib=1,2
          if (ib.eq.1) then
            ni=nnrm
          else
            ni=nshr
          endif
          do jb=1,2
            if (jb.eq.1) then
              nj=nnrm
            else
              nj=nshr
            endif
            do is=1,ni
              i =(ib-1)*nnrm+is
              i3=(ib-1)*3   +is
              do js=1,nj
                j =(jb-1)*nnrm+js
                j3=(jb-1)*3   +js
                delast(i,j)=delast(i,j)
     &                     -delast3d(i3,3)*delast3d(3,j3)/d33
              enddo
            enddo
          enddo
        enddo
c                         ---- elastic strain in thickness direction e_t
c                                             ---- e_t=SUM(d33d(i)*e(i))
        do i=1,nttl
          if ( i.le.nnrm ) then
            id=i
          else
            id=i-nnrm+3
          endif
          d33d(i)=-delast3d(3,id)/d33
        enddo
      endif
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     print parameters for elastic info
c
      subroutine jancae_elast_print ( prela,ndela )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension prela(ndela)
c
      ntela=nint(prela(1))
      write(6,*)
      write (6,*) '*** Elastic Properties',ntela
      select case ( ntela )
      case ( 0 )
        write (6,*) 'Hooke Isotropic Elasticity'
        write (6,*) 'Youngs modulus=',prela(1+1)
        write (6,*) 'Poissons ratio=',prela(1+2)
      case ( 1 )
        write (6,*) 'Hooke Isotropic Elasticity'
        write (6,*) 'Bulk modulus  =',prela(1+1)
        write (6,*) 'Shear modulus =',prela(1+2)
      end select
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
      if ( npbs.eq.0 ) then
        isvtnsr=nttl
      else
        isvtnsr=nttl*(1+npbs)
      endif
      isvttl=isvrsvd+isvsclr+isvtnsr
      if ( nisv.lt.isvttl ) then
        write (6,*) 'check number of internal state variables (isv)'
        write (6,*) 'nisv must be larger than',isvttl
        write (6,*) 'nisv=',nisv
        write (6,*) 'isv : required       ',isvttl
        write (6,*) 'isv : system reserved',isvrsvd
        write (6,*) 'isv : for scaler     ',isvsclr
        write (6,*) 'isv : for tensor     ',isvtnsr
        call jancae_exit ( 9000 )
      endif
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
      p=stv(isvrsvd+1)
c                                              ---- plastic strain comp.
      do i=1,nttl
        pe(i)=stv(isvrsvd+isvsclr+i)
      enddo
c                                         ---- partial back stress comp.
      if ( npbs.ne.0 ) then
        do nb=1,npbs
          do i=1,nttl
            it=isvrsvd+isvsclr+nttl*nb+i
            x(nb,i)=stv(it)
          enddo
        enddo
      endif
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
      do i=1,nttl
        xt(i)=0.0
      enddo
      if ( npbs.eq.0 ) return
c
      do i=1,nttl
        do j=1,npbs
          xt(i)=xt(i)+x(j,i)
        enddo
      enddo
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
      if ( npbs.eq.0 ) return
c
      do i=1,nttl
        do j=1,npbs
          xx(j,i)=x(j,i)
        enddo
      enddo
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
     &                             npbs )
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      dimension prop(mxprop)
c
      n=0
      p=prop(n+1)
      if ( p.ge.1000.0d0 ) p=p-1000.d0
      nela=nint(p)
      select case (nela)
        case (0) ; nd= 2
        case (1) ; nd= 2     ! ht180110
        case default
          write (6,*) 'error elastic property id :',nela
          call jancae_exit ( 9000 )
      end select
      ndela=nd+1
c
      n=ndela
      nyld=nint(prop(n+1))
      select case (nyld)
        case ( 0) ; nd=  0   ! Mises
        case ( 1) ; nd=  6   ! Hill48
        case ( 2) ; nd= 19   ! Yld2004
        case ( 3) ; nd= 14   ! Cazacu2006
        case ( 4) ; nd=  8   ! Karafillis-Boyce
        case ( 5) ; nd= 10   ! Hu2005
        case ( 6) ; nd= 16   ! Yoshida2011
        case (-1) ; nd=  9   ! Gotoh bi-quad
        case (-2) ; nd=  9   ! Yld2000-2d
        case (-3)           ! Vegter
          nd= 3+4*nint(prop(n+2))
        case (-4) ; nd= 9   ! BBC2005
        case (-5) ; nd= 4   ! Yld89
        case (-6)           ! BBC2008
          nd= 2+8*nint(prop(n+2))
        case default
          write (6,*) 'error yield function id :',nyld
          call jancae_exit ( 9000 )
      end select
      ndyld=nd+1
c
      n=ndela+ndyld
      nihd=nint(prop(n+1))
      select case (nihd)
        case ( 0) ; nd= 1   ! Perfecty Plastic
        case ( 1) ; nd= 2   ! Linear
        case ( 2) ; nd= 3   ! Swift
        case ( 3) ; nd= 3   ! Ludwick
        case ( 4) ; nd= 3   ! Voce
        case ( 5) ; nd= 4   ! Voce + Linear
        case ( 6) ; nd= 7   ! Voce + Swift
        case default
          write (6,*) 'error work hardening curve id :',nihd
          call jancae_exit ( 9000 )
      end select
      ndihd=nd+1
c
      n=ndela+ndyld+ndihd
      nkin=nint(prop(n+1))
      select case (nkin)
        case ( 0) ; nd= 0 ; npbs= 0   ! No Kinematic Hardening
        case ( 1) ; nd= 1 ; npbs= 1   ! Prager
        case ( 2) ; nd= 1 ; npbs= 1   ! Ziegler
        case ( 3) ; nd= 2 ; npbs= 1   ! Armstrong & Frederick (1966)
        case ( 4)                     ! Chaboche (1979)
          npbs= nint(prop(n+2))
          nd= 2*npbs + 1
        case ( 5)                     ! Chaboche (1979) - Ziegler Model
          npbs= nint(prop(n+2))
          nd= 2*npbs + 1
        case ( 6) ; nd= 5 ; npbs= 2   ! Yoshida-Uemori
        case default
          write (6,*) 'error kinematic hardening id :',nkin
          call jancae_exit ( 9000 )
      end select
      ndkin=nd+1
c
      return
      end
c
c
c
c***********************************************************************
c     JANCAE/UMMDp : Isotropic Hardening
c***********************************************************************
c
c      0 : Perfectly Plastic
c
c      1 : Linear
c      2 : Swift
c      3 : Ludwick
c      4 : Voce
c      5 : Voce + Linear
c      6 : Voce + Swift
c
c-----------------------------------------------------------------------
c     hardening curve
c
      subroutine jancae_hardencurve ( sy,dsydp,d2sydp2,
     &                                nreq,p,prihd,ndihd )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension prihd(ndihd)
c
      ntihd=nint(prihd(1))
      select case ( ntihd )
c
      case ( 0 )                                     ! Perfectly Plastic
        sy=prihd(1+1)
        if ( nreq.ge.1 ) then
          dsydp=0.0
          if ( nreq.ge.2 ) then
            d2sydp2=0.0
          endif
        endif
c
      case ( 1 )                                                ! Linear
        sy0 =prihd(1+1)
        hard=prihd(1+2)
        sy=sy0+hard*p
        if ( nreq.ge.1 ) then
          dsydp=hard
          if ( nreq.ge.2 ) then
            d2sydp2=0.0
          endif
        endif
c
      case ( 2 )                                                 ! Swift
        c =prihd(1+1)
        e0=prihd(1+2)
        en=prihd(1+3)
        sy=c*(e0+p)**en
        if ( nreq.ge.1 ) then
          dsydp=en*c*(e0+p)**(en-1.0d0)
          if ( nreq.ge.2 ) then
            d2sydp2=en*c*(en-1.0d0)*(e0+p)**(en-2.0d0)
          endif
        endif
c
      case ( 3 )                                               ! Ludwick
        sy0=prihd(1+1)
        c  =prihd(1+2)
        en =prihd(1+3)
        sy=sy0+c*p**en
        if ( nreq.ge.1 ) then
          dsydp=en*c*p**(en-1.0d0)
          if ( nreq.ge.2 ) then
            d2sydp2=en*c*(en-1.0d0)*p**(en-2.0d0)
          endif
        endif
c
      case ( 4 )                                                  ! Voce
        sy0=prihd(1+1)
        q  =prihd(1+2)
        b  =prihd(1+3)
        sy=sy0+q*(1.0d0-exp(-b*p))
        if ( nreq.ge.1 ) then
          dsydp=q*b*exp(-b*p)
          if ( nreq.ge.2 ) then
            d2sydp2=-q*b*b*exp(-b*p)
          endif
        endif
c
      case ( 5 )                                         ! Voce + Linear
        sy0=prihd(1+1)
        q  =prihd(1+2)
        b  =prihd(1+3)
        c  =prihd(1+4)
        sy=sy0+q*(1.0d0-exp(-b*p))+c*p
        if ( nreq.ge.1 ) then
          dsydp=q*b*exp(-b*p)+c
          if ( nreq.ge.2 ) then
            d2sydp2=-q*b*b*exp(-b*p)
          endif
        endif
c
      case ( 6 )                                          ! Voce + Swift
        a  =prihd(1+1)
        sy0=prihd(1+2)
        q  =prihd(1+3)
        b  =prihd(1+4)
        c  =prihd(1+5)
        e0 =prihd(1+6)
        en =prihd(1+7)
        sy=       a  * ( sy0+q*(1.0d0-exp(-b*p)) )+
     &     (1.0d0-a) * ( c*(e0+p)**en            )
        if ( nreq.ge.1 ) then
          dsydp=       a  * ( q*b*exp(-b*p)           )+
     &          (1.0d0-a) * ( en*c*(e0+p)**(en-1.0d0) )
          if ( nreq.ge.2 ) then
            d2sydp2=
     &           a  * ( -q*b*b*exp(-b*p)                   )+
     &    (1.0d0-a) * ( en*c*(en-1.0d0)*(e0+p)**(en-2.0d0) )
          endif
        endif
c
      case default
        write (6,*) 'hardening type error',ntihd
        call jancae_exit (9000)
c
      end select
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     print parameters for isotropic hardening info
c
      subroutine jancae_harden_print ( prihd,ndihd )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension prihd(ndihd)
c
      ntihd=nint(prihd(1))
      write (6,*)
      write (6,*) '*** Isotropic Hardening Law',ntihd
      select case ( ntihd )
      case ( 0 )
        write (6,*) 'Perfect Plasticity'
        write (6,*) 'sy_const=',prihd(1+1)
      case ( 1 )
        write (6,*) 'Linear'
        write (6,*) 'sy = sy0+h*p'
        write (6,*) 'sy0=',prihd(1+1)
        write (6,*) 'h  =',prihd(1+2)
      case ( 2 )
        write (6,*) 'Swift'
        write (6,*) 'sy = c*(e0+p)^en'
        write (6,*) 'c =',prihd(1+1)
        write (6,*) 'e0=',prihd(1+2)
        write (6,*) 'en=',prihd(1+3)
      case ( 3 )
        write (6,*) 'Ludwick'
        write (6,*) 'sy = sy0+c*p^en'
        write (6,*) 'sy0=',prihd(1+1)
        write (6,*) 'c  =',prihd(1+2)
        write (6,*) 'en =',prihd(1+3)
      case ( 4 )
        write (6,*) 'Voce '
        write (6,*) 'sy = sy0+q*(1-exp(-b*p))'
        write (6,*) 'sy0=',prihd(1+1)
        write (6,*) 'q  =',prihd(1+2)
        write (6,*) 'b  =',prihd(1+3)
      case ( 5 )
        write (6,*) 'Voce + Linear'
        write (6,*) 'sy = sy0+q*(1-exp(-b*p))+c*p'
        write (6,*) 'sy0=',prihd(1+1)
        write (6,*) 'q  =',prihd(1+2)
        write (6,*) 'b  =',prihd(1+3)
        write (6,*) 'c  =',prihd(1+4)
      case ( 6 )
        write (6,*) 'Voce+Swift a *( sy0+q*(1-exp(-b*p)) )+'
        write (6,*) 'sy = a *( sy0+q*(1-exp(-b*p)) )+ (1-a)*(
     1 c*(e0+p)^en)'
        write (6,*) 'a  =',prihd(1+1)
        write (6,*) 'sy0=',prihd(1+2)
        write (6,*) 'q  =',prihd(1+3)
        write (6,*) 'b  =',prihd(1+4)
        write (6,*) 'c  =',prihd(1+5)
        write (6,*) 'e0 =',prihd(1+6)
        write (6,*) 'en =',prihd(1+7)
      end select
c
      return
      end
c
c
c
c***********************************************************************
c     JANCAE/UMMDp : Kinematic Hardening
c***********************************************************************
c
c      0 : No Kinematic Hardening
c
c      1 : Prager
c      2 : Ziegler
c      3 : Armstrong & Frederick (1966)
c      4 : Chaboche (1979)
c      5 : Chaboche (1979) - Ziegler Model
c      6 : Yoshida-Uemori
c
c-----------------------------------------------------------------------
c     calc. kinematic hardening function
c
      subroutine jancae_kinematic ( vk,dvkdp,
     &                              dvkds,dvkdx,dvkdxt,
     &                              p,s,x,xt,
     &                              nttl,nnrm,nshr,
     &                              mxpbs,npbs,
     &                              prkin,ndkin,
     &                              pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension vk(    npbs,nttl),
     &          dvkdp( npbs,nttl),
     &          dvkdx( npbs,npbs,nttl,nttl),
     &          dvkds( npbs,nttl,nttl),
     &          dvkdxt(npbs,nttl,nttl),
     &          s(nttl),x(mxpbs,nttl),xt(nttl),
     &          prkin(ndkin),pryld(ndyld)
c
c
      ntkin=nint(prkin(1))
c                                                        ---- initialize
      do i=1,npbs
        do j=1,nttl
          vk(  i,j)=0.0
          dvkdp(i,j)=0.0
          do k=1,nttl
            dvkds( i,j,k)=0.0
            dvkdxt(i,j,k)=0.0
            do l=1,npbs
              dvkdx( i,l,j,k)=0.0
            enddo
          enddo
        enddo
      enddo
c
      select case ( ntkin )
c
      case ( 0 )                                ! No Kinematic Hardening
        return
c
      case ( 1 )                                                ! Prager
        call jancae_kin_prager    ( vk,dvkdp,
     &                              dvkds,dvkdx,dvkdxt,
     &                              p,s,x,xt,
     &                              nttl,nnrm,nshr,
     &                              mxpbs,npbs,
     &                              prkin,ndkin,pryld,ndyld )
c
      case ( 2 )                                               ! Ziegler
        call jancae_kin_ziegler   ( vk,dvkdp,
     &                              dvkds,dvkdx,dvkdxt,
     &                              p,s,x,xt,
     &                              nttl,nnrm,nshr,
     &                              mxpbs,npbs,
     &                              prkin,ndkin,pryld,ndyld )
c
      case ( 3 )                          ! Armstrong & Frederick (1966)
        call jancae_kin_armstrong ( vk,dvkdp,
     &                              dvkds,dvkdx,dvkdxt,
     &                              p,s,x,xt,
     &                              nttl,nnrm,nshr,
     &                              mxpbs,npbs,
     &                              prkin,ndkin,pryld,ndyld )
c
      case ( 4 )                                       ! Chaboche (1979)
        call jancae_kin_chaboche1979 ( vk,dvkdp,
     &                                 dvkds,dvkdx,dvkdxt,
     &                                 p,s,x,xt,
     &                                 nttl,nnrm,nshr,
     &                                 mxpbs,npbs,
     &                                 prkin,ndkin,pryld,ndyld )
c
      case ( 5 )                       ! Chaboche (1979) - Ziegler Model
        call jancae_kin_chaboche1979_ziegler ( vk,dvkdp,
     &                                         dvkds,dvkdx,dvkdxt,
     &                                         p,s,x,xt,
     &                                         nttl,nnrm,nshr,
     &                                         mxpbs,npbs,
     &                                         prkin,ndkin,pryld,ndyld )
c
      case ( 6 )                                        ! Yoshida-Uemori
        call jancae_kin_yoshida_uemori ( vk,dvkdp,
     &                                   dvkds,dvkdx,dvkdxt,
     &                                   p,s,x,xt,
     &                                   nttl,nnrm,nshr,
     &                                   mxpbs,npbs,
     &                                   prkin,ndkin,pryld,ndyld )
c
      case default
        write (6,*) 'still not be supported. ntkin=',ntkin
        call jancae_exit ( 9000 )
      end select
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     dseds and d2seds2 for kinematic hardening
c
      subroutine jancae_dseds_kin ( eta,seta,dseds,d2seds2,
     &                              nttl,nnrm,nshr,
     &                              pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension eta(nttl),
     &          dseds(nttl),d2seds2(nttl,nttl),
     &          pryld(ndyld)
c
c
c                         ---- dseds and d2seds2 for plastic strain inc.
      call jancae_yfunc  ( seta,dseds,d2seds2,2,
     &                     eta,nttl,nnrm,nshr,
     &                     pryld,ndyld )
c
c                   ---- engineering shear strain -> tensor shear strain
      do i=nnrm+1,nttl
        dseds(i)=0.5d0*dseds(i)
        do j=1,nttl
          d2seds2(i,j)=0.5d0*d2seds2(i,j)
        enddo
      enddo
c                                          ---- for plane stress problem
      if ( nnrm.eq.2 ) then
        em1 =        dseds(1)
        em2 =        dseds(2)
        dseds(1)=    dseds(1)    +em1 +em2
        dseds(2)=    dseds(2)    +em2 +em1
        en11=        d2seds2(1,1)
        en12=        d2seds2(1,2)
        en13=        d2seds2(1,3)
        en21=        d2seds2(2,1)
        en22=        d2seds2(2,2)
        en23=        d2seds2(2,3)
        d2seds2(1,1)=d2seds2(1,1)+en11+en21
        d2seds2(1,2)=d2seds2(1,2)+en12+en22
        d2seds2(1,3)=d2seds2(1,3)+en13+en23
        d2seds2(2,1)=d2seds2(2,1)+en21+en11
        d2seds2(2,2)=d2seds2(2,2)+en22+en12
        d2seds2(2,3)=d2seds2(2,3)+en23+en13
      endif
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     print parameters for kinematic hardening
c
      subroutine jancae_kinematic_print ( prkin,ndkin,npbs )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension prkin(ndkin)
c
      ntkin=nint(prkin(1))
      write (6,*)
      write (6,*) '*** Kinematic Hardening Law',ntkin
      select case ( ntkin )
      case ( 0 )
        write (6,*) 'No Kinematic Hardening'
      case ( 1 )
        write (6,*) 'Prager dX=(2/3)*c*{dpe}'
        write (6,*) 'c =',prkin(1+1)
      case ( 2 )
        write (6,*) 'Ziegler dX=dp*c*{{s}-{X}}'
        write (6,*) 'c =',prkin(1+1)
      case ( 3 )
        write (6,*) 'Armstrong-Frederick (1966)'
        write (6,*) 'dX=(2/3)*c*{dpe}-dp*g*{X}'
        write (6,*) 'c =',prkin(1+1)
        write (6,*) 'g =',prkin(1+2)
      case ( 4 )
        write (6,*) 'Chaboche (1979)'
        write (6,*) 'dx(j)=c(j)*(2/3)*{dpe}-dp*g(j)*{x(j)}'
        write (6,*) 'no. of x(j) =',npbs
        do i=1,npbs
          n0=1+(i-1)*2
          write (6,*) 'c(',i,')=',prkin(1+n0+1)
          write (6,*) 'g(',i,')=',prkin(1+n0+2)
        enddo
      case ( 5 )
        write (6,*) 'Chaboche (1979) - Ziegler Model'
        write (6,*) 'dx(j)=((c(j)/se)*{{s}-{X}}-g(j)*{x(j)})*dp'
        write (6,*) 'no. of x(j) =',npbs
        do i=1,npbs
          n0=1+(i-1)*2
          write (6,*) 'c(',i,')=',prkin(1+n0+1)
          write (6,*) 'g(',i,')=',prkin(1+n0+2)
        enddo
      case ( 6 )
        write (6,*) 'Yoshida-Uemori'
        write (6,*) 'no. of x(j) =',npbs
        write (6,*) 'C=',prkin(1+1)
        write (6,*) 'Y=',prkin(1+2)
        write (6,*) 'a=',prkin(1+3)
        write (6,*) 'k=',prkin(1+4)
        write (6,*) 'b=',prkin(1+5)
      end select
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     Prager
c
      subroutine jancae_kin_prager ( vk,dvkdp,
     &                               dvkds,dvkdx,dvkdxt,
     &                               p,s,x,xt,
     &                               nttl,nnrm,nshr,
     &                               mxpbs,npbs,
     &                               prkin,ndkin,
     &                               pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension vk(    npbs,nttl),
     &          dvkdp( npbs,nttl),
     &          dvkdx( npbs,npbs,nttl,nttl),
     &          dvkds( npbs,nttl,nttl),
     &          dvkdxt(npbs,nttl,nttl),
     &          s(nttl),x(mxpbs,nttl),xt(nttl),
     &          prkin(ndkin),pryld(ndyld)
c
      dimension eta(nttl),dseds(nttl),d2seds2(nttl,nttl)
c
      c=prkin(2)/3.0d0*2.0d0
c
      do i=1,nttl
        eta(i)=s(i)-xt(i)
      enddo
c
      call jancae_dseds_kin ( eta,seta,dseds,d2seds2,
     &                        nttl,nnrm,nshr,
     &                        pryld,ndyld )
c
      n=1
c
      do i=1,nttl
        vk(n,i)=c*dseds(i)
      enddo
c
      dcdp=0.0d0
      do i=1,nttl
        dvkdp(n,i)=dcdp*dseds(i)
      enddo
c
      do i=1,nttl
        do j=1,nttl
          dvkds( n,  i,j)= c*d2seds2(i,j)
          dvkdx( n,n,i,j)=-c*d2seds2(i,j)
          dvkdxt(n,  i,j)= 0.0
        enddo
      enddo
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     Ziegler
c
      subroutine jancae_kin_ziegler ( vk,dvkdp,
     &                                dvkds,dvkdx,dvkdxt,
     &                                p,s,x,xt,
     &                                nttl,nnrm,nshr,
     &                                mxpbs,npbs,
     &                                prkin,ndkin,
     &                                pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension vk(    npbs,nttl),
     &          dvkdp( npbs,nttl),
     &          dvkdx( npbs,npbs,nttl,nttl),
     &          dvkds( npbs,nttl,nttl),
     &          dvkdxt(npbs,nttl,nttl),
     &          s(nttl),x(mxpbs,nttl),xt(nttl),
     &          prkin(ndkin),pryld(ndyld)
c
      dimension eta(nttl),am(nttl,nttl)
c
      c=prkin(2)
c
      do i=1,nttl
        eta(i)=s(i)-xt(i)
      enddo
c
      n=1
      do i=1,nttl
        vk(n,i)=c*eta(i)
      enddo
c
      dcdp=0.0
      do i=1,nttl
        dvkdp(n,i)=dcdp*eta(i)
      enddo
c
      call jancae_setunitm ( am,nttl )
      do i=1,nttl
        do j=1,nttl
          dvkds( n,  i,j)= c*am(i,j)
          dvkdx( n,n,i,j)=-c*am(i,j)
          dvkdxt(n,  i,j)= 0.0
        enddo
      enddo
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     Armstrong-Frederick (1966)
c
      subroutine jancae_kin_armstrong ( vk,dvkdp,
     &                                  dvkds,dvkdx,dvkdxt,
     &                                  p,s,x,xt,
     &                                  nttl,nnrm,nshr,
     &                                  mxpbs,npbs,
     &                                  prkin,ndkin,
     &                                  pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension vk(    npbs,nttl),
     &          dvkdp( npbs,nttl),
     &          dvkdx( npbs,npbs,nttl,nttl),
     &          dvkds( npbs,nttl,nttl),
     &          dvkdxt(npbs,nttl,nttl),
     &          s(nttl),x(mxpbs,nttl),xt(nttl),
     &          prkin(ndkin),pryld(ndyld)
c
      dimension eta(nttl),dseds(nttl),d2seds2(nttl,nttl),
     &          am(nttl,nttl)
c
      c=prkin(1+1)/3.0d0*2.0d0
      g=prkin(1+2)
c
      do i=1,nttl
        eta(i)=s(i)-xt(i)
      enddo
c
      call jancae_dseds_kin ( eta,seta,dseds,d2seds2,
     &                        nttl,nnrm,nshr,
     &                        pryld,ndyld )
c
c
      n=1
      do i=1,nttl
        vk(n,i)=c*dseds(i)-g*xt(i)
      enddo
c
      dcdp=0.0d0
      dgdp=0.0d0
      do i=1,nttl
        dvkdp(n,i)=dcdp*dseds(i)-dgdp*xt(i)
      enddo
c
      call jancae_setunitm ( am,nttl )
      do i=1,nttl
        do j=1,nttl
          dvkds( n,  i,j)= c*d2seds2(i,j)
          dvkdx( n,n,i,j)=-c*d2seds2(i,j)-g*am(i,j)
          dvkdxt(n,  i,j)= 0.0
        enddo
      enddo
c
      return
      end
c
c
c-----------------------------------------------------------------------
c     Chaboche (1979)
c
      subroutine jancae_kin_chaboche1979 ( vk,dvkdp,
     &                                     dvkds,dvkdx,dvkdxt,
     &                                     p,s,x,xt,
     &                                     nttl,nnrm,nshr,
     &                                     mxpbs,npbs,
     &                                     prkin,ndkin,
     &                                     pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension vk(    npbs,nttl),
     &          dvkdp( npbs,nttl),
     &          dvkdx( npbs,npbs,nttl,nttl),
     &          dvkds( npbs,nttl,nttl),
     &          dvkdxt(npbs,nttl,nttl),
     &          s(nttl),x(mxpbs,nttl),xt(nttl),
     &          prkin(ndkin),pryld(ndyld)
c
      dimension eta(nttl),dseds(nttl),d2seds2(nttl,nttl),
     &          am(nttl,nttl)
c
      do i=1,nttl
        eta(i)=s(i)-xt(i)
      enddo
c
      call jancae_dseds_kin ( eta,seta,dseds,d2seds2,
     &                        nttl,nnrm,nshr,
     &                        pryld,ndyld )
c
      call jancae_setunitm ( am,nttl )
      do n=1,npbs
        n0=1+(n-1)*2
        c=prkin(1+n0+1)/3.0d0*2.0d0
        g=prkin(1+n0+2)
        do i=1,nttl
          vk(n,i)=c*dseds(i)-g*x(n,i)
        enddo
        dcdp=0.0d0
        dgdp=0.0d0
        do i=1,nttl
          dvkdp(n,i)=dcdp*dseds(i)-dgdp*x(n,i)
        enddo
        do i=1,nttl
          do j=1,nttl
            dvkds( n,  i,j)= c*d2seds2(i,j)
            dvkdx( n,n,i,j)=-g*am(i,j)
            dvkdxt(n,  i,j)=-c*d2seds2(i,j)
          enddo
        enddo
      enddo
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     Chaboche (1979) - Ziegler Model
c
      subroutine jancae_kin_chaboche1979_ziegler ( vk,dvkdp,
     &                                     dvkds,dvkdx,dvkdxt,
     &                                     p,s,x,xt,
     &                                     nttl,nnrm,nshr,
     &                                     mxpbs,npbs,
     &                                     prkin,ndkin,
     &                                     pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension vk(    npbs,nttl),
     &          dvkdp( npbs,nttl),
     &          dvkdx( npbs,npbs,nttl,nttl),
     &          dvkds( npbs,nttl,nttl),
     &          dvkdxt(npbs,nttl,nttl),
     &          s(nttl),x(mxpbs,nttl),xt(nttl),
     &          prkin(ndkin),pryld(ndyld)
c
      dimension eta(nttl),dseds(nttl),d2seds2(nttl,nttl),
     &          am(nttl,nttl)
c
      do i=1,nttl
        eta(i)=s(i)-xt(i)
      enddo
c
      call jancae_dseds_kin ( eta,seta,dseds,d2seds2,
     &                        nttl,nnrm,nshr,
     &                        pryld,ndyld )
c
      call jancae_setunitm ( am,nttl )
      do n=1,npbs
        n0=1+(n-1)*2
        c=prkin(1+n0+1)
        g=prkin(1+n0+2)
        do i=1,nttl
          vk(n,i)=(c/seta)*eta(i)-g*x(n,i)
        enddo
        dcdp=0.0d0
        dgdp=0.0d0
        do i=1,nttl
          dvkdp(n,i)=(dcdp/seta)*eta(i)-dgdp*x(n,i)
        enddo
        do i=1,nttl
          do j=1,nttl
            dvkds( n,  i,j)= (c/seta)*am(i,j)
            dvkdx( n,n,i,j)=-g*am(i,j)
            dvkdxt(n,  i,j)=-(c/seta)*am(i,j)
          enddo
        enddo
      enddo
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     Yoshida_Uemori (***)
c
      subroutine jancae_kin_yoshida_uemori ( vk,dvkdp,
     &                                     dvkds,dvkdx,dvkdxt,
     &                                     p,s,x,xt,
     &                                     nttl,nnrm,nshr,
     &                                     mxpbs,npbs,
     &                                     prkin,ndkin,
     &                                     pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension vk(    npbs,nttl),
     &          dvkdp( npbs,nttl),
     &          dvkdx( npbs,npbs,nttl,nttl),
     &          dvkds( npbs,nttl,nttl),
     &          dvkdxt(npbs,nttl,nttl),
     &          s(nttl),x(mxpbs,nttl),xt(nttl),
     &          prkin(ndkin),pryld(ndyld)
c
      dimension eta(nttl),dseds(nttl),d2seds2(nttl,nttl),
     &          am(nttl,nttl)
c
      pc=prkin(1+1)  ! C
      py=prkin(1+2)  ! Y
      pa=prkin(1+3)  ! a
      pk=prkin(1+4)  ! k
      pb=prkin(1+5)  ! b
c
      do i=1,nttl
        eta(i)=s(i)-xt(i)
      enddo
c
      call jancae_dseds_kin ( eta,seta,dseds,d2seds2,
     &                        nttl,nnrm,nshr,
     &                        pryld,ndyld )
      call jancae_setunitm ( am,nttl )
c
      n=1
      do i=1,nttl
        vk(n,i)=pc*( (pa/py)*eta(i) - sqrt(pa/seta)*x(n,i) )
      enddo
      do i=1,nttl
        dvkdp(n,i)=0.0
      enddo
      do i=1,nttl
        do j=1,nttl
          dvkds( n,  i,j)= pc*pa/py*am(i,j)
          dvkdxt(n,  i,j)=-pc*pa/py*am(i,j)
          dvkdx( n,n,i,j)= pc*sqrt(pa)*
     &                    ( -am(i,j)/sqrt(seta)
     &                      +x(n,i)*dseds(j)/(2.0d0*seta**(1.5d0)) )
        enddo
      enddo
c
      n=2
      do i=1,nttl
        vk(n,i)=pk*( 2.0d0/3.0d0*pb*dseds(i)-x(n,i) )
      enddo
      do i=1,nttl
        dvkdp(n,i)=0.0
      enddo
      do i=1,nttl
        do j=1,nttl
          dvkds( n,  i,j)= 2.0d0/3.0d0*pb+pk*d2seds2(i,j)
          dvkdxt(n,  i,j)=-2.0d0/3.0d0*pb+pk*d2seds2(i,j)
          dvkdx( n,n,i,j)= -pk*am(i,j)
        enddo
      enddo
c
      return
      end
c
c
c
c***********************************************************************
c     JANCAE/MMSM/MPWG : UTILITY SUBROUTINES for UMMDp
c***********************************************************************
c     jancae_clear1 ( a,n )
c       clear 1st order tensor (vector) a(n)
c     jancae_clear2 ( a,n,m )
c       clear 2nd order tensor (matrix) a(n,m)
c     jancae_clear3 ( a,n,m,l )
c       clear 3rd order tensor a(n,m,l)
c     jancae_setunitm ( a,n )
c       set unit 2nd oder tensor [I]
c     jancae_print1 ( text,a,n )
c       print 1st order tensor with title (text)
c     jancae_print2 ( text,a,n,m )
c       print 2nd order tensor with title (text)
c     jancae_mv (v,a,u,nv,nu)
c       mutiply matrix to vector v(nv)=a(nv,nu)*u(nu)
c     jancae_mm (a,b,c,na1,na2,nbc)
c       mutiply matrix and matrix a(na1,na2)=b(na1,nbc)*c(nbc,na2)
c     jancae_vvs ( s,u,v,n )
c       calc. scalar (inner) product of v(n) & u(n)
c     jancae_minv ( b,a,n,d )
c       calc. inverse matrix b(n,n)=a(n,n)^-1 and det(a)
c       branch to following routines
c       jancae_ludcmp( a,n,indx,d ) : LU-decomposition
c       jancae_lubksb(a,n,indx,b)   : backward subsitution
c       jancae_minv2 ( b,a,deta )   : for 2*2 matrix
c       jancae_minv3 ( b,a,deta )   : for 3*3 matrix
c     jancae_eigen_sym3 ( es,ev,a )
c       calc. eigen value and normalized eigen vector (3*3sym)
c     jancae_printinfo.
c       print inc.info. and element info.
c     jancae_printinout
c       print input and output of user subroutine
c
c
c-----------------------------------------------------------------------
c     zero clear vector a(n)
c
      subroutine jancae_clear1 ( a,n )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n)
c
      do i=1,n
        a(i)=0.0d0
      enddo
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     zero clear matrix a(n,m)
c
      subroutine jancae_clear2 ( a,n,m )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n,m)
c
      do i=1,n
        do j=1,m
          a(i,j)=0.0d0
        enddo
      enddo
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     zero clear matrix a(n,m,l)
c
      subroutine jancae_clear3 ( a,n,m,l )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n,m,l)
c
      do i=1,n
        do j=1,m
          do k=1,l
            a(i,j,k)=0.0d0
          enddo
        enddo
      enddo
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     set unit 2nd oder tensor [I]
c
      subroutine jancae_setunitm ( a,n )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n,n)
c
      call jancae_clear2 ( a,n,n )
      do i=1,n
        a(i,i)=1.0d0
      enddo
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     print vector a(n) with text
c
      subroutine jancae_print1 ( text,a,n )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n)
      character text*32
c
      write (6,*) text
      write (6,9000) (a(i),i=1,n)
 9000 format (6e16.8)
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     print matrix a(n,m) with text
c
      subroutine jancae_print2 ( text,a,n,m )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n,m)
      character text*32
c
      write (6,*) text
      do i=1,n
        write (6,9000) (a(i,j),j=1,m)
      enddo
 9000 format (6e16.8)
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     calculate multiplication of matrix and vector
c     {v}=[a]{u}
c
      subroutine jancae_mv (v,a,u,nv,nu)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension v(nv),a(nv,nu),u(nu)
c
      call jancae_clear1 ( v,nv )
      do i=1,nv
        do j=1,nu
          v(i)=v(i)+a(i,j)*u(j)
        enddo
      enddo
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     calculate multiplication of matrix and matrix
c     [a]=[b][c]
c
      subroutine jancae_mm (a,b,c,na1,na2,nbc)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(na1,na2),b(na1,nbc),c(nbc,na2)
c
      call jancae_clear2 ( a,na1,na2 )
      do i=1,na1
        do j=1,na2
          do k=1,nbc
            a(i,j)=a(i,j)+b(i,k)*c(k,j)
          enddo
        enddo
      enddo
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     calculate scaler product of vectors
c     s={u}T{v}
c
      subroutine jancae_vvs ( s,u,v,n )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension v(n),u(n)
c
      s=0.0
      do i=1,n
        s=s+u(i)*v(i)
      enddo
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     calculate inverse matrix using LU decomposition
c     [b]=[a]-1
c
      subroutine jancae_minv ( b,a,n,d )
c
c     Ref. http://astr-www.kj.yamagata-u.ac.jp/~shibata/kbg/
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n,n),b(n,n)
      dimension indx(n),y(n),c(n,n),aorg(n,n)
      character text*32
      logical   check
c
      check=.false.
      eps=1.0d-36
c
      do i=1,n
        do j=1,n
          aorg(i,j)=a(i,j)
        enddo
      enddo
c
      anorm=0.0
      do i=1,n
        do j=1,n
          if ( anorm.lt.abs(a(i,j)) ) anorm=abs(a(i,j))
        enddo
      enddo
      do i=1,n
        do j=1,n
          a(i,j)=a(i,j)/anorm
        enddo
      enddo
c
      if ( n.eq.2 ) then
        call jancae_minv2 ( b,a,d,eps )
        goto 100
      elseif ( n.eq.3 ) then
        call jancae_minv3 ( b,a,d,eps )
        goto 100
      endif
c
      call jancae_ludcmp ( a,n,indx,d,eps )
c                                                 ---- check determinant
      if ( abs(d).le.eps ) then
         write (6,*) 'determinant det[a] error',d
         write (6,*) 'stop in minv'
         call jancae_exit(9000)
      endif
c                                                            ---- B=A^-1
      do j=1,n
        call jancae_clear1 ( y,n )
        y(j)=1.0d0
        call jancae_lubksb ( a,n,indx,y,eps )
        do i=1,n
          b(i,j)=y(i)
        enddo
      enddo
c
  100 continue
      ani=1.0d0/anorm
      do i=1,n
        do j=1,n
          a(i,j)=aorg(i,j)
          b(i,j)=b(i,j)*ani
        enddo
      enddo
c                                                             ---- check
      if ( check ) then
        write (6,*) 'check inverse matrix',n
        text='original matrix [A]'
        call jancae_print2 ( text,a,n,n )
        text='inversed matrix [A]^-1'
        call jancae_print2 ( text,b,n,n )
        call jancae_clear2 ( c,n,n )
        do i=1,n
          do j=1,n
            do k=1,n
              c(i,j)=c(i,j)+b(i,k)*a(k,j)
            enddo
          enddo
        enddo
        text='[A]^-1*[A]=[I] ?'
        call jancae_print2 ( text,c,n,n )
      endif
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     LU decomposition
c
      subroutine jancae_ludcmp( a,n,indx,d,eps )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n,n),indx(n)
      dimension vtemp(n)
      character text*32
c
      d=1.0d0
      do i=1,n
        aamax=0.0d0
        do j=1,n
          if ( abs(a(i,j)).gt.aamax ) aamax=abs(a(i,j))
        enddo
        if ( aamax.le.eps ) then
          write (6,*) 'singular matrix in jancae_ludcmp'
          text='matrix detail'
          call jancae_print2 ( text,a,n,n )
          call jancae_exit ( 9000 )
        endif
        vtemp(i)=1.0d0/aamax
      enddo
c
      do j=1,n
        do i=1,j-1
          sum=a(i,j)
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
        enddo
        aamax=0.0d0
        do i=j,n
          sum=a(i,j)
          do k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
          dum=vtemp(i)*abs(sum)
          if ( dum.ge.aamax ) then
            imax=i
            aamax=dum
          endif
        enddo
        if ( j.ne.imax ) then
          do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo
          d=-d
          vtemp(imax)=vtemp(j)
        endif
        indx(j)=imax
c       if ( abs(a(i,j)).le.eps ) a(i,j)=eps     !2010.07.02 c.out
        if ( j.ne.n ) then
          ajj=a(j,j)                             !2010.07.02 add
          if ( abs(ajj).le.eps ) ajj=eps         !2010.07.02 add
          dum=1.0d0/ajj                          !2010.07.02 mod
          do i=j+1,n
            a(i,j)=a(i,j)*dum
          enddo
        endif
      enddo
c                                                 ---- get the det. of A
      do j=1,n
        d=d*a(j,j)
      enddo
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     LU backward substitution
c
      subroutine jancae_lubksb( a,n,indx,b,eps )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n,n),b(n),indx(n)
c
      ii=0
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if ( ii.ne.0 ) then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
          enddo
        else if ( abs(sum).ge.eps ) then
          ii=i
        endif
        b(i)=sum
      enddo
      do i=n,1,-1
        sum=b(i)
        do j=i+1,n
          sum=sum-a(i,j)*b(j)
        enddo
        b(i)=sum/a(i,i)
      enddo
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     calculate inverse matrix[2,2]
c
      subroutine jancae_minv2 ( b,a,deta,eps )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension b(2,2),a(2,2)
c
c
      deta=a(1,1)*a(2,2)-a(1,2)*a(2,1)
      if ( abs(deta).le.eps ) then
         write (6,*) 'determinant det[a] error',deta
         write (6,*) 'stop in minv2'
         call jancae_exit(9000)
      endif
c
      detai=1.0d0/deta
      b(1,1)= a(2,2)*detai
      b(1,2)=-a(1,2)*detai
      b(2,1)=-a(2,1)*detai
      b(2,2)= a(1,1)*detai
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     calculate inverse matrix[3,3]
c
      subroutine jancae_minv3 ( b,a,deta,eps )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension b(3,3),a(3,3)
c
      deta=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))+
     &     a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3))+
     &     a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
      if ( abs(deta).le.eps ) then
         write (6,*) 'determinant det[a] error',deta
         write (6,*) 'stop in minv3'
         call jancae_exit(9000)
      endif
c
      detai=1.0d0/deta
      b(1,1)=(a(2,2)*a(3,3)-a(2,3)*a(3,2))*detai
      b(1,2)=(a(1,3)*a(3,2)-a(1,2)*a(3,3))*detai
      b(1,3)=(a(1,2)*a(2,3)-a(1,3)*a(2,2))*detai
      b(2,1)=(a(2,3)*a(3,1)-a(2,1)*a(3,3))*detai
      b(2,2)=(a(1,1)*a(3,3)-a(1,3)*a(3,1))*detai
      b(2,3)=(a(1,3)*a(2,1)-a(1,1)*a(2,3))*detai
      b(3,1)=(a(2,1)*a(3,2)-a(2,2)*a(3,1))*detai
      b(3,2)=(a(1,2)*a(3,1)-a(1,1)*a(3,2))*detai
      b(3,3)=(a(1,1)*a(2,2)-a(1,2)*a(2,1))*detai
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     print informations for debug (intro)
c
      subroutine jancae_printinfo ( inc,nnrm,nshr )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /jancae1/ne,ip,lay
c
      nttl=nnrm+nshr
c
      write (6,*) '----- JANCAE.UMMDp debug info. -----'
      write (6,*) 'increment=',inc
      write (6,*) 'elem,ip,lay=',ne,ip,lay
      write (6,*) 'nttl,nnrm,nshr=',nttl,nnrm,nshr
      nerr=0
      if ( nnrm.eq.3 ) then
        if ( nshr.eq.3 ) then
          write (6,*) '3d solid element'
        else if ( nshr.eq.1 ) then
          write (6,*) 'plane strain or axi-sym solid element'
        else
          nerr=nerr+1
        endif
      else if ( nnrm.eq.2 ) then
        if ( nshr.eq.1 ) then
          write (6,*) 'plane stress or thin shell element'
        else if ( nshr.eq.3 ) then
          write (6,*) 'thick shell element'
        else
          nerr=nerr+1
        endif
      else
        nerr=nerr+1
      endif
      if ( nerr.ne.0 ) then
        write (6,*) 'no supported element type',nnrm,nshr
        call jancae_exit (9000)
      endif
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     print informations for debug (input/output)
c
      subroutine jancae_printinout ( io,s,de,d,nttl,
     &                               stv,nstv )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension s(nttl),de(nttl),d(nttl,nttl),stv(nstv)
      character text*32
c
      if ( io.eq.0 ) then
        text='initial stresses'
      else
        text='updated stresses'
      endif
      call jancae_print1 ( text,s,nttl )
c
      if ( io.eq.0 ) then
        text='initial internal state var.'
      else
        text='updated internal state var.'
      endif
      call jancae_print1 ( text,stv,nstv )
c
      if ( io.eq.0 ) then
        text='driving strain increment'
        call jancae_print1 ( text,de,nttl )
      else
        text='tangent modulus matrix'
        call jancae_print2 ( text,d,nttl,nttl )
      endif
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c    calculate eigen value and eigen vector by jacobi method
c
c    ( this subroutine is copied from http://www.flagshyp.com/ )
c
c    input  :  a(3,3)  : symmetric matrix to be analyzed
c    output :  es(i)   : i-th eigen value
c              ev(i,3) : normalized eigen vector for i-th eigen value
c
      subroutine jancae_eigen_sym3 ( es,ev,a )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension   a(3,3),es(3),ev(3,3)
      dimension   w(3,3),prc(3,3)
c
      msweep=100
      eps=1.0d-8
c
c                                                       ---- preparation
      ax=0.0d0
      er=0.0d0
      do i=1,3
        do j=1,3
          if ( abs(a(i,j)).gt.ax ) ax=abs(a(i,j))
          er=er+abs(a(i,j)-a(j,i))
        enddo
      enddo
      if ( er/ax.gt.eps ) then
        write (6,*) 'a is not symmetric'
        write (6,*) 'stop in jancae_eigen_sym3'
        call jancae_exit (9000)
      endif
      do i=1,3
        do j=1,3
          w(i,j)=a(i,j)/ax
        enddo
      enddo
c                                    ---- initialise prc to the identity
      do i=1,3
        do j=1,3
          prc(i,j)=0.0d0
        end do
        prc(i,i)=1.0d0
        es(i)=w(i,i)
      enddo
c
c                                                   ---- starts sweeping
c
      do is=1,msweep
c
        sum=0.0d0
        do ip=1,2
          do iq=ip+1,3
            sum=sum+abs( w(ip,iq) )
          enddo
        enddo
c       write (6,*) 'ite',is,sum,eps
c
c            ---- if the sum of off-diagonal terms is zero evaluates the
c                                                     esches and returns
c
        if ( abs(sum).lt.eps ) then
          do i=1,3
            do j=1,3
              ev(i,j)=prc(j,i)
            enddo
            es(i)=es(i)*ax
          enddo
          return
        endif
c
c                             ---- performs the sweep in three rotations
c                                         ---- one per off diagonal term
c
        do ip=1,2
          do iq=ip+1,3
            od=100.0d0*abs( w(ip,iq) )
            if ( abs(od).gt.eps ) then
              hd=es(iq)-es(ip)
c
c                                      ---- evaluates the rotation angle
c
              theta=0.5d0*hd/w(ip,iq)
              t=1.0d0/(abs(theta)+sqrt(1.0d0+theta**2))
              if ( theta.lt.0.0d0 ) t=-t
c
c                                   ---- re-evaluates the diagonal terms
c
              c=1.0d0/sqrt(1.0d0+t**2)
              s=t*c
              tau=s/(1.0d0+c)
              h=t*w(ip,iq)
              es(ip)=es(ip)-h
              es(iq)=es(iq)+h
c
c                     ---- re-evaluates the remaining off-diagonal terms
c
              ir=6-ip-iq
              g=w( min(ir,ip),max(ir,ip) )
              h=w( min(ir,iq),max(ir,iq) )
              w( min(ir,ip),max(ir,ip) )=g-s*(h+g*tau)
              w( min(ir,iq),max(ir,iq) )=h+s*(g-h*tau)
c
c                                          ---- rotates the eigenvectors
c
              do ir=1,3
                g=prc(ir,ip)
                h=prc(ir,iq)
                prc(ir,ip)=g-s*(h+g*tau)
                prc(ir,iq)=h+s*(g-h*tau)
              enddo
            endif
            w(ip,iq)=0.0d0
          enddo
        enddo
      enddo
c
c                              ---- if convergence is not achieved stops
c
      write (6,*) 'did not converge in eigen calculation.'
      write (6,*) 'msweep=',msweep
      write (6,*) 'eps=',eps
      write (6,*) 'sum=',sum
      write (6,*) 'stop in jancae_eigen_sym3'
      call jancae_exit(9000)
c
      end
c
c
c
c-----------------------------------------------------------------------
c     checking existence of file named "flname"
c
       logical function jancae_file_exist ( flname )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character   flname*16
c
      nio=616
      open  ( nio,file=flname,status='old',err=10 )
c
      close ( nio,            status='keep' )
      jancae_file_exist=.true.
      return
c
   10 jancae_file_exist=.false.
      return
c
      end
c
c
c
c--------------------------------------------------------------(bbc2005)
c     BBC2005 Yield Function
c
      subroutine jancae_bbc2005 ( s,se,dseds,d2seds2,nreq,
     &                            pryld,ndyld )
c-----------------------------------------------------------------------
      implicit none

c                                                    ---- input & output
      integer, intent(in) :: nreq , ndyld
      real*8, intent(in) :: s(3), pryld(ndyld)
      real*8, intent(out) :: se, dseds(3), d2seds2(3,3)
c                                                   ---- local variables
      real*8 a, b, L, M, N, P, Q, R, nn, fir, las, oo,
     &  th(3),lth(3),pp,
     &  phi, Al, AA, BB, CC, DD, kk, lth1_2, lth12,
     &  lth2_3, lth23, d2phi_a, d2phi_b,
     &  dlthds(3,3), dphidlth(3), dsedphi, d2sedphi2,
     &  d2pdg2, k_1, k_2,
     &  d2pdp2, k_a, k_b, k_c,
     &  d2phidlth2(3,3),d2gds2(3,3),d2lds2(3,3),d2pds2(3,3),
     &  d2lthds2(3,3,3)
      integer i, j, d, e, k, mm, fact, ii, jj
c
c     local variables : symbols in BBC2005 Yield Functions document
c
c     In this subroutine, intermidiate function th(1), th(2) and th(3)
c     of the original paper are represented with Theta(i).
c     th(1) =L*s(1)*M*s(2)
c     th(2) =sqrt((N*s(1)-P*s(2))**2+s(3)**2)
c     th(3) =sqrt((Q*s(1)-R*s(2))**2+s(3)**2)
c     k,a,b,L,M,N,P,Q,R    :material parameter
c     phi                  :phi
c     Al                   :A
c     lth(1) = th(1)**2
c     lth(2) = th(2)**2
c     lth(3) = th(3)**2
c     dlthds(3,3)            :d(lth(i))/d(sigma_x,x,txy)
c     dphidlth(3)            :d(phi)/d(lth(i))
c     dsephi                 :d(se)/d(phi)
c     d2phidlth2(3,3)        :d2(phi)/dg(lth(i))2 ,i =1,2,3
c     d2lthds2(i,3,3)        :d2(lth(i))/d(sigma_x,y,txy)2, i=1,2,3
c
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
      do i = 0, mm
        nn = nn + fact(k)/(fact(k-2*i)*fact(2*i))
     &       *AA**(k-2*i)*BB**(2*i)
        oo = oo + fact(k)/(fact(k-2*i)*fact(2*i))
     &       *CC**(k-2*i)*DD**(2*i)
      enddo
c
      fir = 2.0d0*a*nn
      las = 2.0d0*b*oo
c
      phi = fir + las
c
      kk = dble(k)
c
      Al=(a*(N+L)**(2*k)+a*(N-L)**(2*k)
     &   +b*(N+Q)**(2*k)+b*(N-Q)**(2*k))**(1.0d0/(2.0d0*kk))
c
      se = phi**(1/(2*kk))/Al
c
      dseds(:) = 0
      d2seds2(:,:) = 0
c
c
c                                           ----  1st order differential
      if ( nreq.ge.1 ) then
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
      endif
c
      if ( lth(2) < 1e-15 * se**2) then
        lth(2) = 1e-15 * se**2
      endif
c
      if ( lth(3) < 1e-15 * se**2) then
        lth(3) = 1e-15 * se**2
      endif
c
      lth1_2 = lth(1)+lth(2)
      lth12  = lth(1)*lth(2)
      lth2_3 = lth(2)+lth(3)
      lth23  = lth(2)*lth(3)
c
      if (lth1_2 < 1e-15 * se**2) then
        lth1_2 = 1e-15 * se**2
      endif
c
      if (lth2_3 < 1e-15 * se**2) then
        lth2_3 = 1e-15 * se**2
      endif
c
      dphidlth(:) = 0.0d0
c
      do i = 0, mm
        dphidlth(1) = dphidlth(1) + fact(k)/(fact(k-2*i)*fact(2*i))*
     &     (2*a*(i*4**i*lth(2)**i*lth(1)**(i-1)*lth1_2
     &     **(k-2*i)+(k-2*i)*4**i*lth12**i*lth1_2**(-2*i+k-1)))
c
        dphidlth(2) = dphidlth(2) + fact(k)/(fact(k-2*i)*fact(2*i))*
     &    (2*a*(i*4**i*lth(2)**(i-1)*lth(1)**i*lth1_2
     &    **(k-2*i)+(k-2*i)*4**i*lth12**i*
     &    lth1_2**(-2*i+k-1))
     &    +2*b*(i*4**i*lth(2)**(i-1)*lth(3)**i*lth2_3
     &    **(k-2*i)+(k-2*i)*4**i*lth23**i*
     &    lth2_3**(-2*i+k-1)))
c
        dphidlth(3) = dphidlth(3) + fact(k)/(fact(k-2*i)*fact(2*i))*
     &     (2*b*(i*4**i*lth(2)**i*lth(3)**(i-1)*lth2_3
     &     **(k-2*i)+(k-2*i)*4**i*lth23**i*lth2_3**(-2*i+k-1)))
      enddo
c
      dlthds(1,1) = 2*L*(M*s(2)+L*s(1))
      dlthds(1,2) = 2*M*(M*s(2)+L*s(1))
      dlthds(1,3) = 0d0
      dlthds(2,1) = 2*N*(N*s(1)-P*s(2))
      dlthds(2,2) = -2*P*(N*s(1)-P*s(2))
      dlthds(2,3) = 2*s(3)
      dlthds(3,1) = 2*Q*(Q*s(1)-R*s(2))
      dlthds(3,2) = -2*R*(Q*s(1)-R*s(2))
      dlthds(3,3) = 2*s(3)
c
c
      do i = 1, 3
        dseds(i) = 0.0d0
        do j = 1, 3
          dseds(i) = dseds(i) + dsedphi*dphidlth(j)*dlthds(j,i)
        enddo
c        write (150,*) s(1), s(2), s(3), i, dseds(i)
      enddo
c
      endif
c
c
c                                            ---- 2nd order differential
c
      if ( nreq.ge.2 ) then
c
      d2sedphi2 = (1/kk/2.0d00-1)*phi**(1/kk/2.0d0-2)/
     &            (kk*Al)/2.0d0
c
      d2phidlth2(:,:) = 0.0d0
c
      do i = 0, mm
c
      k_a = (i-1)*i*4**i
      k_b = (k-2*i)*i*4**i
      k_c = (k-2*i-1)*(k-2*i)*4**i
c
      d2phidlth2(1,1) = d2phidlth2(1,1)+fact(k)/(fact(k-2*i)*fact(2*i))
     &     *2*a*(k_a*lth(1)**(i-2)*lth(2)**i*lth1_2**(k-2*i)
     &     +2*(k_b*lth(2)**i*lth(1)**(i-1)*lth1_2**(k-2*i-1))
     &     +k_c*lth12**i*lth1_2**(k-2*i-2))
c
      d2phidlth2(1,2) = d2phidlth2(1,2)+fact(k)/(fact(k-2*i)*fact(2*i))
     &     *2*a*(i**2*4**i*lth(1)**(i-1)*lth(2)**(i-1)*lth1_2**(k-2*i)
     &     +k_b*lth(1)**i*lth(2)**(i-1)*lth1_2**(k-2*i-1)
     &     +k_b*lth(1)**(i-1)*lth(2)**i*lth1_2**(k-2*i-1)
     &     +(k-2*i-1)*(k-2*i)*4**i*lth12**i*lth1_2**(k-2*i-2))
c
      d2phidlth2(2,2) = d2phidlth2(2,2)+fact(k)/(fact(k-2*i)*fact(2*i))
     &     *(2*b*(k_a*lth(2)**(i-2)*lth(3)**i*lth2_3**(k-2*i)
     &     +2*k_b*lth(3)**i*lth(2)**(i-1)*lth2_3**(k-2*i-1)
     &     +k_c*lth23**i*lth2_3**(k-2*i-2))
     &     +2*a*(k_a*lth(2)**(i-2)*lth(1)**i*lth1_2**(k-2*i)
     &     +2*k_b*lth(1)**i*lth(2)**(i-1)*lth1_2**(k-2*i-1)
     &     +k_c*lth12**i*lth1_2**(k-2*i-2)))
c
      d2phidlth2(2,3) = d2phidlth2(2,3)+fact(k)/(fact(k-2*i)*fact(2*i))
     &     *2*b*(i**2*4**i*lth(2)**(i-1)*lth(3)**(i-1)*lth2_3**(k-2*i)
     &     +k_b*lth(2)**(i-1)*lth(3)**i*lth2_3**(k-2*i-1)
     &     +k_b*lth(2)**i*lth(3)**(i-1)*lth2_3**(k-2*i-1)
     &     +k_c*lth23**i*lth2_3**(k-2*i-2))
c
      d2phidlth2(3,3) = d2phidlth2(3,3)+fact(k)/(fact(k-2*i)*fact(2*i))
     &     *2*b*(k_a*lth(2)**i*lth(3)**(i-2)*lth2_3**(k-2*i)
     &     +2*k_b*lth(2)**i*lth(3)**(i-1)*lth2_3**(k-2*i-1)
     &     +k_c*lth23**i*lth2_3**(k-2*i-2))
c
      enddo
c
      d2phidlth2(2,1) = d2phidlth2(1,2)
      d2phidlth2(3,2) = d2phidlth2(2,3)
c
c
      d2lthds2(:,:,:) = 0.0d0
c
      d2lthds2(1,1,1) = 2*L**2
      d2lthds2(1,1,2) = 2*L*M
      d2lthds2(1,2,1) = d2lthds2(1,1,2)
      d2lthds2(1,2,2) = 2*M**2
c
      d2lthds2(2,1,1) = 2*N**2
      d2lthds2(2,1,2) = -2*N*P
      d2lthds2(2,2,1) = d2lthds2(2,1,2)
      d2lthds2(2,2,2) = 2*P**2
      d2lthds2(2,3,3) = 2.0d0
c
      d2lthds2(3,1,1) = 2*Q**2
      d2lthds2(3,1,2) = -2*Q*R
      d2lthds2(3,2,1) = d2lthds2(3,1,2)
      d2lthds2(3,2,2) = 2*R**2
      d2lthds2(3,3,3) = 2.0d0
c
      do i=1,3
        do j=1,3
          d2seds2(i,j)=0.0d0
          do ii=1,3
            do jj=1,3
              d2seds2(i,j) = d2seds2(i,j)
     &+d2sedphi2*dphidlth(jj)*dlthds(jj,j) *dphidlth(ii)*dlthds(ii,i)
     &+dsedphi*     d2phidlth2(ii,jj)*dlthds(jj,j)     *dlthds(ii,i)
            enddo
            d2seds2(i,j)= d2seds2(i,j)
     &+dsedphi*dphidlth(ii)                   *d2lthds2(ii,i,j)
          enddo
        enddo
      enddo
c
      endif
c
      end subroutine jancae_bbc2005
c
c
c
c----  function to solve factorial ----------------------------(bbc2005)
      integer function fact(n) result(m)
      integer, intent(in) :: n
      m=1
      if ( n.ge.1 ) then
          do i = 1,n
            m = m * i
          enddo
      endif
      end function fact
c
c
cc--------------------------------------------------------------(bbc2008)
c
c     BBC2008 yield function & its differentials for JANCAE UMMDp
c
c     Title: Plane-Stress Yield Criterion For Highly-Anisotropic
c            Sheet Metals
c     Authors: Dan-Sorin Comsa, Dorel Banabic
c              (Technical University of Cluj-Napoca, 15C. Daicoviciu,
c              400020 Cluj-Napoca, Romania)
c     Paper: Numisheet 2008, Sep. 1-5, 2008 - Interlaken, Switzerland
c
c
c-----------------------------------------------------------------------
c     coded by Kunio TAKEKOSHI (Terrabyte Co., Ltd)
c      - January  7, 2011 Initiall release
c      - January 15, 2011 Rewrite explanations.
c      - February 24, 2015 clean up (deleted some comments etc.)
c
c     modified by H.Takizawa (JANCAE)
c      - February 23, 2015
c--------------------------------------------------------------(bbc2008)
      subroutine jancae_bbc2008 (s,se,dseds,d2seds2,nreq,
     &                           pryld,ndyld)
c
      implicit none
c
c                                                    ---- input & output
      integer, intent(in) :: nreq , ndyld
      real*8, intent(in) :: s(3), pryld(ndyld)
      real*8, intent(out) :: se, dseds(3), d2seds2(3,3)

c                                                   ---- local variables
c
c     nds and ndk are used as s and k in the literature of BBC2008,
c       respectively.
c
      integer :: nds, ndk
c
      nds=nint(pryld(2))
      ndk=nint(pryld(3))
c
      call jancae_bbc2008_core (s,se,dseds,d2seds2,nreq,
     &                          pryld,ndyld,nds,ndk)
c
      return
      end subroutine jancae_bbc2008
c
c
c     <<local variables>>
c
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
c--------------------------------------------------------------(bbc2008)
      subroutine jancae_bbc2008_core ( s , se , dseds , d2seds2 , nreq ,
     &                                 pryld , ndyld , sp , kp )
c
      implicit none
c
c                                                    ---- input & output
      integer, intent(in) :: nreq , ndyld , sp , kp
      real*8, intent(in) :: s(3), pryld(ndyld)
      real*8, intent(out) :: se, dseds(3), d2seds2(3,3)
c
c                                                   ---- local variables
      integer csp, m, eta
c
c                                       ---- local variables for BBC2008
      real*8  :: se1, wp
      real*8  :: Lp(sp,3,3), Mp(sp,3,3), Np(sp,3,3)
      real*8  :: kCm(0:kp)
c
c                                               ---- temporary variables
      real*8 wpi(2)
      real*8 dFds(3), d2Fds2(3,3)
      real*8 phiL, phiM, phiN
      real*8 dphiLds(3), dphiMds(3), dphiNds(3)
      real*8 d2phiLds2(3,3), d2phiMds2(3,3), d2phiNds2(3,3)
c                                    ---- phiL_m = phiL^m, se2k = se^2kp
      real*8 phiL_m, phiM_kp_m, phiN_m, se2k
c
c                                                   ---- local functions
      real*8 jancae_bbc2008_get_se
c
c     ----------------
c        parameters
c     ----------------
c
      call jancae_bbc2008_setup (pryld, ndyld, sp, kp, wp,
     &                           Lp, Mp, Np, kCm, se1)
c
c ----------------
c  se section
c ----------------
c
c                             ---- The unit of this se is (stress)^(2kp)
      se = jancae_bbc2008_get_se (sp, kp, wp, s, Lp, Mp, Np, kCm)
c
c      ---- see eq.(x.y.2b) and eq.(x.y.2) for se2k and se, respectively
      se2k = se / se1
      se = se2k**(1.0d0 / (2.0d0*kp))
c
c   ---- If a main routine requests this subroutine to calculate se only
      if (nreq.eq.0) then
        return
      endif
c
c
c --------------------------
c  dseds & d2seds2 section
c --------------------------
      call jancae_clear1(dFds,3)
      call jancae_clear2(d2Fds2,3,3)

c                              ---- long-long-long do loops starts here.
      do csp=1,sp
c
        call jancae_bbc2008_get_w_phi (wpi,phiL,phiM,phiN,
     &                                 csp,sp,wp,Lp,Mp,Np,s)
c
        do m=0,kp
c
c     phiM^m, phiL^(kp-m) and phiN^m terms sometimes become 0**0.
c     To get consistency with the yield function and its differentials
c       disscussed by Banabic et al., 0**0 is required to be 1.
c
          phiL_m = 1.0d0
          phiN_m = 1.0d0
          if (m.ne.0) then
            phiL_m = phiL**m
            phiN_m = phiN**m
          endif
c
          phiM_kp_m = 1.0d0
          if ((kp-m).ne.0) then
            phiM_kp_m = phiM**(kp-m)
          endif
c
c
          call jancae_bbc2008_get_dphiXds (dphiLds,Lp,s,csp,m,sp)
          call jancae_bbc2008_get_dphiXds (dphiMds,Mp,s,csp,(kp-m),sp)
          call jancae_bbc2008_get_dphiXds (dphiNds,Np,s,csp,m,sp)
c
c                                             ---- <dseds>, see (x.y.2f)
          dFds(1:3) = dFds(1:3) + kCm(m) *
     &                    (  wpi(1)*(  dphiMds(1:3) * phiL_m
     &                               + dphiLds(1:3) * phiM_kp_m )
     &                     + wpi(2)*(  dphiMds(1:3) * phiN_m
     &                               + dphiNds(1:3) * phiM_kp_m ))
c
c
c                     ---- <d2seds2>, see (x.y.2g), d2F/ds(eta)ds(gamma)
          if (nreq.eq.2) then
c
            call jancae_bbc2008_get_d2phiXds2 (d2phiLds2,Lp,s,csp,m,sp)
            call jancae_bbc2008_get_d2phiXds2 (d2phiMds2,Mp,s,csp,
     &                                                       (kp-m),sp)
            call jancae_bbc2008_get_d2phiXds2 (d2phiNds2,Np,s,csp,m,sp)
c
            do eta=1,3
              d2Fds2(eta,1:3) = d2Fds2(eta,1:3) + kCm(m)
     &             * (  wpi(1) * (  d2phiMds2(eta,1:3) * phiL_m
     &                            + dphiMds(1:3) * dphiLds(eta)
     &                            + dphiLds(1:3) * dphiMds(eta)
     &                            + d2phiLds2(eta,1:3) * phiM_kp_m )
     &                + wpi(2) * (  d2phiMds2(eta,1:3) * phiN_m
     &                            + dphiMds(1:3) * dphiNds(eta)
     &                            + dphiNds(1:3) * dphiMds(eta)
     &                            + d2phiNds2(eta,1:3) * phiM_kp_m ))
            enddo
c
          endif
c
c                                                ---- end of m=0,kp loop
        enddo
c
c                                                ---- end of i=1,sp loop
      enddo
c
c
c                            ---- < dseds >, see (x.y.2f), se2k = se^2kp
      dseds(1:3) =  dFds(1:3) * se / (se1 * 2.0d0 * kp * se2k)

c                  ---- < d2seds2 >, see (x.y.2g), d2se/ds(eta)ds(gamma)
      if (nreq.eq.2) then
        do eta=1,3
          d2seds2(eta,1:3) =
     &            d2Fds2(eta,1:3) * se / (se1 * 2.0d0 * kp * se2k)
     &            - (2.0d0*kp-1.0d0) * dseds(eta) * dseds(1:3) / se
        enddo
      endif
c
      return
      end subroutine jancae_bbc2008_core
c
c
c
c--------------------------------------------------------------(bbc2008)
c     jancae_bbc2008_get_w_phi ()
c     A subroutine to get w^(i-1), w^(s-i) and phiX variables
c-----------------------------------------------------------------------
      subroutine jancae_bbc2008_get_w_phi
     &                  (wpi, phiL, phiM, phiN,
     &                   csp, sp, wp, Lp, Mp, Np, s)
c
      implicit none
c
      real*8, intent(out) :: wpi(2), phiM, phiL, phiN
      integer, intent(in) :: csp, sp
      real*8, intent(in) :: wp
      real*8, intent(in) :: s(3), Lp(sp,3,3)
      real*8, intent(in) :: Mp(sp,3,3), Np(sp,3,3)
c
      real*8 jancae_bbc2008_get_phiX
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
c--------------------------------------------------------------(bbc2008)
c     jancae_bbc2008_get_se (sp, kp, wp, s, Lp, Mp, Np, kCm)
c     get equiv. stress
c
c     Caution!
c      The unit of get_se is (stress)^(2*kp).
c-----------------------------------------------------------------------
      real*8 function jancae_bbc2008_get_se (sp,kp,wp,s,Lp,Mp,Np,kCm)
c
      implicit none
c
      integer, intent(in) :: sp, kp
      real*8, intent(in) :: wp, s(3)
      real*8, intent(in) :: Lp(sp,3,3)
      real*8, intent(in) :: Mp(sp,3,3), Np(sp,3,3)
      real*8, intent(in) :: kCm(0:kp)
c
      integer csp, m
      real*8 phiM, phiL, phiN
      real*8 wpi(2), phiL_m, phiM_kp_m, phiN_m
c
      jancae_bbc2008_get_se = 0.0d0
c
      do csp=1,sp
c
        call jancae_bbc2008_get_w_phi
     &                 (wpi, phiL, phiM, phiN,
     &                  csp, sp, wp, Lp, Mp, Np, s)
c
        do m=0,kp
c
          phiL_m = 1.0d0
          phiN_m = 1.0d0
          if (m.ne.0) then
            phiL_m = phiL**m
            phiN_m = phiN**m
          endif
c
          phiM_kp_m = 1.0d0
          if ((kp-m).ne.0) then
            phiM_kp_m = phiM**(kp-m)
          endif
c
          jancae_bbc2008_get_se =
     &    jancae_bbc2008_get_se
     &     + kCm(m) * phiM_kp_m * ( wpi(1) * phiL_m + wpi(2) * phiN_m )
c
        enddo
c
      enddo
c
      return
      end
c
c
c
c--------------------------------------------------------------(bbc2008)
c     jancae_bbc2008_get_phiX (Xp, s, csp , sp)
c     A function to calculate s(a)*X(a,b)*s(b) (summation convention)
c
c     Xp: coefficient tensor
c     s : stress vector
c     csp: CURRENT sp value, 1<=csp<=sp
c
c     local variables
c     nc: the number of components.
c     XXp: = Xp(csp, nc, nc)
c-----------------------------------------------------------------------
c
      real*8 function jancae_bbc2008_get_phiX (Xp, s, csp , sp)
c
      implicit none
      integer, parameter :: nc = 3
c
      integer, intent(in) :: csp , sp
      real*8, intent(in) :: Xp(sp, nc, nc), s(nc)
c
      integer i
      real*8 XXp(nc, nc), v(nc)
c
c                                  ---- convert 3rd tensor to 2nd tensor
      do i=1,nc
        XXp(i,1:nc) = Xp(csp, i, 1:nc)
      enddo
c
      call jancae_mv (v, XXp, s, nc, nc)
      call jancae_vvs (jancae_bbc2008_get_phiX, v, s, nc)
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     jancae_bbc2008_get_dphiXds (dphiXds, Xp, s, csp, lambda,sp)
c     A subroutine to calculate d(phiX^(lambda))/ds.
c     It returns dphiXds(nc).
c
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
c--------------------------------------------------------------(bbc2008)
c
      subroutine jancae_bbc2008_get_dphiXds (dphiXds,Xp,s,csp,lambda,sp)
c
      implicit none
      integer, parameter :: nc = 3
c
      integer, intent(in) :: csp, lambda , sp
      real*8, intent(in) :: Xp(sp, nc, nc), s(nc)
      real*8, intent(out) :: dphiXds(nc)
c
      integer i
      real*8 XXp(nc, nc), v(nc), phi
c
      call jancae_clear1(dphiXds,nc)
c
c                              ---- If lambda is 0, return dphiXds = {0}
      if (lambda.eq.0) then
        return
      endif
c
c                                  ---- convert 3rd tensor to 2nd tensor
      do i=1,nc
        XXp(i,1:nc) = Xp(csp, i, 1:nc)
      enddo
c
c     In the bbc2008 section of the document "User subroutines for
c       Metalic Plasticity model?", expression (x.y.2d) has
c       "chi(gamma, beta)*s(beta)" (summation convention).
c     In this routine, corresponding term is v(nc) obtained
c       from jancae_mv().
c
      call jancae_mv (v, XXp, s, nc, nc)
c
      if (lambda.eq.1) then
        dphiXds(1:nc) = 2.0d0 * v(1:nc)
      else
        call jancae_vvs (phi, v, s, nc)
        dphiXds(1:nc) = 2.0d0 * lambda * phi**(lambda-1) * v(1:nc)
      endif
c
      return
      end subroutine jancae_bbc2008_get_dphiXds
c
c
c
c-----------------------------------------------------------------------
c     jancae_bbc2008_get_d2phiXds2 (d2phiXds2, Xp, s, csp, lambda,sp)
c     A subroutine to calculate d2(phiX^(lambda))/(dsds').
c     It returns d2phiXdsds(nc,nc).
c
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
c--------------------------------------------------------------(bbc2008)
c
      subroutine jancae_bbc2008_get_d2phiXds2
     &                     (d2phiXds2,Xp,s,csp,lambda,sp)
c
      implicit none
      integer, parameter :: nc = 3
c
      integer, intent(in) :: csp, lambda ,sp
      real*8, intent(in) :: Xp(sp, nc, nc), s(nc)
      real*8, intent(out) :: d2phiXds2(nc, nc)
c
      integer i
      real*8 XXp(nc, nc), v(nc), phi, phi_lambda2
c
c                             ---- see eq.(x.y.2e), the case lambda <= 1
      if (lambda.le.1) then
        do i=1,nc
          d2phiXds2(i,1:nc) = 2.0d0 * lambda * Xp(csp, i, 1:nc)
        enddo
        return
      endif
c
c
      do i=1,nc
        XXp(i,1:nc) = Xp(csp, i, 1:nc)
      enddo
c
      call jancae_mv (v, XXp, s, nc, nc)
      call jancae_vvs (phi, v, s, nc)
c
      phi_lambda2 = 1.0d0
      if (lambda.ne.2) then
        phi_lambda2 = phi**(lambda-2)
      endif
c
      call jancae_clear2(d2phiXds2, nc, nc)
c
c                                            ---- d2phiX/(ds(i)ds(1:nc))
      do i=1,nc
        d2phiXds2(i, 1:nc) = 2.0d0 * lambda * phi_lambda2 *
     &   ( 2.0d0 * (lambda - 1) * v(1:nc) * v(i)
     &    + phi * XXp(1:nc, i) )
      enddo
c
      end subroutine jancae_bbc2008_get_d2phiXds2
c
c
c
c--------------------------------------------------------------(bbc2008)
c     setup_bbc2008_parameters()
c     A routine to setup local variables.
c-----------------------------------------------------------------------
      subroutine jancae_bbc2008_setup (pryld, ndyld, sp, kp, wp,
     &                                 Lp, Mp, Np, kCm, se1)
c
      implicit none
c
      integer, intent(in) :: sp, kp , ndyld
      real*8, intent(in) :: pryld(ndyld)
      real*8, intent(inout) :: wp, se1
      real*8, intent(inout) :: Lp(sp,3,3)
      real*8, intent(inout) :: Mp(sp,3,3), Np(sp,3,3)
      real*8, intent(inout) :: kCm(0:kp)
c
c                                                   ---- local variables
      integer csp, k, l, m, n
      real*8  Comb(kp*2, 0:kp*2)
c
c                                  ---- dummy variable to calculate se1.
      real*8 dummy_s(3)
c
c                                                          ---- function
      real*8 jancae_bbc2008_get_se
c
      wp = 1.5d0 ** (1.0d0/sp)
c
c                                        ---- Combination variables, kCm
      kCm(0) = 1.0d0
      kCm(kp) = 1.0d0
c
      do k=1,2*kp
c                                    ---- caution: Comb has zero origin.
        Comb(k,0) = 1.0d0
        Comb(k,k) = 1.0d0
        do m=1,k-1
          Comb(k, m) = Comb(k-1, m-1) + Comb(k-1, m)
        enddo
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
        if (k.eq.(2*kp)) then
          n = 1
          do m=2,k-2,2
            kCm(n) = Comb(k, m)
            n = n + 1
          enddo
        endif

      enddo
c
c
c                                                           ---- tensors
      do csp=1,sp
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
      enddo
c
c
c     equiv. stress in uniaxial stress state.
c     dummy_s = (1.0d0, 0.0d0, 0.0d0)
c     ** The unit of this se1 is (stress)^(2kp)
c
      call jancae_clear1 (dummy_s, 3)
      dummy_s(1) = 1.0d0
      se1 = jancae_bbc2008_get_se (sp, kp, wp, dummy_s, Lp, Mp, Np, kCm)

      return
      end subroutine jancae_bbc2008_setup
c
c
c
c--------------------------------------------------------------(bbc2008)
c     setup_MN_tensors (ic,csp,pryld,ndyld, Xp,sp)
c       ic: initial counter
c
c       This routine returns Mp or Np tensor.
c       Mp and Np tensors are the same style,
c       thus this subroutine has been created.
c-----------------------------------------------------------------------
      subroutine jancae_bbc2008_setup_MN_tensors (ic,csp,
     &                                      pryld,ndyld,Xp,sp)
c
      implicit none
c
      integer, intent(in) :: ic,csp,sp,ndyld
      real*8, intent(in) :: pryld(ndyld)
      real*8, intent(inout) :: Xp(sp,3,3)
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
cc-----------------------------------------------------------(cazacu2006)
c     Cazacu 2006 yield function and its dfferentials
c     ( IJP v.22(2006) p1171-1194. )
c
c     !!! CAUTION !!!
c     Plane stress condition is NOT implemented in this code.
c
      subroutine jancae_cazacu2006 ( s,se,dseds,d2seds2,nreq,
     &                               pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension s(6),dseds(6),d2seds2(6,6),pryld(ndyld)
c
      dimension sigma(6), psigma(3)
      dimension phi(3), psi(3), omega(3)
      dimension c(6,6), ct(6,6)
c
      dimension DFDH(3), DFDpsigma(3)
      dimension DpsigmaDH(3,3), DHdsigma(3,6), DsigmaDs(6,6)
      dimension DFDs(6)
c
      dimension D2FDpsigma2(3,3), D2psigmaDH2(3,3,3)
      dimension D2FDH2(3,3)
      dimension D2HDsigma2(3,6,6), D2FDs2(6,6)
c
      dimension s0(6)
      dimension dummat(3,6)
c
      parameter( pi = 3.141592653589793d0 )
      parameter( eps= 1.0d-5 )
c
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
c
c
c                                        ---- set anisotropic parameters
c
      call jancae_clear2 ( c,6,6 )
c
      c(1,1)= pryld(1+ 1)           ! C11
      c(1,2)= pryld(1+ 2)           ! C12
      c(1,3)= pryld(1+ 3)           ! C13
      c(2,1)= pryld(1+ 4)           ! C21
      c(2,2)= pryld(1+ 5)           ! C22
      c(2,3)= pryld(1+ 6)           ! C23
      c(3,1)= pryld(1+ 7)           ! C31
      c(3,2)= pryld(1+ 8)           ! C32
      c(3,3)= pryld(1+ 9)           ! C33
      c(4,4)= pryld(1+10)           ! C44     ! tau_xy
      c(5,5)= pryld(1+11)           ! C55     ! tau_xz
      c(6,6)= pryld(1+12)           ! C66     ! tau_yz
      a     = pryld(1+13)           ! a
      ck    = pryld(1+14)           ! k
c
      ai= 1.0d0/a
c
c                                         ---- Calculate phi, psi, omega
c
      phi(1)=(2.0d0*c(1,1)-c(1,2)-c(1,3))/3.0d0
      phi(2)=(2.0d0*c(1,2)-c(2,2)-c(2,3))/3.0d0
      phi(3)=(2.0d0*c(1,3)-c(2,3)-c(3,3))/3.0d0
c
      psi(1)=(-c(1,1)+2.0d0*c(1,2)-c(1,3))/3.0d0
      psi(2)=(-c(1,2)+2.0d0*c(2,2)-c(2,3))/3.0d0
      psi(3)=(-c(1,3)+2.0d0*c(2,3)-c(3,3))/3.0d0
c
      omega(1)=(c(1,1)+c(1,2)-2.0d0*c(1,3))/3.0d0
      omega(2)=(c(1,2)+c(2,2)-2.0d0*c(2,3))/3.0d0
      omega(3)=(c(1,3)+c(2,3)-2.0d0*c(3,3))/3.0d0
c
c      ---- Calculate 4th order orthotropic tensor "L" ( named ct here )
c
      call jancae_clear2( ct,6,6 )
c
      ct(1,1)= phi(1)
      ct(1,2)= psi(1)
      ct(1,3)=-omega(1)
      ct(2,1)= phi(2)
      ct(2,2)= psi(2)
      ct(2,3)=-omega(2)
      ct(3,1)= phi(3)
      ct(3,2)= psi(3)
      ct(3,3)=-omega(3)
      ct(4,4)= c(4,4)
      ct(5,5)= c(5,5)
      ct(6,6)= c(6,6)
c
c                               ---- Calculate linear transformed stress
c
      call jancae_mv( sigma,ct,s,6,6 )
c
c            ---- Calculate principal values of transformed stress sigma
c                                                     by Cardan's method
c
c                                       ---- 1st, 2nd and 3rd invariants
      H1=(sigma(1)+sigma(2)+sigma(3))/3.0d0
      H2=(sigma(5)**2.0d0+sigma(6)**2.0d0+sigma(4)**2.0d0
     1    -sigma(2)*sigma(3)-sigma(3)*sigma(1)-sigma(1)*sigma(2))/3.0d0
      H3=(2.0d0*sigma(5)*sigma(6)*sigma(4)+sigma(1)*sigma(2)*sigma(3)
     1    -sigma(1)*sigma(6)**2.0d0-sigma(2)*sigma(5)**2.0d0
     2    -sigma(3)*sigma(4)**2.0d0)/2.0d0 ! sigma(5) <-> sigma(6)
c
      p=H1**2.0d0+H2
      q=(2.0d0*H1**3.0d0+3.0d0*H1*H2+2.0d0*H3)/2.0d0
      if ( abs(p).ge.1.0d-16 ) then
        theta=q/(p**1.5d0)
        if ( theta >  1.0d0 ) theta= 1.0d0
        if ( theta < -1.0d0 ) theta=-1.0d0
        theta=acos(theta)
      else
        theta=0.0
      endif
c
c                               ---- Calculate principal values of sigma
      psigma(1)=2.0d0*sqrt(p)*cos(theta/3.0d0)+H1
      psigma(2)=2.0d0*sqrt(p)*cos((theta+4.0d0*pi)/3.0d0)+H1
      psigma(3)=2.0d0*sqrt(p)*cos((theta+2.0d0*pi)/3.0d0)+H1
c
c                                                 ---- Equivalent stress
c
c                                          ---- Calculate yield function
      F=(abs(psigma(1))-ck*psigma(1))**a
     1 +(abs(psigma(2))-ck*psigma(2))**a
     2 +(abs(psigma(3))-ck*psigma(3))**a
c                                           ---- Denominator coefficient
      D=(abs(phi(1))-ck*phi(1))**a
     1 +(abs(phi(2))-ck*phi(2))**a
     2 +(abs(phi(3))-ck*phi(3))**a
c
      se=(F/D)**ai
c
c                                            ---- 1st order differential
c
      if ( nreq >= 1 ) then
c                                              ---- D(se)/D(F) -> Scalar
c
        DseDF=(1.0d0/D)**ai*ai*F**(ai-1.0d0)
c
c                                    ---- D(F)/D(psigma) -> 1 x 3 Vector
c
        do i = 1,3
          DFDpsigma(i)=a*(psigma(i)/abs(psigma(i))-ck)
     1                 *(abs(psigma(i))-ck*psigma(i))**(a-1.0d0)
        enddo
c
c                                 ---- D(F)/D(H) by using D(psigma)/D(H)
c                                         D(F)/D(H)      -> 1 x 3 Vector
c                                         D(psigma)/D(H) -> 3 x 3 Matrix
c
        call jancae_clear1( DFDH,3 )
        call jancae_clear2( DpsigmaDH,3,3 )
c
        if ( abs(psigma(2)-psigma(3))/se > eps .and.
     1       abs(psigma(2)-psigma(1))/se > eps ) then
c                                                 ---- Not Singular case
          do i=1,3
            denom=psigma(i)**2.0d0-2.0d0*H1*psigma(i)-H2
            DpsigmaDH(i,1)=psigma(i)**2.0d0/denom
            DpsigmaDH(i,2)=psigma(i)/denom
            DpsigmaDH(i,3)=2.0d0/3.0d0/denom
          enddo
c
          do i=1,3
            do j=1,3
              DFDH(i)=DFDH(i)+DFDpsigma(j)*DpsigmaDH(j,i)
            enddo
          enddo
        else
c                                                     ---- Singular case
          if ( abs(psigma(2)-psigma(3))/se <= eps) then
c
c                                           ---- Case1 S2=S3 ( theta=0 )
            denom=psigma(1)**2.0d0-2.0d0*H1*psigma(1)-H2
            DpsigmaDH(1,1)=psigma(1)**2.0d0/denom
            DpsigmaDH(1,2)=psigma(1)/denom
            DpsigmaDH(1,3)=2.0d0/3.0d0/denom
c
            DFDH(1)=DpsigmaDH(1,1)*(DFDpsigma(1)-DFDpsigma(2))
     1              +3.0d0*DFDpsigma(2)
            DFDH(2)=DpsigmaDH(1,2)*(DFDpsigma(1)-DFDpsigma(2))
            DFDH(3)=DpsigmaDH(1,3)*(DFDpsigma(1)-DFDpsigma(2))
          elseif ( abs(psigma(2)-psigma(1))/se <= eps ) then
c
c                                          ---- Case2 S2=S1 ( theta=pi )
            denom=psigma(3)**2.0d0-2.0d0*H1*psigma(3)-H2
            DpsigmaDH(3,1)=psigma(3)**2.0d0/denom
            DpsigmaDH(3,2)=psigma(3)/denom
            DpsigmaDH(3,3)=2.0d0/3.0d0/denom
c
            DFDH(1)=DpsigmaDH(3,1)*(DFDpsigma(3)-DFDpsigma(2))
     1              +3.0d0*DFDpsigma(2)
            DFDH(2)=DpsigmaDH(3,2)*(DFDpsigma(3)-DFDpsigma(2))
            DFDH(3)=DpsigmaDH(3,3)*(DFDpsigma(3)-DFDpsigma(2))
          else

          endif
        endif
c
c                                     ---- D(H)/D(sigma) -> 3 x 6 Matrix
c
         call jancae_clear2 ( DHDsigma,3,6 )
c
         DHDsigma(1,1)=1.0d0/3.0d0
         DHDsigma(1,2)=1.0d0/3.0d0
         DHDsigma(1,3)=1.0d0/3.0d0
c
         DHDsigma(2,1)=-1.0d0/3.0d0*(sigma(2)+sigma(3))
         DHDsigma(2,2)=-1.0d0/3.0d0*(sigma(3)+sigma(1))
         DHDsigma(2,3)=-1.0d0/3.0d0*(sigma(1)+sigma(2))
         DHDsigma(2,4)=2.0d0/3.0d0*sigma(4)
         DHDsigma(2,5)=2.0d0/3.0d0*sigma(5)
         DHDsigma(2,6)=2.0d0/3.0d0*sigma(6)
c
c                                            !!-sigma(5)**2.0d0)
         DHDsigma(3,1)=0.5d0*(sigma(2)*sigma(3)-sigma(6)**2.0d0)
c                                            !!-sigma(6)**2.0d0)
         DHDsigma(3,2)=0.5d0*(sigma(3)*sigma(1)-sigma(5)**2.0d0)
         DHDsigma(3,3)=0.5d0*(sigma(1)*sigma(2)-sigma(4)**2.0d0)
         DHDsigma(3,4)=sigma(5)*sigma(6)-sigma(3)*sigma(4)
         DHDsigma(3,6)=sigma(6)*sigma(4)-sigma(1)*sigma(5) !!...(3,5)=
         DHDsigma(3,5)=sigma(4)*sigma(5)-sigma(2)*sigma(6) !!...(3,6)=
c
c                                     ---- D(sigma)/D(s) -> 6 x 6 Matrix
c
        do i=1,6
          do j=1,6
            DsigmaDs(i,j)=ct(i,j)
          end do
        end do
c
c                                        ---- D(se)/D(s) -> 1 x 3 Vector
c
        call jancae_clear1( DFDs,6 )
        call jancae_clear2( dummat,3,6 )
c
        do i=1,6
          do j=1,3
            do k=1,6
              dummat(j,i)=dummat(j,i)+DHDsigma(j,k)*DsigmaDs(k,i)
            enddo
            DFDs(i)=DFDs(i)+DFDH(j)*dummat(j,i)
          enddo
        enddo
      endif
c
      do i=1,6
        dseds(i)=DseDF*DFDs(i)
      enddo
c
c                                            ---- 2nd order differential
c
      if ( nreq >= 2 ) then
c                                              ---- D(se)/D(F) -> Scalar
c
        D2seDF2=(1.0d0/D)**ai*ai*(ai-1.0d0)*F**(ai-2.0d0)
c
c                                  ---- D2(F)/D(psigma)2 -> 3 x 3 Matrix
c
        call jancae_clear2( D2FDpsigma2,3,3 )
c
        do i=1,3
          D2FDpsigma2(i,i)=a*(psigma(i)/abs(psigma(i))-ck)**2.0d0
     1                     *(abs(psigma(i))-ck*psigma(i))**(a-2.0d0)
        enddo
c
c                              ---- D2(psigma)/D(H)2 -> 3 x 3 x 3 Matrix
c
        call jancae_clear3 ( D2psigmaDH2,3,3,3 )
c
        if ( abs(psigma(2)-psigma(3))/se>eps .and.
     1       abs(psigma(2)-psigma(1))/se>eps ) then
c                                                 ---- Not Singular case
c
          do i=1,3
            denom=(psigma(i)**2.0d0-2.0d0*H1*psigma(i)-H2)**3.0d0
c
            D2psigmaDH2(i,1,1)=2.0d0*psigma(i)**3.0d0
     1                         *(psigma(i)**2.0d0-3.0d0*H1*psigma(i)
     2                         -2.0d0*H2)/denom
            D2psigmaDH2(i,2,2)=-2.0d0*psigma(i)*(H1*psigma(i)+H2)/denom
            D2psigmaDH2(i,3,3)=-8.0d0/9.0d0*(psigma(i)-H1)/denom
c
            D2psigmaDH2(i,1,2)=psigma(i)**2.0d0*(psigma(i)**2.0d0
     1                         -4.0d0*H1*psigma(i)-3.0d0*H2)/denom
            D2psigmaDH2(i,2,3)=-2.0d0/3.0d0*(psigma(i)**2.0d0+H2)/denom
            D2psigmaDH2(i,3,1)=-4.0d0/3.0d0*psigma(i)
     1                         *(H1*psigma(i)+H2)/denom
c
            D2psigmaDH2(i,2,1)=D2psigmaDH2(i,1,2)
            D2psigmaDH2(i,3,2)=D2psigmaDH2(i,2,3)
            D2psigmaDH2(i,1,3)=D2psigmaDH2(i,3,1)
          enddo
c
c                                       ---- D2(F)/D(H)2 -> 3 x 3 Matrix
c
          call jancae_clear2 ( D2FDH2,3,3 )
c
          do iq=1,3
            do m=1,3
              do ip=1,3
                do l=1,3
                  D2FDH2(iq,m)=D2FDH2(iq,m)+D2FDpsigma2(ip,l)
     1                         *DpsigmaDH(l,m)*DpsigmaDH(ip,iq)
                enddo
                D2FDH2(iq,m)=D2FDH2(iq,m)+DFDpsigma(ip)
     1                       *D2psigmaDH2(ip,iq,m)
              enddo
            enddo
          enddo
c
c                               ---- D2(H)/D(sigma)2 -> 3 x 6 x 6 Matrix
c
          call jancae_clear3 ( D2HDsigma2,3,6,6 )
c
          D2HDsigma2(2,1,2)=-1.0d0/3.0d0
          D2HDsigma2(2,2,3)=-1.0d0/3.0d0
          D2HDsigma2(2,3,1)=-1.0d0/3.0d0
c
          D2HDsigma2(2,2,1)=D2HDsigma2(2,1,2)
          D2HDsigma2(2,3,2)=D2HDsigma2(2,2,3)
          D2HDsigma2(2,1,3)=D2HDsigma2(2,3,1)
c
          D2HDsigma2(2,4,4)=2.0d0/3.0d0
          D2HDsigma2(2,5,5)=2.0d0/3.0d0
          D2HDsigma2(2,6,6)=2.0d0/3.0d0
c
          D2HDsigma2(3,1,2)=sigma(3)/2.0d0
          D2HDsigma2(3,2,3)=sigma(1)/2.0d0
          D2HDsigma2(3,3,1)=sigma(2)/2.0d0
c
          D2HDsigma2(3,2,1)=D2HDsigma2(3,1,2)
          D2HDsigma2(3,3,2)=D2HDsigma2(3,2,3)
          D2HDsigma2(3,1,3)=D2HDsigma2(3,3,1)
c
          D2HDsigma2(3,4,4)=-sigma(3)
          D2HDsigma2(3,6,6)=-sigma(1)         !!...(3,5,5)
          D2HDsigma2(3,5,5)=-sigma(2)         !!...(3,6,6)
c
          D2HDsigma2(3,4,5)=sigma(6)
          D2HDsigma2(3,5,6)=sigma(4)
          D2HDsigma2(3,6,4)=sigma(5)
c
          D2HDsigma2(3,5,4)=D2HDsigma2(3,4,5)
          D2HDsigma2(3,6,5)=D2HDsigma2(3,5,6)
          D2HDsigma2(3,4,6)=D2HDsigma2(3,6,4)
c
          D2HDsigma2(3,1,6)=-sigma(6)         !!...(3,1,5)=-sigma(5)
          D2HDsigma2(3,6,1)=D2HDsigma2(3,1,6) !!...(3,5,1)=...(3,1,5)
c
          D2HDsigma2(3,2,5)=-sigma(5)         !!...(3,2,6)=-sigma(6)
          D2HDsigma2(3,5,2)=D2HDsigma2(3,2,5) !!...(3,6,2)=...(3,2,6)
c
          D2HDsigma2(3,3,4)=-sigma(4)
          D2HDsigma2(3,4,3)=D2HDsigma2(3,3,4)
c
c                                       ---- D2(F)/D(s)2 -> 6 x 6 Matrix
c
          call jancae_clear2( D2FDs2,6,6 )
          call jancae_clear2( dummat,3,6 )
c
          do i=1,3
            do j=1,6
              do ip=1,6
                dummat(i,j)=dummat(i,j)+DHDsigma(i,ip)*DsigmaDs(ip,j)
              end do
            end do
          end do
c
          do i=1,6
            do j=1,6
              do iq=1,3
                do m=1,3
                  D2FDs2(i,j)=D2FDs2(i,j)+D2FDH2(iq,m)
     1                        *dummat(iq,i)*dummat(m,j)
                end do
c
                do n=1,6
                  do ir=1,6
                    D2FDs2(i,j)=D2FDs2(i,j)+D2HDsigma2(iq,ir,n)
     1                          *DFDH(iq)*DsigmaDs(ir,i)*DsigmaDs(n,j)
                  end do
                end do
              end do
            end do
          end do
c
c                                      ---- D2(se)/D(s)2 -> 6 x 6 Matrix
c
          do i=1,6
            do j=1,6
              d2seds2(i,j)=D2seDF2*DfDs(i)*DfDs(j)+DseDF*D2FDs2(i,j)
            end do
          end do
        else
c
c                                                     ---- Singular case
c
          del= eps
c
          do i=1,6
            s0(i)= s(i)
          enddo
c
          do i=1,6
            do j=1,6
              if (i == j) then
                s0(i)=s(i)-del
                sea= cazacu2006_seND(s0,ct,phi,ck,a,ai)
                s0(i)=s(i)+del
                seb= cazacu2006_seND(s0,ct,phi,ck,a,ai)
c
                s0(i)=s(i)
                abc1=(se-sea)/del
                abc2=(seb-se)/del
                d2seds2(i,j) = (abc2-abc1)/del
              else
                s0(i)=s(i)-del
                s0(j)=s(j)-del
                seaa= cazacu2006_seND(s0,ct,phi,ck,a,ai)
c
                s0(i)=s(i)+del
                s0(j)=s(j)-del
                seba= cazacu2006_seND(s0,ct,phi,ck,a,ai)
c
                s0(i)=s(i)-del
                s0(j)=s(j)+del
                seab= cazacu2006_seND(s0,ct,phi,ck,a,ai)
c
                s0(i)=s(i)+del
                s0(j)=s(j)+del
                sebb= cazacu2006_seND(s0,ct,phi,ck,a,ai)
c
                s0(i)=s(i)
                s0(j)=s(j)
                abc1 = (seba - seaa)/(2.0d0*del)
                abc2 = (sebb - seab)/(2.0d0*del)
                d2seds2(i,j) = (abc2 - abc1)/(2.0d0*del)
              end if
            end do
          end do
        end if
      endif
c
      return
      end
c
c
c
c-----------------------------------------------------------(cazacu2006)
      double precision function cazacu2006_seND(s,ct,phi,ck,a,ai)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension s(6),sigma(6),psigma(3)
      dimension ct(6,6), phi(3)
      parameter( pi = 3.141592653589793d0 )
c
c                               ---- Calculate linear transformed stress
c
      call jancae_mv( sigma,ct,s,6,6 )
c                                       ---- 1st, 2nd and 3rd invariants
      H1=(sigma(1)+sigma(2)+sigma(3))/3.0d0
      H2=(sigma(5)**2.0d0+sigma(6)**2.0d0+sigma(4)**2.0d0
     1    -sigma(2)*sigma(3)-sigma(3)*sigma(1)-sigma(1)*sigma(2))/3.0d0
      H3=(2.0d0*sigma(5)*sigma(6)*sigma(4)+sigma(1)*sigma(2)*sigma(3)
     2    -sigma(1)*sigma(6)**2.0d0-sigma(2)*sigma(5)**2.0d0
     3    -sigma(3)*sigma(4)**2.0d0)/2.0d0 ! sigma(5) <-> sigma(6)
c
      p=H1**2.0d0+H2
      q=(2.0d0*H1**3.0d0+3.0d0*H1*H2+2.0d0*H3)/2.0d0
      theta=q/(p**1.5d0)
      if ( theta >  1.0d0 ) theta= 1.0d0
      if ( theta < -1.0d0 ) theta=-1.0d0
      theta=acos(theta)
c                               ---- Calculate principal values of sigma
      psigma(1)=2.0d0*sqrt(p)*cos(theta/3.0d0)+H1
      psigma(2)=2.0d0*sqrt(p)*cos((theta+4.0d0*pi)/3.0d0)+H1
      psigma(3)=2.0d0*sqrt(p)*cos((theta+2.0d0*pi)/3.0d0)+H1
c
c                                                 ---- Equivalent stress
c                                          ---- Calculate yield function
      F = (abs(psigma(1))-ck*psigma(1))**a
     1   +(abs(psigma(2))-ck*psigma(2))**a
     2   +(abs(psigma(3))-ck*psigma(3))**a
c                                           ---- Denominator coefficient
      D = (abs(phi(1))-ck*phi(1))**a
     1   +(abs(phi(2))-ck*phi(2))**a
     2   +(abs(phi(3))-ck*phi(3))**a
c
      cazacu2006_seND=(F/D)**ai
c
      return
      end
c
c
c
c***********************************************************************
c     JANCAE/UMMDp : Yield Criteria
c***********************************************************************
c
c      0 : von Mises isotropic (1913)
c
c      1 : Hill quadratic (1948)
c      2 : Barlat yld2004 (2005)
c      3 : Cazacu (2006)
c      4 : Karafillis-Boyce (1993)
c      5 : Hu (2005)
c      6 : Yohsida (2011)
c
c     -1 : Gotoh biquadratic (1978)
c     -2 : Barlat YLD2000-2d (2000)
c     -3 : Vegter
c     -4 : Banabic BBC2005
c     -5 : Barlat YLD89
c     -6 : Banabic BBC2008
c     -7 : Hill 1990
c
c-----------------------------------------------------------------------
c     yield function and its dfferentials
c
      subroutine jancae_yfunc  ( se,cdseds,cd2seds2,nreq,
     &                           cs,nttl,nnrm,nshr,
     &                           pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension cs(nttl),cdseds(nttl),cd2seds2(nttl,nttl),
     &          pryld(ndyld)
      dimension s(6),dseds(6),d2seds2(6,6),indx(6)
c
      ntyld=nint(pryld(1))
c
      if ( ntyld.lt.0 ) then
        if ( (nnrm.ne.2).or.(nshr.ne.1) ) then
          write (6,*) 'error in jancae_yfunc'
          write (6,*) 'ntyld<0 for plane stress'
          write (6,*) 'nnrm,nshr,ntyld:',nnrm,nshr,ntyld
          call jancae_exit (9000)
        endif
        goto 100
      endif
c
      ss=0.0
      do i=1,nttl
        ss=ss+cs(i)**2
      enddo
      if ( (ss.le.0.0).and.(nreq.eq.0) ) then
        se=0.0
        return
      endif
c
c                                                ---- 3D yield functions
c
c                                        ---- set index to s(i) to cs(i)
      do i=1,6
        indx(i)=0
      enddo
      if ( nnrm.eq.3 ) then
        do i=1,nttl
          indx(i)=i
        enddo
      else if ( nnrm.eq.2 ) then
        indx(1)=1
        indx(2)=2
        indx(3)=0
        do i=1,nshr
          indx(3+i)=2+i
        enddo
      endif
c                                                          ---- set s(i)
      call jancae_clear1 ( s,6 )
      do i=1,6
        if ( indx(i).ne.0 ) then
          s(i)=cs(indx(i))
        endif
      enddo
c
      select case ( ntyld )
      case ( 0 )                                      ! von Mises (1913)
        call jancae_mises ( s,se,dseds,d2seds2,nreq )
      case ( 1 )                                !  Hill quadratic (1948)
        call jancae_hill_1948 ( s,se,dseds,d2seds2,nreq,
     &                          pryld,ndyld )
      case ( 2 )                                 ! Barlat yld2004 (2005)
        call jancae_yld2004_18p ( s,se,dseds,d2seds2,nreq,
     &                            pryld,ndyld )
      case ( 3 )                                         ! Cazacu (2006)
        call jancae_cazacu2006 ( s,se,dseds,d2seds2,nreq,
     &                           pryld,ndyld )
      case ( 4 )                               ! Karafillis-Boyce (1993)
        call jancae_KarafillisBoyce ( s,se,dseds,d2seds2,nreq,
     &                                pryld,ndyld )
      case ( 5 )                                             ! Hu (2005)
        call jancae_hu_2005 ( s,se,dseds,d2seds2,nreq,
     &                       pryld,ndyld )
      case ( 6 )                                      ! F.Yoshida (2011)
        call jancae_yoshida_2011 ( s,se,dseds,d2seds2,nreq,
     &                             pryld,ndyld )
c
      case default
        write (6,*) 'error in jancae_yfunc'
        write (6,*) 'ntyld error :',ntyld
        call jancae_exit (9000)
      end select
c
c                                                        ---- set dse/ds
      if ( nreq.ge.1 ) then
        do i=1,6
          if ( indx(i).ne.0 ) cdseds(indx(i))=dseds(i)
        enddo
      endif
c                                                      ---- set d2se/ds2
      if ( nreq.ge.2 ) then
        do i=1,6
          if ( indx(i).ne.0 ) then
            do j=1,6
              if ( indx(j).ne.0 ) then
                cd2seds2(indx(i),indx(j))=d2seds2(i,j)
              endif
            enddo
          endif
        enddo
      endif
c
      return
c
c
  100 continue
c                                      ---- plane stress yield functions
c
      select case ( ntyld )
      case ( -1 )                             ! Gotoh biquadratic (1978)
        call jancae_gotoh ( cs,se,cdseds,cd2seds2,nreq,
     &                      pryld,ndyld )
      case ( -2 )                             ! Barlat yld2000-2d (2000)
        call jancae_yld2000 ( cs,se,cdseds,cd2seds2,nreq,
     &                        pryld,ndyld )
      case ( -3 )                                        ! Vegter (2006)
        call jancae_vegter ( cs,se,cdseds,cd2seds2,nreq,
     &                       pryld,ndyld )
      case ( -4 )                                      ! Banabic BBC2005
        call jancae_bbc2005 ( cs,se,cdseds,cd2seds2,nreq,
     &                        pryld,ndyld )
      case ( -5 )                                         ! Barlat yld89
        call jancae_yld89 ( cs,se,cdseds,cd2seds2,nreq,
     &                      pryld,ndyld )
      case ( -6 )                                      ! Banabic BBC2008
        call jancae_bbc2008 ( cs,se,cdseds,cd2seds2,nreq,
     &                        pryld,ndyld )
      case ( -7 )                                            ! Hill 1990
        call jancae_hill90  ( cs,se,cdseds,cd2seds2,nreq,
     &                        pryld,ndyld )

      case default
        write (6,*) 'error in jancae_yfunc'
        write (6,*) 'ntyld error :',ntyld
        call jancae_exit (9000)
      end select
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     print type and parameters for yield functions
c
      subroutine jancae_yfunc_print ( pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension pryld(ndyld)
c
      ntyld=pryld(1)
      write (6,*)
      write (6,*) '*** Yield Criterion',ntyld
      select case ( ntyld )
      case ( 0 )
        write (6,*) 'von Mises'
      case ( 1 )
        write (6,*) 'Hill 1948'
        write (6,*) 'F=',pryld(2)
        write (6,*) 'G=',pryld(3)
        write (6,*) 'H=',pryld(4)
        write (6,*) 'L=',pryld(5)
        write (6,*) 'M=',pryld(6)
        write (6,*) 'N=',pryld(7)
      case ( 2 )
        write (6,*) 'Barlat yld2004-18p'
        n0=1
        do i=1,18
           n0=n0+1
           write (6,*) 'a(',i,')=',pryld(n0)
        enddo
        write (6,*) 'M=',pryld(1+18+1)
      case ( 3 )
        write (6,*) 'Cazacu 2006'
        n0=1
        do i=1,3
          do j=1,3
            n0=n0+1 ; write (6,*) 'c(',i,',',j,')=',pryld(n0)
          enddo
        enddo
        do i=4,6
          n0=n0+1 ; write (6,*) 'c(',i,',',i,')=',pryld(n0)
        enddo
        n0=n0+1 ; write (6,*) 'a =',pryld(n0)
        n0=n0+1 ; write (6,*) 'ck=',pryld(n0)
      case ( 4 )
        write (6,*) 'Karafillis-Boyce'
        n0=1
        do i=1,6
          do j=i,6
            n0=n0+1
            write (6,*) 'L(',i,',',j,') =',pryld(n0)
          end do
        end do
        n0=n0+1 ; write (6,*) 'k =',pryld(n0)
        n0=n0+1 ; write (6,*) 'c =',pryld(n0)
      case ( 5 )
        write (6,*) 'Weilong Hu 4th order (2005)'
        n0=1
        do i=1,5
          n0=n0+1
          write (6,*) 'X(',i,')=',pryld(n0)
        enddo
        n0=n0+1 ; write (6,*) 'X(',7,')=',pryld(n0)
        do i=1,3
          n0=n0+1
          write (6,*) 'C(',i,')=',pryld(n0)
        enddo
      case ( 6 )
        write (6,*) 'F.Yoshida 6th order (2011)'
        n0=1
        do i=1,16
          n0=n0+1
          write (6,*) 'c(',i,')=',pryld(n0)
        enddo
c
      case ( -1 )
        write (6,*) 'Gotoh biquadratic'
        do i=1,9
          write (6,*) 'A(',i,')=',pryld(i+1)
        enddo
      case ( -2 )
        write (6,*) 'Barlat YLD2000-2d'
        do i=1,8
          write (6,*) 'a(',i,')=',pryld(i+1)
        enddo
        write (6,*) 'M=',pryld(9+1)
      case ( -3 )
        write (6,*) 'Vegter '
        write (6,*) 'nf=',nint(pryld(2))
        write (6,*) 'f_bi0=',pryld(3)
        write (6,*) 'r_bi0=',pryld(4)
        do i=0,nint(pryld(2))
          write (6,*) 'test angle=',90.0d0*float(i)/pryld(2)
          write (6,*) 'phi_un(',i,')=',pryld(4+i*4+1)
          write (6,*) 'phi_sh(',i,')=',pryld(4+i*4+2)
          write (6,*) 'phi_ps(',i,')=',pryld(4+i*4+3)
          write (6,*) 'omg(   ',i,')=',pryld(4+i*4+4)
        enddo
c       do i=1,7
c         write (6,*) 'phi_un(',i-1,')=',pryld(1+i   )
c         write (6,*) 'phi_sh(',i-1,')=',pryld(1+i+ 7)
c         write (6,*) 'phi_ps(',i-1,')=',pryld(1+i+14)
c         write (6,*) 'omg   (',i-1,')=',pryld(1+i+23)
c       enddo
c       write (6,*)   'f_bi0=',pryld(1+22)
c       write (6,*)   'r_bi0=',pryld(1+23)
c       write (6,*)   'nf   =',nint(pryld(1+31))
      case ( -4 )
        write (6,*) 'BBC2005'
        write (6,*) 'k of order 2k',pryld(1+1)
        write (6,*) 'a=',pryld(1+2)
        write (6,*) 'b=',pryld(1+3)
        write (6,*) 'L=',pryld(1+4)
        write (6,*) 'M=',pryld(1+5)
        write (6,*) 'N=',pryld(1+6)
        write (6,*) 'P=',pryld(1+7)
        write (6,*) 'Q=',pryld(1+8)
        write (6,*) 'R=',pryld(1+9)
      case ( -5 )
        write (6,*) 'YLD89'
        write (6,*) 'order M=',pryld(1+1)
        write (6,*) 'a      =',pryld(1+2)
        write (6,*) 'h      =',pryld(1+3)
        write (6,*) 'p      =',pryld(1+4)
      case ( -6 )
        write (6,*) 'BBC2008'
        write (6,*) 's      =',nint(pryld(1+1))
        write (6,*) 'k      =',nint(pryld(1+2))
        do i=1,nint(pryld(1+1))
          write (6,*) 'i=',i
          n=2+(i-1)*8
          write (6,*) 'l_1=',pryld(n+1)
          write (6,*) 'l_2=',pryld(n+2)
          write (6,*) 'm_1=',pryld(n+3)
          write (6,*) 'm_2=',pryld(n+4)
          write (6,*) 'm_3=',pryld(n+5)
          write (6,*) 'n_1=',pryld(n+6)
          write (6,*) 'n_2=',pryld(n+7)
          write (6,*) 'n_3=',pryld(n+8)
        enddo
      case ( -7 )
        write (6,*) 'Hill90'
        write (6,*) 'a   =',pryld(1+1)
        write (6,*) 'b   =',pryld(1+2)
        write (6,*) 'tau =',pryld(1+3)
        write (6,*) 'sigb=',pryld(1+4)
        write (6,*) 'M   =',pryld(1+5)
      end select
c
      return
      end
c
c
c
c----------------------------------------------------------------(gotoh)
c     Gotoh biquadratic yield function and its dfferentials
c     (J.JSTP vol.19 no.208 p377-385 (1978-5) )
c
      subroutine jancae_gotoh ( s,se,dseds,d2seds2,nreq,
     &                          pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension s(3),dseds(3),d2seds2(3,3),pryld(ndyld)
      dimension a(9),c(4,4),v(4),
     &          t(4),dtds(4,3),d2tds2(4,3,3)
c
c     a(i)      : coef.s of Gotoh's 4th order function
c     t(i)      : stress^2 vector
c                 ={ sx^2, sx*sy, sy^2, txy^2 }
c     c(i,j)    : matrix to calc. se
c              Gotoh's function se = ( {t}^T*[c]*{t} )^(1/4)
c
c     dtds(i,j) : diff. of t(i) with respect to s(j)
c     d2tds2(i,j,k)
c               : 2nd order diff. of t(i) w.r.t s(j) & s(k)
c
c                                            ---- anisotropic parameters
      do i=1,9
        a(i)=pryld(1+i)
      enddo
c                                                      ---- coef. matrix
      c(1,1)=a(1)
      c(1,2)=a(2)*0.5d0
      c(1,3)=0.0
      c(1,4)=a(6)*0.5d0
      c(2,2)=a(3)
      c(2,3)=a(4)*0.5d0
      c(2,4)=a(7)*0.5d0
      c(3,3)=a(5)
      c(3,4)=a(8)*0.5d0
      c(4,4)=a(9)
      do i=2,4
        do j=1,i-1
          c(i,j)=c(j,i)
        enddo
      enddo
c                                                    ---- t-vector (s^2)
      t(1)=s(1)*s(1)
      t(2)=s(1)*s(2)
      t(3)=s(2)*s(2)
      t(4)=s(3)*s(3)
c                                                 ---- equivalent stress
      call jancae_mv  ( v,c,t,4,4 )
      call jancae_vvs ( phi,t,v,4 )
c
      if ( phi.le.0.0 ) phi=0.0
      se=sqrt(sqrt(phi))
c                                            ---- 1st order differential
      if ( nreq.ge.1 ) then
        call jancae_clear2 ( dtds,4,3 )
        dtds(1,1)=s(1)*2.0d0
        dtds(2,1)=s(2)
        dtds(2,2)=s(1)
        dtds(3,2)=s(2)*2.0d0
        dtds(4,3)=s(3)*2.0d0
        call jancae_clear1 ( v,4 )
        do i=1,3
          do j=1,4
            do k=1,4
              v(i)=v(i)+2.0d0*t(j)*c(j,k)*dtds(k,i)
            enddo
          enddo
        enddo
        q=0.25d0*phi**(-0.75d0)
        do i=1,3
          dseds(i)=q*v(i)
        enddo
      endif
c                                            ---- 2nd order differential
      if ( nreq.ge.2 ) then
        call jancae_clear3 ( d2tds2,4,3,3 )
        d2tds2(1,1,1)=2.0d0
        d2tds2(2,1,2)=1.0d0
        d2tds2(2,2,1)=1.0d0
        d2tds2(3,2,2)=2.0d0
        d2tds2(4,3,3)=2.0d0
        call jancae_clear2 ( d2seds2,3,3 )
        do i=1,3
          do j=1,3
            do m=1,4
              do n=1,4
                d2seds2(i,j)=d2seds2(i,j)+
     &                       2.0d0*c(m          ,n       )*
     &                        ( dtds(m,i)*dtds(  n  ,j)+
     &                           t(  m)  *d2tds2(n,i,j)  )
              enddo
            enddo
          enddo
        enddo
        do i=1,3
          do j=1,3
            d2seds2(i,j)=q*(d2seds2(  i,   j)
     &                      -0.75d0*v(i)*v(j)/phi)
          enddo
        enddo
      endif
c
      return
      end
c
c
c
cc---------------------------------------------------------------(hill48)
c     Hill(1948) yield function and its dfferentials
c     (  Proc. Roy. Soc. A193(1948) p281-297 )
c
c     ( flow curve must be defined in uniaxial sx vs ex )
c
      subroutine jancae_hill_1948 ( s,se,dseds,d2seds2,nreq,
     &                              pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension s(6),dseds(6),d2seds2(6,6),pryld(ndyld)
      dimension c(6,6),v(6)
c
c     c(i,j)    : matrix to calc. se
c     Hill's 1948 function se = ({s}^T*[c]*{s})^(1/2)
c
c                                            ---- anisotropic parameters
c                              pf means "F", pg means "G"....
c                         F=G=H=1 and L=M=N=3 means von Mises
      pf=pryld(1+1)
      pg=pryld(1+2)
      ph=pryld(1+3)
      pl=pryld(1+4)
      pm=pryld(1+5)
      pn=pryld(1+6)
c                                                      ---- coef. matrix
      call jancae_clear2 ( c,6,6 )
      c(1,1)=    pg+ph
      c(1,2)=      -ph
      c(1,3)=   -pg
      c(2,1)=      -ph
      c(2,2)= pf   +ph
      c(2,3)=-pf
      c(3,1)=   -pg
      c(3,2)=-pf
      c(3,3)= pf+pg
      c(4,4)=2.0d0*pn
      c(5,5)=2.0d0*pm
      c(6,6)=2.0d0*pl
      do i=1,6
        do j=1,6
          c(i,j)=c(i,j)/(pg+ph)
        enddo
      enddo
c
      call jancae_mv  ( v,c,s,6,6 )
      call jancae_vvs ( phi,s,v,6 )
c                                                 ---- equivalent stress
      if ( phi.le.0.0 ) phi=0.0
      se=sqrt(phi)
c                                            ---- 1st order differential
      if ( nreq.ge.1 ) then
        do i=1,6
          dseds(i)=v(i)/se
        enddo
      endif
c                                            ---- 2nd order differential
      if ( nreq.ge.2 ) then
        do i=1,6
          do j=1,6
            d2seds2(i,j)=( -v(i)*v(j)/phi+c(i,j) )/se
          enddo
        enddo
      endif
c
      return
      end
c
c
c
c---------------------------------------------------------------(hill90)
c     Hill-1990 anisotropic yield function and its dfferentials
c
c     ---( 1990 ) Pergamon Press Plc.
c      by R.Hill (Department of Applied Mathematics and Theoretical
c                 Physics, University of Cambridge, Cambridge) EW, U.K)
c
c     "Constitutive Modelling of Orthotropic Plasticity in Sheet Metals"
c      Received in 18 July 1989)
c      J. Mech. Physics. Solids Vol.38, No.3, pp-405-417,  1990
c      Printed in Great Britain
c-----------------------------------------------------------------------
c     coded by Tatsuhiko Ine ( at NIED ), 22/11/2010
c-----------------------------------------------------------------------
c
      subroutine jancae_hill90 ( s,se,dseds,d2seds2,nreq,
     &                           pryld,ndyld )
c-----------------------------------------------------------------------
c                                                   ---- input arguments
c
c     s(3) :stress 6-components -> in this program only use 3 stress
c           components but we define that stress array has 6 components.
c     pryld(ndyld)  : material parameter for yield function Hill's 1990
c     pryld(1+1)= a
c     pryld(1+2)= b
c     pryld(1+3)= tau
c     pryld(1+4)= sigb
c     pryld(1+5)= m  (=>am)
c     nreq : request code = 0 ( output equivalent stress only)
c          : request code = 1 ( output equivalent stress=se , dseds(3) )
c          : request code = 2 ( output equivalent stress=se , dseds(3),
c                           d2seds2(3,3) )
c
c-----------------------------------------------------------------------
c                                                  ---- output arguments
c
c     se       : equivalent stress
c     dseds(3) : differential coefficient of first order
c              : equivalent stress differentiated by stress (6 components)
c     d2seds2(3,3) : differential coefficient of second order
c              : equivalent stress differentiated by stress (6 components)
c-----------------------------------------------------------------------
c
c     a1(3)        : vector to calc for equivalent stress
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
c     d2fds1(3,3) : d2f1/ds2 =fai1 2nd order derivative by stress
c     d2fds2(3,3) : d2f2/ds2 =fai2 2nd order derivative by stress
c     d2fds3(3,3) : d2f3/ds2 =fai3 2nd order derivative by stress
c     d2fds4(3,3) : d2f4/ds2 =fai4 2nd order derivative by stress
c     d2fds_t(3,3) : d2f/ds2 =fai  2nd order derivative by stress
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension s(3),dseds(3),d2seds2(3,3),pryld(ndyld)
      dimension c(3,3),v(3)
c
      dimension A1(3),A2(3,3),A3(3,3)
     &         ,A4(3,3)
c
      dimension dfds1(3),dfds2(3),dfds3(3)
     &         ,dfds4(3),dfds_t(3)
      dimension dxds1(3),dxds2(3),dxds3(3)
     &         ,dxds4(3)
c
      dimension d2fds1(3,3),d2fds2(3,3)
     &         ,d2fds3(3,3)
     &         ,d2fds4(3,3)
     &         ,d2fds_t(3,3)
c
      dimension dx1dx1(3,3),dx2dx2(3,3)
     &         ,dx3dx3(3,3)
     &         ,dx4dx4(3,3)
c
      dimension df4df3(3,3),df3df4(3,3)
c
      dimension d2xds1(3,3),d2xds2(3,3)
     &         ,d2xds3(3,3)
     &         ,d2xds4(3,3)
      character text*32
c
c                                      ---- define a1-matrix, a2-matrix
      data a1/ 1.0d0 , 1.0d0 , 0.0d0 /,
     &     a2/ 1.0d0 ,-1.0d0 , 0.0d0 ,
     &        -1.0d0 , 1.0d0 , 0.d00 ,
     &         0.0d0 , 0.0d0 , 4.0d0 /,
     &     a3/ 1.0d0 , 0.0d0 , 0.0d0 ,
     &         0.0d0 , 1.0d0 , 0.0d0 ,
     &         0.0d0 , 0.0d0 , 2.0d0 /
c
c                                            ---- anisotropic parameters
       a      =  pryld(1+1)
       b      =  pryld(1+2)
       tau    =  pryld(1+3)
       sigb   =  pryld(1+4)
       am      = pryld(1+5)
c
       syini   =  1.0d0
       sigbtm =(sigb/tau)**am
       alarge = 1.0d0 + sigbtm -2.0d0*a + b
c
c
c                               ---- coef. matrix of material parameters
c                                ---- define a4-matrix consists of a & b
       call jancae_clear2( a4,3,3 )
       a4(1,1) = -2.0d0*a +  b
       a4(2,2) =  2.0d0*a +  b
       a4(1,2) =          - b
       a4(2,1) =          - b
c
c                                                              ---- fai1
      x1=s(1)+s(2)
      fai1=(abs( x1 ))**am
c
c                                                              ---- fai2
      call jancae_mv  ( v,a2,s,3,3 )
      call jancae_vvs ( x2,s,v,3 )
      fai2= sigbtm * ( x2 )**(am/2.0d0)
c
c                                                              ---- fai3
      call jancae_mv  ( v,a3,s,3,3 )
      call jancae_vvs ( x3,s,v,3 )
      fai3=  ( x3 )**(am/2.0d0-1.0d0)
c
c                                                              ---- fai4
      call jancae_mv  ( v,a4,s,3,3 )
      call jancae_vvs ( x4,s,v,3 )
      fai4=  ( x4 )
c
c                                             ---- yield fuction : fyild
      fyild= fai1 + fai2 + fai3*fai4
c
c                                                 ---- equivalent stress
      se= (fyild/alarge)**(1.0d0/am)
c
      if ( nreq.eq.0 ) return

c               ---- 1st order differential coefficient of yield fuction
c     dfdsi(i) : diff. of fai-i(i) with respect to s(j)
c     dxdsi(i) : diff. of x-i(i) with respect to s(j)
c                         1st order differential of x_number_i
c
c
c                                                          ---- dfai1/ds
      dxds1(1)=1.0
      dxds1(2)=1.0
      dxds1(3)=0.0
c
      wrk= am * ( abs(x1)**(am-2) )*x1
        do i=1,3
         dfds1(i)=wrk*dxds1(i)
        enddo
c
c                                                          ---- dfai2/ds
      wrk= sigbtm * (am/2.0) * ( x2 )**(am/2.0-1.0)
      call jancae_mv( dxds2,a2,s,3,3 )
        do i=1,3
         dxds2(i) = 2.0 *  dxds2(i)
         dfds2(i) = wrk * dxds2(i)
        enddo
c
c                                                          ---- dfai3/ds
      wrk= (am/2.0-1.0) * ( x3 )**(am/2.0-2.0)
      call jancae_mv( dxds3,a3,s,3,3 )
        do i=1,3
         dxds3(i) = 2.0 * dxds3(i)
         dfds3(i) = wrk * dxds3(i)
        enddo
c
c                                                          ---- dfai4/ds
      call jancae_mv( dxds4,a4,s,3,3 )
        do i=1,3
         dxds4(i)=  2.0 * dxds4(i)
         dfds4(i)=  dxds4(i)
        enddo
c
c
c        ---- 1st order differential coefficient of yield fuction result
c                                     ---- dfai/ds()= result = dfds_t(i)
c
        do i=1,3
          dfds_t(i)=dfds1(i)+dfds2(i)+dfds3(i)*fai4+fai3*dfds4(i)
        enddo
c
c           ---- 1st order differential coefficient of equivalent stress
c
        wrk= (abs(fyild/alarge))**(1.0/am-1.0) / ( am*alarge )
        do i=1,3
          dseds(i)= wrk * dfds_t(i)
        enddo
c
c
        if ( nreq.eq.1 ) return
c
c            --- 2st order differential coefficient of equivalent stress
c                                                   with respect to s(j)
c-----------------------------------------------------------------------
c     dfds_t(3,3) : 2nd order differ of fai by s(j) & s(k)
c     df2ds1(3,3) : 2nd order diff. of  s(j) & s(k)
c-----------------------------------------------------------------------
c
c                                                        ---- d2fai1/ds2
       wrk = am*(am-1.0)*( abs(x1) )**(am-2.0)
        do i=1,3
         do j=1,3
          d2fds1(i,j) = wrk * dxds1(i)*dxds1(j)
         enddo
        enddo
c
c                                                        ---- d2fai2/ds2
      wrk1= sigbtm * (am/2.0)
        if( abs(x2).lt.1e-10) then
           x2=1e-10
        end if
      wrk2=(am/2.0-1.0) * ( x2**(am/2.0-2.0) )
      wrk3=  x2**(am/2.0-1.0)
      wrk2 = wrk1 * wrk2
      wrk3 = wrk1 * wrk3
c
c             ---- make [ dx2 * dx2(t) ] & [d2x/ds2] & make [d2fai2/ds2]
        do i=1,3
         do j=1,3
         dx2dx2(i,j)=dxds2(i)*dxds2(j)
         d2xds2(i,j)=2.0*a2(j,i)
         d2fds2(i,j)= wrk2 * dx2dx2(i,j) + wrk3 * d2xds2(i,j)
         enddo
        enddo
c
c
c                                   ---- d2fai3/ds2   make   d2fds3(i,j)
      wrk1 =(am/2.0-1.0)
      wrk2 =(am/2.0-2.0) * ( x3**(am/2.0-3.0) )
      wrk3 = x3**(am/2.0-2.0)
      wrk2 =  wrk1 * wrk2
      wrk3 =  wrk1 * wrk3
c
c                                   ---- [d2x3/ds2] &  make [d2fai3/ds2]
        do i=1,3
         do j=1,3
          dx3dx3(i,j)= dxds3(i)*dxds3(j)
          d2xds3(i,j)= 2.0 * a3(j,i)
          d2fds3(i,j)= wrk2 * dx3dx3(i,j) + wrk3 * d2xds3(i,j)
         enddo
        enddo
c
c                                                 ---- [d2fai3/ds2]*fai4
        do i=1,3
         do j=1,3
          d2fds3(i,j)= d2fds3(i,j) * fai4
         enddo
        enddo
c
c                                          ---- [dfai4/ds]*[dfai3/ds](T)
        do i=1,3
         do j=1,3
          df4df3(i,j)= dfds4(i)*dfds3(j)
         enddo
        enddo
c
c
c                                                        ---- d2fai4/ds2
c                                                 ---- make [d2fai3/ds2]
        do i=1,3
         do j=1,3
          d2fds4(i,j)=  2.0*a4(i,j)
         enddo
        enddo
c
c        ---- 2nd order differential coefficient of yield fuction result
c                                  ---- d2fai/ds2()= result = d2fds_t(i)
c
        do i=1,3
         do j=1,3
           d2fds_t(i,j) =    d2fds1(i,j)  + d2fds2(i,j)
     &               +  d2fds3(i,j)*fai4  + df4df3(i,j)
     &               +  df4df3(j,i)       + fai3 * d2fds4(i,j)
         enddo
        enddo
c
c
c            ---- 2n order differential coefficient of equivalent stress
c                                                              by stress
c
        wrk1 = 1.0/(am*alarge)
        wrk2 = (1.0/am-1.0)/alarge
        wrk3 = (fyild/alarge)**(1.0/am-2.0)
        wrk4 = (fyild/alarge)**(1.0/am-1.0)
        wrk2 =  wrk1  * wrk2* wrk3
        wrk4 =  wrk1  * wrk4
c
        do i=1,3
         do j=1,3
          d2seds2(i,j)=  wrk2 * dfds_t(i)* dfds_t(j)
     &                 + wrk4 * d2fds_t(i,j)
         enddo
        enddo
c
c
c
      return
      end
c
c
c
c---------------------------------------------------------------(hu2005)
c     WeiLong Hu (2005) yield function and its dfferentials
c
c     "An orthotropic yield criterion in a 3-D general stress state"
c      IJP,v.21(2005), pp.1771-1796. )
c
      subroutine jancae_hu_2005 ( s,se,dseds,d2seds2,nreq,
     &                            pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (maxa=100)
      dimension s(6),dseds(6),d2seds2(6,6),pryld(ndyld)
      dimension a(maxa),ipow(maxa,3)
c
c
      nd0=2
c
      n=0
      do it=0,nd0
        n=n+(nd0-it)*2+1
      enddo
      nterms=n
      if ( maxa.lt.nterms ) then
        write (6,*) 'increase maxa :',maxa,nterms
        call jancae_exit ( 9000 )
      endif
c
      n=0
      ipow=0
      do it=0,nd0
        ndmax=nd0*2-it*2
        do jy=0,ndmax
          jx=ndmax-jy
          n=n+1
          ipow(n,1)=jx
          ipow(n,2)=jy
          ipow(n,3)=it
        enddo
      enddo
c
      a     = 0.0
      a( 1)= 1.0d0*pryld(1+ 1)    ! X1
      a( 2)= 1.0d0*pryld(1+ 2)    ! X2
      a( 3)= 1.0d0*pryld(1+ 3)    ! X3
      a( 4)= 1.0d0*pryld(1+ 4)    ! X4
      a( 5)= 1.0d0*pryld(1+ 5)    ! X5
      a( 6)= 1.0d0*pryld(1+ 7)    ! C1 <-
      a( 7)=-1.0d0*pryld(1+ 9)    ! C3 <-
      a( 8)= 1.0d0*pryld(1+ 8)    ! C2 <-
      a( 9)= 1.0d0*pryld(1+ 6)    ! X7
c
      call jancae_hy_polytype ( s,se,dseds,d2seds2,nreq,
     &                          nd0,a,ipow,maxa,nterms )
c
      return
      end
c
c
cc------------------------------------------------------(karafillisboyce)
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
cc----------------------------------------------------------------(mises)
c     von Mises yield function and its dfferentials
c     ( 1913 )
c
      subroutine jancae_mises ( s,se,dseds,d2seds2,nreq )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension s(6),dseds(6),d2seds2(6,6)
      dimension c(6,6),v(6)
c                                                      ---- coef. matrix
      call jancae_clear2 ( c,6,6 )
      do i=1,3
        do j=1,3
          c(i,j)=-0.5d0
        enddo
        c(i,i)=1.0d0
      enddo
      do i=4,6
        c(i,i)=3.0d0
      enddo
c
      call jancae_mv  ( v,c,s,6,6 )
      call jancae_vvs ( phi,s,v,6 )
c                                                 ---- equivalent stress
      se=sqrt(phi)
c                                            ---- 1st order differential
      if ( nreq.ge.1 ) then
        do i=1,6
          dseds(i)=v(i)/se
        enddo
      endif
c                                            ---- 2nd order differential
      if ( nreq.ge.2 ) then
        do i=1,6
          do j=1,6
            d2seds2(i,j)=( -v(i)*v(j)/phi+c(i,j) )/se
          enddo
        enddo
      endif
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     Vegter(2006) yield function and its differentials
c     (International Journal of Plasticity 22(2006) P-557-580)
c
c     This is plane stress yield function.
c
c     Version1.5 (Build20121224)
c        coded by Tomokage Inoue(Aisin AW Industries Co.,Ltd.)
c     Version1.5a 20150225
c        modified by H.Takizawa
c
      subroutine jancae_vegter ( s,se,dseds,d2seds2,nreq,
     &                           pryld,ndyld )
c
      implicit none
      integer, intent(in) :: nreq , ndyld
      real*8,  intent(in) :: s(3), pryld(ndyld)
      real*8,  intent(out) :: se, dseds(3), d2seds2(3,3)
c
      integer  :: nf
c
      nf=nint(pryld(2))-1
      call jancae_vegter_core ( s,se,dseds,d2seds2,nreq,
     &                          pryld,ndyld,nf )
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c
      subroutine jancae_vegter_core ( s,se,dseds,d2seds2,nreq,
     &                                pryld,ndyld,nf )
c-----------------------------------------------------------------------
c
      implicit none
c                                                    ---- input & output
      integer, intent(in) :: nreq , ndyld , nf
      real*8, intent(in) :: s(3), pryld(ndyld)
      real*8, intent(out) :: se, dseds(3), d2seds2(3,3)
c
c                                                   ---- local variables
      real*8, parameter :: pi=3.141592653589793d0
      real*8, parameter :: TOL0=1.0d-8
      real*8, parameter :: TOL2=1.0d-2
      real*8, parameter :: TOL =1.0d-4
c
c        TOL0    : exception treatment tolerance of stress state
c        TOL2    : se tolerance change f(1) to f(2)
c        TOL     : exception treatment tolerance of vsin2t
c
c ---------------------------------
      real*8 x(4),vsqrt,vcos2t,vsin2t,theta,theta_rv
      integer i, j, k, m, n,isflag,iareaflag
      real*8 phi_un(0:nf),phi_sh(0:nf),phi_ps(0:nf),f_bi0,
     &       omg(0:nf),r_bi0
      real*8 fun1,fun2,run,fsh1,fsh2,rsh,fps1,fps2,rps,fbi1,fbi2,rbi,
     &       fun1r,fun2r,runr,fps1r,fps2r,rpsr
      real*8 alfa,a(2),b(2),c(2),mm(2),nn(2),mu,f(2)
      real*8 beta,aa,bb,cc,dd
c
c                            over this line for equivalent stress
c                            ---------------------------------------
c
      real*8 dxds(3,3),dphidx(3),dfdmu(2),dadc(2),dcdc(2),dbdc(2)
      real*8 dndc(2),dmdc(2),dmdctmp,dndctmp,P(2),nnmm, dfdc(2)
c  -------------------------------------------
c
c                            over this line for 1st order differential
c                            -----------------------------------------
c
      real*8 d2adc2(2),d2bdc2(2),d2cdc2(2),d2ndc2(2),d2mdc2(2),
     &       vvtmp(0:6),vvtmp_rv(0:6),d2mdc2tmp,d2ndc2tmp,
     &       d2phidx2(3,3)
      integer ithetaflag
c
c                            over this line for 2nd order differential
c-----------------------------------------------------------------------
c                             ---- vegter's parameters from main routine
c     pryld(2) : nf           ! fourier term
c     pryld(3) : f_bi0        ! eq. bi-axial
c     pryld(4) : r_bi0        ! eq. bi-axial
c     i : 1~nf
c     n0=4+i*4
c     test angle  : 90/(nf)*i
c     pryld(n0+1) : phi_un(i) ! uniaxial
c     pryld(n0+2) : phi_sh(i) ! pure shear
c     pryld(n0+3) : phi_ps(i) ! plane strain
c     pryld(n0+4) : omg(   i)
c
c     this program assumes 2 type of tests
c       type 1 : 3 tests (0,45,90 deg)
c                set :phi_un,phi_sh,phi_ps,omg (1 to 3)
c                     phi_un,phi_sh,phi_ps,omg (4 to 6) =0.0
c       type 2 : 6 tests (0,15,30,45,60,75,90 deg)
c                set :phi_un,phi_sh,phi_ps,omg (1 to 6)
c
c
c
c     2015.2.25 new storage format : H.Takizawa
      f_bi0=     pryld(3)
      r_bi0=     pryld(4)
      do i=0,nf
        phi_un(i)=pryld(4+i*4+1)
        phi_sh(i)=pryld(4+i*4+2)
        phi_ps(i)=pryld(4+i*4+3)
        omg(   i)=pryld(4+i*4+4)
      enddo
c
c
      se=0.0
      do i=1,3
        se=se+s(i)**2
      enddo
      if ( se.le.0.0 ) then
        se=0.0
        return
      endif

c                                       ---- calc x(i) from eq.(14)~(17)
c     x(i)    : Principal stress(i=1^2)
c             : cos2theta(i=3), sin2theta(i=4)
c     isflag  : 0 s(i,i=1,3)=0
c             : 1 not s(i)=0
c                                 ---- exception treatment if all s(i)=0
c
      if(abs(s(1)).le.TOL0.and.abs(s(2)).le.TOL0.and.
     &                                         abs(s(3)).le.TOL0)then
      isflag=0
      call jancae_clear1 ( x,4 )
      goto 100
c
      else
      isflag=1
      call jancae_clear1 ( x,4 )
c
c                           ---- exception treatment if s(1)=s(2),s(3)=0
c
      if(abs(s(1)-s(2)).le.TOL0.and.abs(s(3)).le.TOL0) then
        theta=0.0d0
        theta_rv=0.5d0*pi
        x(1)=0.5d0*(s(1)+s(2))
        x(2)=0.5d0*(s(1)+s(2))
        x(3)=cos(2.0d0*theta)
        x(4)=sin(2.0d0*theta)
        vcos2t=x(3)
        vsin2t=x(4)
c
c                                                    ---- normal process
      else
        vsqrt=sqrt((s(1)-s(2))**2+4.0d0*s(3)**2)
        x(1)=0.5d0*(s(1)+s(2)+vsqrt)
        x(2)=0.5d0*(s(1)+s(2)-vsqrt)
        x(3)=(s(1)-s(2))/vsqrt
        x(4)=2.0d0*s(3)/vsqrt
        vcos2t=x(3)
        vsin2t=x(4)
        theta=0.5d0*acos(vcos2t)
        theta_rv=0.5d0*pi-theta
        endif
      endif
c
c                        ---- calc fk(k=un,sh,ps,bi) , rk(k=un,sh,ps,bi)
c
c                                                          ! un=uniaxial
        fun1=0.0d0
        fun1r=0.0d0
        run=0.0d0
        runr=0.0d0
      do m=0,nf
        fun1=fun1+phi_un(m)*cos(2.0d0*dble(m)*theta)
        fun1r=fun1r+phi_un(m)*cos(2.0d0*dble(m)*theta_rv)
        run=run+omg(m)*cos(2.0d0*dble(m)*theta)
        runr=runr+omg(m)*cos(2.0d0*dble(m)*theta_rv)
      enddo
        fun2=0.0d0
        fun2r=0.0d0
c                                                        ! sh=pure shear
        fsh1=0.0d0
        fsh2=0.0d0
      do m=0,nf
        fsh1=fsh1+phi_sh(m)*cos(2.0d0*dble(m)*theta)
        fsh2=fsh2-phi_sh(m)*cos(2.0d0*dble(m)*theta_rv)
      enddo
        rsh=-1.0d0
c                                                      ! ps=plane strain
        fps1=0.0d0
        fps1r=0.0d0
        rps=0.0d0
        rpsr=0.0d0
      do m=0,nf
        fps1=fps1+phi_ps(m)*cos(2.0d0*dble(m)*theta)
        fps2=0.5d0*fps1
        fps1r=fps1r+phi_ps(m)*cos(2.0d0*dble(m)*theta_rv)
        fps2r=0.5d0*fps1r
      enddo
        rps=-0.0d0
        rpsr=0.0d0
c                                                      ! bi=equi-biaxial
        fbi1=f_bi0
        fbi2=fbi1
        rbi=((r_bi0+1.0d0)+(r_bi0-1.0d0)*vcos2t)/
     &           ((r_bi0+1.0d0)-(r_bi0-1.0d0)*vcos2t)
c
c                                 ---- case distribution by stress state
      if(x(1).ne.0.0d0)then
        alfa=x(2)/x(1)
      endif
      if(x(2).ne.0.0d0)then
        beta=x(1)/x(2)
      endif
c
c     iareaflag    :stress state flag(i=0~6)
c
      if(x(1).gt.0.0d0.and.alfa.lt.0.0d0.and.alfa.ge.fsh2/fsh1) then
        iareaflag=1
      else if(x(1).gt.0.0d0.and.alfa.ge.0.0d0
     &                           .and.alfa.lt.fps2/fps1) then
        iareaflag=2
      else if(x(1).gt.0.0d0.and.alfa.ge.fps2/fps1
     &                           .and.alfa.le.1.0d0) then
        iareaflag=3
c
      else if(x(1).lt.0.0d0.and.alfa.ge.1.0d0
     &                           .and.alfa.lt.fps1r/fps2r) then
        iareaflag=4
      else if(x(1).lt.0.0d0.and.beta.le.fps2r/fps1r
     &                                      .and.beta.gt.0.0d0) then
        iareaflag=5
      else if(x(1).ge.0.0d0.and.beta.le.0.0d0
     &                           .and.beta.gt.fsh1/fsh2) then
        iareaflag=6
c
      else
        go to 100
      endif
c
c                                       ---- calc. hingepoint b(i,i=1~2)
      select case ( iareaflag )
c
      case ( 1 )                                           ! iareaflag=1
         a(1)=fsh1
         a(2)=fsh2
         c(1)=fun1
         c(2)=fun2
         nn(1)=1.0d0
         nn(2)=rsh
         mm(1)=1.0d0
         mm(2)=run
         call jancae_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
c
      case ( 2 )                                           ! iareaflag=2
         a(1)=fun1
         a(2)=fun2
         c(1)=fps1
         c(2)=fps2
         nn(1)=1.0d0
         nn(2)=run
         mm(1)=1.0d0
         mm(2)=rps
         call jancae_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
c
      case ( 3 )                                           ! iareaflag=3
         a(1)=fps1
         a(2)=fps2
         c(1)=fbi1
         c(2)=fbi2
         nn(1)=1.0d0
         nn(2)=rps
         mm(1)=1.0d0
         mm(2)=rbi
         call jancae_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
c
      case ( 4 )                                           ! iareaflag=4
         a(1)=-fbi1
         a(2)=-fbi2
         c(1)=-fps2r
         c(2)=-fps1r
         nn(1)=-1.0d0
         nn(2)=-rbi
         mm(1)=-rpsr
         mm(2)=-1.0d0
         call jancae_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
c
      case ( 5 )                                           ! iareaflag=5
         a(1)=-fps2r
         a(2)=-fps1r
         c(1)=-fun2r
         c(2)=-fun1r
         nn(1)=-rpsr
         nn(2)=-1.0d0
         mm(1)=-runr
         mm(2)=-1.0d0
         call jancae_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
c
      case ( 6 )                                           ! iareaflag=6
         a(1)=-fun2r
         a(2)=-fun1r
         c(1)=fsh1
         c(2)=fsh2
         nn(1)=-runr
         nn(2)=-1.0d0
         mm(1)=1.0d0
         mm(2)=rsh
         call jancae_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
c
      case default
         write (6,*) 'iareaflag error :',iareaflag
         call jancae_exit (9000)
      end select

c                            ---- calc. fourier coefficient mu(0<=mu<=1)
      call jancae_vegter_calc_mu ( x,a,b,c,mu,iareaflag,s,theta
     &                                            ,aa,bb,cc,dd )
c
c                            ---- calc. normalized yield locus f(i)i=1~2
      call jancae_vegter_calc_fi ( x,a,b,c,mu,f )
      go to 200
c
c                                                 ---- equivalent stress
  100 continue
      se=0.0d0
      go to 300
  200 continue
      if(f(1).le.TOL2) then
         se=x(2)/f(2)
      else
         se=x(1)/f(1)
      endif
c
      go to 300
c
  300 continue
c
c                                            ---- 1st order differential
c
      if ( nreq.ge.1 ) then
c                        ---- set dadc,dcdc,dndc,dmdc for eq.(A.7)^(A.9)
c
      call jancae_clear1 ( dadc,2 )
      call jancae_clear1 ( dbdc,2 )
      call jancae_clear1 ( dcdc,2 )
      call jancae_clear1 ( dndc,2 )
      call jancae_clear1 ( dmdc,2 )
c
      select case ( iareaflag )
c
      case ( 1 )                                           ! iareaflag=1
       if(abs(vsin2t).ge.TOL) then
        do m=0,nf
          dadc(1)=dadc(1)+phi_sh(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta)
     &                             /sin(2.0d0*theta)
          dadc(2)=dadc(2)+phi_sh(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta_rv)
     &                             /sin(2.0d0*theta)
        enddo
       else
          do m=0,nf
              dadc(1)=dadc(1)+phi_sh(m)*dble(m)**2
              dadc(2)=dadc(2)+phi_sh(m)*dble(m)**2
          enddo
       endif
c
       if(abs(vsin2t).ge.TOL) then
         do m=0,nf
          dcdc(1)=dcdc(1)+phi_un(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta)
     &                             /sin(2.0d0*theta)
         enddo
        else
          do m=0,nf
              dcdc(1)=dcdc(1)+phi_un(m)*dble(m)**2
          enddo
        endif
          dcdc(2)=0.0d0
c
          dndc(1)=0.0d0
          dndc(2)=0.0d0
c
          dmdc(1)=0.0d0
       if(abs(vsin2t).ge.TOL) then
         do m=0,nf
          dmdc(2)=dmdc(2)+omg(m)*dble(m)
     &                          *sin(2.0d0*dble(m)*theta)
     &                          /sin(2.0d0*theta)
         enddo
        else
          do m=0,nf
              dmdc(2)=dmdc(2)+omg(m)*dble(m)**2
          enddo
        endif
      call jancae_vegter_calc_dbdc ( a,b,c,dadc,dbdc,dcdc,mm,nn,
     &                               dndc,dmdc,iareaflag,P,nnmm)
c
c
      case ( 2 )                                           ! iareaflag=2
       if(abs(vsin2t).ge.TOL) then
         do m=0,nf
          dadc(1)=dadc(1)+phi_un(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta)
     &                             /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dadc(1)=dadc(1)+phi_un(m)*dble(m)**2
         enddo
       endif
          dadc(2)=0.0d0
c
       if(abs(vsin2t).ge.TOL) then
         do m=0,nf
          dcdc(1)=dcdc(1)+phi_ps(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta)
     &                             /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dcdc(1)=dcdc(1)+phi_ps(m)*dble(m)**2
         enddo
       endif
          dcdc(2)=0.5d0*dcdc(1)
c
          dndc(1)=0.0d0
       if(abs(vsin2t).ge.TOL) then
        do m=0,nf
          dndc(2)=dndc(2)+omg(m)*dble(m)
     &                           *sin(2.0d0*dble(m)*theta)
     &                           /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dndc(2)=dndc(2)+omg(m)*dble(m)**2
         enddo
       endif
c
          dmdc(1)=0.0d0
          dmdc(2)=0.0d0
      call jancae_vegter_calc_dbdc ( a,b,c,dadc,dbdc,dcdc,mm,nn,
     &                               dndc,dmdc,iareaflag,P,nnmm)
c
c
      case ( 3 )                                           ! iareaflag=3
       if(abs(vsin2t).ge.TOL) then
        do m=0,nf
          dadc(1)=dadc(1)+phi_ps(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta)
     &                             /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dadc(1)=dadc(1)+phi_ps(m)*dble(m)**2
         enddo
       endif
          dadc(2)=0.5d0*dadc(1)
c
          dcdc(1)=0.0d0
          dcdc(2)=0.0d0
c
          dndc(1)=0.0d0
          dndc(2)=0.0d0
c
          dmdc(1)=0.0d0
          dmdctmp=r_bi0+1.0d0-(r_bi0-1.0d0)*vcos2t
          dmdc(2)=2.0d0*(r_bi0*r_bi0-1.0d0)/(dmdctmp*dmdctmp)
      call jancae_vegter_calc_dbdc ( a,b,c,dadc,dbdc,dcdc,mm,nn,
     &                               dndc,dmdc,iareaflag,P,nnmm)
c
c
      case ( 4 )                                           ! iareaflag=4
          dadc(1)=0.0d0
          dadc(2)=0.0d0
c
       if(abs(vsin2t).ge.TOL) then
        do m=0,nf
          dcdc(2)=dcdc(2)+phi_ps(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta_rv)
     &                             /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dcdc(2)=dcdc(2)+phi_ps(m)*dble(m)**2
         enddo
       endif
          dcdc(1)=0.5d0*dcdc(2)
c
          dndc(1)=0.0d0
          dndctmp=r_bi0+1.0d0-(r_bi0-1.0d0)*vcos2t
          dndc(2)=-2.0d0*(r_bi0*r_bi0-1.0d0)/(dndctmp*dndctmp)
c
          dmdc(1)=0.0d0
          dmdc(2)=0.0d0
      call jancae_vegter_calc_dbdc ( a,b,c,dadc,dbdc,dcdc,mm,nn,
     &                               dndc,dmdc,iareaflag,P,nnmm)
c
c
      case ( 5 )                                           ! iareaflag=5
       if(abs(vsin2t).ge.TOL) then
        do m=0,nf
          dadc(2)=dadc(2)+phi_ps(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta_rv)
     &                             /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dadc(2)=dadc(2)+phi_ps(m)*dble(m)**2
         enddo
       endif
          dadc(1)=0.5d0*dadc(2)
c
          dcdc(1)=0.0d0
       if(abs(vsin2t).ge.TOL) then
        do m=0,nf
          dcdc(2)=dcdc(2)+phi_un(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta_rv)
     &                             /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dcdc(2)=dcdc(2)+phi_un(m)*dble(m)**2
         enddo
       endif
c
          dndc(1)=0.0d0
          dndc(2)=0.0d0
c
       if(abs(vsin2t).ge.TOL) then
        do m=0,nf
          dmdc(1)=dmdc(1)+omg(m)*dble(m)
     &                          *sin(2.0d0*dble(m)*theta_rv)
     &                          /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dmdc(1)=dmdc(1)+omg(m)*dble(m)**2
         enddo
       endif
          dmdc(2)=0.0d0
      call jancae_vegter_calc_dbdc ( a,b,c,dadc,dbdc,dcdc,mm,nn,
     &                               dndc,dmdc,iareaflag,P,nnmm)
c
c
      case ( 6 )                                           ! iareaflag=6
          dadc(1)=0.0d0
       if(abs(vsin2t).ge.TOL) then
        do m=0,nf
          dadc(2)=dadc(2)+phi_un(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta_rv)
     &                             /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dadc(2)=dadc(2)+phi_un(m)*dble(m)**2
         enddo
       endif
c
       if(abs(vsin2t).ge.TOL) then
        do m=0,nf
          dcdc(1)=dcdc(1)+phi_sh(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta)
     &                             /sin(2.0d0*theta)
          dcdc(2)=dcdc(2)+phi_sh(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta_rv)
     &                             /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dcdc(1)=dcdc(1)+phi_sh(m)*dble(m)**2
          dcdc(2)=dcdc(2)+phi_sh(m)*dble(m)**2
         enddo
       endif
c
       if(abs(vsin2t).ge.TOL) then
        do m=0,nf
          dndc(1)=dndc(1)+omg(m)*dble(m)
     &                          *sin(2.0d0*dble(m)*theta_rv)
     &                          /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dndc(1)=dndc(1)+omg(m)*dble(m)**2
         enddo
       endif
          dndc(2)=0.0d0
c
          dmdc(1)=0.0d0
          dmdc(2)=0.0d0
      call jancae_vegter_calc_dbdc ( a,b,c,dadc,dbdc,dcdc,mm,nn,
     &                               dndc,dmdc,iareaflag,P,nnmm)
c
      case default
        write (6,*) 'iareaflag error(dseds) :',iareaflag
        call jancae_exit (9000)
      end select
c
c
c                             ---- calc. dphidx(i) (i=1~3)  eq.(21)~(23)
      call jancae_vegter_calc_dphidx ( dphidx,se,a,b,c,dadc,dbdc,dcdc,
     &                      mm,nn,dndc,dmdc,f,iareaflag,mu,x,dfdmu,dfdc)
c
c
c                                   ---- calc. dseds(i) (i=1~3)  eq.(20)
      call jancae_vegter_calc_dseds ( dseds,x,dphidx,vcos2t,vsin2t,
     &                                iareaflag,isflag,dxds)
c
      endif
c
c
c                                            ---- 2nd order differential
      if ( nreq.ge.2 ) then
c                     ---- set d2adc2,d2cdc2,d2ndc2,d2mdc2 for d2bdc2(2)
c
      call jancae_clear1 ( d2adc2,2 )
      call jancae_clear1 ( d2bdc2,2 )
      call jancae_clear1 ( d2cdc2,2 )
      call jancae_clear1 ( d2ndc2,2 )
      call jancae_clear1 ( d2mdc2,2 )
c
c                     ---- define exception treatment condition of theta
c                        ---- if theta<=0.002865deg then apply exception
c
      if(abs(vsin2t).le.TOL) then
         ithetaflag=1
      else
         ithetaflag=0
      endif
c
      if(ithetaflag.eq.1) then
        vvtmp(0)=0.0d0
        vvtmp(1)=0.0d0
        vvtmp(2)=2.0d0
        vvtmp(3)=8.0d0
        vvtmp(4)=20.0d0
        vvtmp(5)=40.0d0
        vvtmp(6)=70.0d0
c
        vvtmp_rv(0)=0.0d0
        vvtmp_rv(1)=0.0d0
        vvtmp_rv(2)=-2.0d0
        vvtmp_rv(3)=8.0d0
        vvtmp_rv(4)=-20.0d0
        vvtmp_rv(5)=40.0d0
        vvtmp_rv(6)=-70.0d0
c
      else
       do m=0,nf
        vvtmp(m)=cos(2.0d0*theta)*sin(2.0d0*dble(m)*theta)/
     &          (sin(2.0d0*theta)**3)-dble(m)*cos(2.0d0*dble(m)*theta)/
     &                                           (sin(2.0d0*theta)**2)
        vvtmp_rv(m)=cos(2.0d0*theta)*sin(2.0d0*dble(m)*theta_rv)/
     &       (sin(2.0d0*theta)**3)+dble(m)*cos(2.0d0*dble(m)*theta_rv)/
     &                                           (sin(2.0d0*theta)**2)
       enddo
c
      endif
c
      select case ( iareaflag )
c
      case ( 1 )                                           ! iareaflag=1
         do m=0,nf
          d2adc2(1)=d2adc2(1)+phi_sh(m)*dble(m)*vvtmp(m)
          d2adc2(2)=d2adc2(2)+phi_sh(m)*dble(m)*vvtmp_rv(m)
         enddo
c
         do m=0,nf
          d2cdc2(1)=d2cdc2(1)+phi_un(m)*dble(m)*vvtmp(m)
         enddo
          d2cdc2(2)=0.0d0
c
          d2ndc2(1)=0.0d0
          d2ndc2(2)=0.0d0
c
          d2mdc2(1)=0.0d0
         do m=0,nf
          d2mdc2(2)=d2mdc2(2)+omg(m)*dble(m)*vvtmp(m)
         enddo
      call jancae_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,nn,dndc,
     &     dmdc,iareaflag,d2adc2,d2bdc2,d2cdc2,d2ndc2,d2mdc2,P,nnmm,s)
c
c
      case ( 2 )                                           ! iareaflag=2
         do m=0,nf
          d2adc2(1)=d2adc2(1)+phi_un(m)*dble(m)*vvtmp(m)
         enddo
          d2adc2(2)=0.0d0
c
         do m=0,nf
          d2cdc2(1)=d2cdc2(1)+phi_ps(m)*dble(m)*vvtmp(m)
         enddo
          d2cdc2(2)=0.5d0*d2cdc2(1)
c
          d2ndc2(1)=0.0d0
         do m=0,nf
          d2ndc2(2)=d2ndc2(2)+omg(m)*dble(m)*vvtmp(m)
         enddo
c
          d2mdc2(1)=0.0d0
          d2mdc2(2)=0.0d0
      call jancae_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,nn,dndc,
     &     dmdc,iareaflag,d2adc2,d2bdc2,d2cdc2,d2ndc2,d2mdc2,P,nnmm,s)
c
c
      case ( 3 )                                           ! iareaflag=3
         do m=0,nf
          d2adc2(1)=d2adc2(1)+phi_ps(m)*dble(m)*vvtmp(m)
         enddo
          d2adc2(2)=0.5d0*d2adc2(1)
c
          d2cdc2(1)=0.0d0
          d2cdc2(2)=0.0d0
c
          d2ndc2(1)=0.0d0
          d2ndc2(2)=0.0d0
c
          d2mdc2(1)=0.0d0
             d2mdc2tmp=r_bi0+1.0d0-(r_bi0-1.0d0)*vcos2t
          d2mdc2(2)=4.0d0*(r_bi0**2-1.0d0)*(r_bi0-1.0d0)/
     &                                              (d2mdc2tmp**3)
      call jancae_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,nn,dndc,
     &     dmdc,iareaflag,d2adc2,d2bdc2,d2cdc2,d2ndc2,d2mdc2,P,nnmm,s)
c
c
      case ( 4 )                                          !  iareaflag=4
          d2adc2(1)=0.0d0
          d2adc2(2)=0.0d0
c
         do m=0,nf
          d2cdc2(2)=d2cdc2(2)+phi_ps(m)*dble(m)*vvtmp_rv(m)
         enddo
          d2cdc2(1)=0.5d0*d2cdc2(2)
c
          d2ndc2(1)=0.0d0
             d2ndc2tmp=r_bi0+1.0d0-(r_bi0-1.0d0)*vcos2t
          d2ndc2(2)=-4.0d0*(r_bi0**2-1.0d0)*(r_bi0-1.0d0)/
     &                                              (d2ndc2tmp**3)
c
          d2mdc2(1)=0.0d0
          d2mdc2(2)=0.0d0
      call jancae_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,nn,dndc,
     &     dmdc,iareaflag,d2adc2,d2bdc2,d2cdc2,d2ndc2,d2mdc2,P,nnmm,s)
c
c
      case ( 5 )                                           ! iareaflag=5
         do m=0,nf
          d2adc2(2)=d2adc2(2)+phi_ps(m)*dble(m)*vvtmp_rv(m)
         enddo
          d2adc2(1)=0.5d0*d2adc2(2)
c
          d2cdc2(1)=0.0d0
         do m=0,nf
          d2cdc2(2)=d2cdc2(2)+phi_un(m)*dble(m)*vvtmp_rv(m)
         enddo
c
          d2ndc2(1)=0.0d0
          d2ndc2(2)=0.0d0
c
         do m=0,nf
          d2mdc2(1)=d2mdc2(1)+omg(m)*dble(m)*vvtmp_rv(m)
         enddo
          d2mdc2(2)=0.0d0
      call jancae_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,nn,dndc,
     &     dmdc,iareaflag,d2adc2,d2bdc2,d2cdc2,d2ndc2,d2mdc2,P,nnmm,s)
c
c
      case ( 6 )                                           ! iareaflag=6
          d2adc2(1)=0.0d0
         do m=0,nf
          d2adc2(2)=d2adc2(2)+phi_un(m)*dble(m)*vvtmp_rv(m)
         enddo
c
         do m=0,nf
          d2cdc2(1)=d2cdc2(1)+phi_sh(m)*dble(m)*vvtmp(m)
          d2cdc2(2)=d2cdc2(2)+phi_sh(m)*dble(m)*vvtmp_rv(m)
         enddo
cn
         do m=0,nf
          d2ndc2(1)=d2ndc2(1)+omg(m)*dble(m)*vvtmp_rv(m)
         enddo
          d2ndc2(2)=0.0d0
c
          d2mdc2(1)=0.0d0
          d2mdc2(2)=0.0d0
      call jancae_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,nn,dndc,
     &     dmdc,iareaflag,d2adc2,d2bdc2,d2cdc2,d2ndc2,d2mdc2,P,nnmm,s)
c
      case default
        write (6,*) 'iareaflag error(d2seds2) :',iareaflag
        call jancae_exit (9000)
      end select
c
c
c                                     ---- calc. d2phidx2(k,l) (k,l=1~3)
      call jancae_vegter_calc_d2phidx2 (d2phidx2,se,a,b,c,dadc,dbdc,
     &             dcdc,mm,nn,dndc,dmdc,f,iareaflag,mu,x,d2adc2,d2bdc2,
     &            d2cdc2,d2ndc2,d2mdc2,dfdmu,dfdc,s,aa,bb,cc,dd,dphidx)
c
c
c                                      ---- calc. d2seds2(i,j) (i,j=1~3)
      call jancae_vegter_calc_d2seds2 (d2seds2,d2phidx2,se,a,b,c,mu,x,
     &         vcos2t,vsin2t,iareaflag,dxds,dphidx,isflag,s,dseds,
     &                                 pryld,ndyld)
c
      endif
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c                                      under this line is branch_routine
c---------------------------------------------------------------(vegter)
c     calc. hingepoint b(i,i=1~2)
c
      subroutine jancae_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
c-----------------------------------------------------------------------
      implicit none
c
      real*8, intent(in) :: a(2),c(2),mm(2),nn(2),s(3)
      real*8, intent(out) :: b(2)
      real*8  bb, b1u,b2u
      real*8, parameter :: TOL=1.0d-8
      integer, intent(in) ::iareaflag
c
      b1u=mm(2)*(nn(1)*a(1)+nn(2)*a(2))-nn(2)*(mm(1)*c(1)+mm(2)*c(2))
      b2u=nn(1)*(mm(1)*c(1)+mm(2)*c(2))-mm(1)*(nn(1)*a(1)+nn(2)*a(2))
c
      bb=nn(1)*mm(2)-mm(1)*nn(2)
      if(abs(bb).le.TOL) then
         write (6,*) 'hingepoint singular error! '
         call jancae_exit (9000)
      endif
c
      b(1)=b1u/bb
      b(2)=b2u/bb
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     calc. fourier coefficient mu(0<=mu<=1)
c
      subroutine jancae_vegter_calc_mu ( x,a,b,c,mu,iareaflag,s,theta,
     &                                                     aa,bb,cc,dd)
c-----------------------------------------------------------------------
      implicit none
c
      real*8, intent(in) :: x(4),a(2),b(2),c(2),s(3),theta
      real*8, intent(out) :: mu
      real*8  aa,bb,cc,dd,xx(2)
      real*8, parameter :: TOL1=1.0d-8
      real*8, parameter :: TOL2=1.0d-8
      integer, intent(in) ::iareaflag
      integer imuflag
c
      aa=x(2)*(a(1)+c(1)-2.0d0*b(1))-x(1)*(a(2)+c(2)-2.0d0*b(2))
      bb=2.0d0*x(2)*(b(1)-a(1))-2.0d0*x(1)*(b(2)-a(2))
      cc=x(2)*a(1)-x(1)*a(2)
c
      if(abs(aa).le.TOL1) then
         write (6,*) 'calc. mu singular error! ',abs(aa),iareaflag
         call jancae_exit (9000)
      endif
c
      dd=bb*bb-4.0d0*aa*cc
      if(dd.ge.0.0d0) then
        xx(1)=0.5d0*(-bb+sign(sqrt(dd),-bb))/aa
        xx(2)=cc/(aa*xx(1))
c
      else
         write (6,*) 'negative dd ! ',dd,iareaflag
         call jancae_exit (9000)
      endif
c
      if(xx(1).ge.0.0d0.and.xx(1).le.1.0000005d0) then
          mu=xx(1)
           imuflag=1
        else if(xx(2).ge.0.0d0.and.xx(2).le.1.0000005d0) then
          mu=xx(2)
           imuflag=2
        else if(abs(xx(1)).le.TOL2.or.abs(xx(2)).le.TOL2)then
          mu=0.0d0
        else
         write (6,*) 'can not find mu ! solve error ',iareaflag,xx(1)
     & ,xx(2)
         call jancae_exit (9000)
      endif
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     calc. normalized yield locus f(i)i=1~2
c
      subroutine jancae_vegter_calc_fi ( x,a,b,c,mu,f )
c-----------------------------------------------------------------------
      implicit none
c
      real*8, intent(in) :: x(4),a(2),b(2),c(2),mu
      real*8, intent(out) :: f(2)
      integer i
c
      do i=1,2
        f(i)=a(i)+2.0d0*mu*(b(i)-a(i))+mu*mu*(a(i)+c(i)-2.0d0*b(i))
      enddo
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     calc. dbdc(i) (i=1~2) eq.(A.7)
c
      subroutine jancae_vegter_calc_dbdc ( a,b,c,dadc,dbdc,
     &                 dcdc, mm,nn,dndc,dmdc,iareaflag,P,nnmm)
c-----------------------------------------------------------------------
      implicit none
c
      real*8, intent(in) :: a(2),b(2),c(2),dadc(2),dcdc(2),mm(2),
     &                      nn(2),dndc(2),dmdc(2)
      real*8, intent(out) :: dbdc(2),P(2),nnmm
      integer, intent(in) :: iareaflag
      real*8, parameter :: TOL=1.0d-8
      real*8 nminv
      integer i
c
      P(1)=nn(1)*dadc(1)+dndc(1)*(a(1)-b(1))+nn(2)*dadc(2)+
     &                                        dndc(2)*(a(2)-b(2))
c
      P(2)=mm(1)*dcdc(1)+dmdc(1)*(c(1)-b(1))+mm(2)*dcdc(2)+
     &                                        dmdc(2)*(c(2)-b(2))
c
      nnmm=nn(1)*mm(2)-mm(1)*nn(2)
         if(abs(nnmm).lt.TOL) then
            write (6,*) 'nnmm too small! ',nnmm
            call jancae_exit (9000)
         endif
      nminv=1.0d0/nnmm
c
      dbdc(1)=nminv*(P(1)*mm(2)-P(2)*nn(2))
      dbdc(2)=nminv*(P(2)*nn(1)-P(1)*mm(1))
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     calc. dphidx(i) (i=1~3)  eq.(21)~(23)
c
      subroutine jancae_vegter_calc_dphidx ( dphidx,se,a,b,c,dadc,
     &            dbdc,dcdc,mm,nn,dndc,dmdc,f,iareaflag,mu,x,dfdmu,dfdc)
c-----------------------------------------------------------------------
      implicit none
c
      real*8, intent(in) :: se,a(2),b(2),c(2),dadc(2),dbdc(2),dcdc(2),
     &                      mm(2),nn(2),dndc(2),dmdc(2),f(2),mu,x(4)
      real*8, intent(out) :: dphidx(3),dfdmu(2),dfdc(2)
      integer, intent(in) :: iareaflag
      real*8, parameter :: TOL1=0.996515d0       ! =44.9deg
      real*8, parameter :: TOL2=1.003497d0       ! =45.1deg
      real*8, parameter :: TOL3=1.0d-8
      real*8 dphidxtmp(3),dphidxcoe,dphidxcinv,
     &       tmp_u,tmp_b,vtan
      integer i
c
c                                    ---- calc. dfdc(i) (i=1~2)  eq.(23)
      do i=1,2
        dfdc(i)=dadc(i)+2.0d0*mu*(dbdc(i)-dadc(i))+mu*mu*(dadc(i)
     &                                       +dcdc(i)-2.0d0*dbdc(i))
      enddo
c
c                                   ---- calc. dfdmu(i) (i=1~2)  eq.(22)
      do i=1,2
        dfdmu(i)=2.0d0*(b(i)-a(i))+2.0d0*mu*(a(i)+c(i)-2.0d0*b(i))
      enddo
c
c                            ---- calc. dphidx(i) (i=1~3)  eq.(21),(C.1)
      dphidxcinv=f(1)*dfdmu(2)-f(2)*dfdmu(1)
         if(abs(dphidxcinv).lt.TOL3) then
            write (6,*) 'eq.(21) too small! ',dphidxcinv
            call jancae_exit (9000)
         endif
             dphidxcoe=1.0d0/dphidxcinv
c
c                            ---- if condition to avoid singular eq.(20)
c                                            ---- apply 44.9 to 45.1 deg
c
      if(iareaflag.eq.3.or.iareaflag.eq.4) then
      vtan=x(2)/x(1)
      end if
c
      if(iareaflag.eq.4.and.vtan.ge.TOL1.and.vtan.le.TOL2) then
c
      tmp_u=1.0d0*(2.0d0*(1.0d0-mu)*dbdc(2)+mu*dcdc(2))*dfdmu(1)
     &     -1.0d0*(2.0d0*(1.0d0-mu)*dbdc(1)+mu*dcdc(1))*dfdmu(2)
c
      tmp_b=2.0d0*(1.0d0-mu)*(b(1)-b(2))+mu*(c(1)-c(2))

      dphidxtmp(1)=dfdmu(2)
      dphidxtmp(2)=-dfdmu(1)
      dphidxtmp(3)=tmp_u/tmp_b
c
      else if(iareaflag.eq.3.and.vtan.ge.TOL1.and.vtan.le.TOL2) then
c
      tmp_u=1.0d0*(2.0d0*mu*dbdc(2)+(1.0d0-mu)*dadc(2))*dfdmu(1)
     &     -1.0d0*(2.0d0*mu*dbdc(1)+(1.0d0-mu)*dadc(1))*dfdmu(2)
c
      tmp_b=2.0d0*mu*(b(1)-b(2))+(1.0d0-mu)*(a(1)-a(2))

      dphidxtmp(1)=dfdmu(2)
      dphidxtmp(2)=-dfdmu(1)
      dphidxtmp(3)=tmp_u/tmp_b
c
      else
c
      dphidxtmp(1)=dfdmu(2)
      dphidxtmp(2)=-dfdmu(1)
      dphidxtmp(3)=se*(dfdc(2)*dfdmu(1)-dfdc(1)*dfdmu(2))
c
      endif
      do i=1,3
        dphidx(i)=dphidxcoe*dphidxtmp(i)
      enddo
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     calc. dseds(i) (i=1~3)  eq.(20)
c
      subroutine jancae_vegter_calc_dseds ( dseds,x,dphidx,vcos2t,
     &                                     vsin2t,iareaflag,isflag,dxds)
c-----------------------------------------------------------------------
      implicit none
c
      real*8, intent(in) :: x(4),dphidx(3),vcos2t,vsin2t
      real*8, intent(out) :: dseds(3),dxds(3,3)
      integer, intent(in) :: iareaflag,isflag
      real*8, parameter :: TOL1=0.996515d0       ! =44.9deg
      real*8, parameter :: TOL2=1.003497d0       ! =45.1deg
      real*8 dxds_t(3,3),vtan
      integer i,j
c
c                  ---- set linear transformation matrix dxds(3,3) eq.18
c
c                            ---- if condition to avoid singular eq.(20)
c                                            ---- apply 44.9 to 45.1 deg
c
      if(iareaflag.eq.3.or.iareaflag.eq.4) then
      vtan=x(2)/x(1)
      end if
c
      if(iareaflag.eq.3.and.vtan.ge.TOL1.and.vtan.le.TOL2) then
        dxds(1,1)=0.5d0*(1.0d0+vcos2t)
        dxds(2,1)=0.5d0*(1.0d0-vcos2t)
        dxds(3,1)=vsin2t*vsin2t
        dxds(1,2)=0.5d0*(1.0d0-vcos2t)
        dxds(2,2)=0.5d0*(1.0d0+vcos2t)
        dxds(3,2)=-vsin2t*vsin2t
        dxds(1,3)=vsin2t
        dxds(2,3)=-vsin2t
        dxds(3,3)=-2.0d0*vsin2t*vcos2t
        dxds_t=transpose(dxds)
c
      else if(iareaflag.eq.4.and.vtan.ge.TOL1.and.vtan.le.TOL2) then
        dxds(1,1)=0.5d0*(1.0d0+vcos2t)
        dxds(2,1)=0.5d0*(1.0d0-vcos2t)
        dxds(3,1)=vsin2t*vsin2t
        dxds(1,2)=0.5d0*(1.0d0-vcos2t)
        dxds(2,2)=0.5d0*(1.0d0+vcos2t)
        dxds(3,2)=-vsin2t*vsin2t
        dxds(1,3)=vsin2t
        dxds(2,3)=-vsin2t
        dxds(3,3)=-2.0d0*vsin2t*vcos2t
        dxds_t=transpose(dxds)
c
      else
        dxds(1,1)=0.5d0*(1.0d0+vcos2t)
        dxds(2,1)=0.5d0*(1.0d0-vcos2t)
        dxds(3,1)=vsin2t*vsin2t/(x(1)-x(2))
        dxds(1,2)=0.5d0*(1.0d0-vcos2t)
        dxds(2,2)=0.5d0*(1.0d0+vcos2t)
        dxds(3,2)=-vsin2t*vsin2t/(x(1)-x(2))
        dxds(1,3)=vsin2t
        dxds(2,3)=-vsin2t
        dxds(3,3)=-2.0d0*vsin2t*vcos2t/(x(1)-x(2))
        dxds_t=transpose(dxds)
      endif
c
c                                   ---- calc. dseds(i) (1=1~3)  eq.(20)
      call jancae_clear1 ( dseds,3 )
c
      do i=1,3
        do j=1,3
          dseds(i)=dseds(i)+dxds_t(i,j)*dphidx(j)
        enddo
      enddo
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     calc. d2bdc2(i) (i=1~2)
c
      subroutine jancae_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,
     $      nn,dndc,dmdc,iareaflag,d2adc2,d2bdc2,d2cdc2,d2ndc2,d2mdc2,
     &      P,nnmm,s)
c-----------------------------------------------------------------------
      implicit none
c
      real*8, intent(in) :: a(2),b(2),c(2),dadc(2),dcdc(2),mm(2),
     &                      nn(2),dndc(2),dmdc(2),dbdc(2),
     &                      d2adc2(2),d2cdc2(2),d2ndc2(2),d2mdc2(2),
     &                      P(2),nnmm,s(3)
      real*8, intent(out) :: d2bdc2(2)
      integer, intent(in) :: iareaflag
      real*8  dnnmmdc,dPdc(2),dp1m2p2n2,dp2n1p1m1
      integer i
c
      dnnmmdc=dndc(1)*mm(2)+nn(1)*dmdc(2)-dmdc(1)*nn(2)-mm(1)*dndc(2)
c
      dPdc(1)=dndc(1)*dadc(1)+nn(1)*d2adc2(1)+d2ndc2(1)*(a(1)-b(1))
     &                                   +dndc(1)*(dadc(1)-dbdc(1))
     &       +dndc(2)*dadc(2)+nn(2)*d2adc2(2)+d2ndc2(2)*(a(2)-b(2))
     &                                   +dndc(2)*(dadc(2)-dbdc(2))
c
      dPdc(2)=dmdc(1)*dcdc(1)+mm(1)*d2cdc2(1)+d2mdc2(1)*(c(1)-b(1))
     &                                   +dmdc(1)*(dcdc(1)-dbdc(1))
     &       +dmdc(2)*dcdc(2)+mm(2)*d2cdc2(2)+d2mdc2(2)*(c(2)-b(2))
     &                                   +dmdc(2)*(dcdc(2)-dbdc(2))
c
      dp1m2p2n2=dPdc(1)*mm(2)+P(1)*dmdc(2)-dPdc(2)*nn(2)-P(2)*dndc(2)
      dp2n1p1m1=dPdc(2)*nn(1)+P(2)*dndc(1)-dPdc(1)*mm(1)-P(1)*dmdc(1)
c
      d2bdc2(1)=-1.0d0*dnnmmdc*(P(1)*mm(2)-P(2)*nn(2))/(nnmm*nnmm)
     &         +dp1m2p2n2/nnmm
      d2bdc2(2)=-1.0d0*dnnmmdc*(P(2)*nn(1)-P(1)*mm(1))/(nnmm*nnmm)
     &         +dp2n1p1m1/nnmm
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     calc. d2phidx2(k,l) (k,l=1~3)
c
      subroutine jancae_vegter_calc_d2phidx2 (d2phidx2,se,a,b,c,dadc,
     &              dbdc,dcdc,mm,nn,dndc,dmdc,f,iareaflag,mu,x,d2adc2,
     &      d2bdc2,d2cdc2,d2ndc2,d2mdc2,dfdmu,dfdc,s,aa,bb,cc,dd,dphidx)
c-----------------------------------------------------------------------
      implicit none
c
      real*8, intent(in) :: se,a(2),b(2),c(2),dadc(2),dbdc(2),dcdc(2),
     &                      mm(2),nn(2),dndc(2),dmdc(2),f(2),mu,x(4),
     &                      d2adc2(2),d2cdc2(2),d2ndc2(2),d2mdc2(2),
     &                              d2bdc2(2),dfdmu(2),dfdc(2),s(3),
     &                              aa,bb,cc,dd,dphidx(3)
      real*8, intent(out) :: d2phidx2(3,3)
      integer, intent(in) :: iareaflag
c     real*8, parameter :: TOL1=0.996515d0       ! =44.9deg
c     real*8, parameter :: TOL2=1.003497d0       ! =45.1deg
c
      integer i,j
      real*8 daadx(3),dbbdx(3),dccdx(3),ddddx(3),dmudx(3),
     &       d2fdmu2(2),d2fdmudc(2),d2fdcdmu(2),d2fdc2(2)
      real*8 vcommon,vtmp1,vtmp2,vtmp3,vtmp4
      real*8 va,vc,dvadx(3),dsedx(3),dvcdx(3)
c
c                                           ---- calc.  dmudx(i) (i=1~3)
      daadx(1)=-a(2)-c(2)+2.0d0*b(2)
      dbbdx(1)=-2.0d0*(b(2)-a(2))
      dccdx(1)=-a(2)
c
      daadx(2)=a(1)+c(1)-2.0d0*b(1)
      dbbdx(2)=2.0d0*(b(1)-a(1))
      dccdx(2)=a(1)
c
      daadx(3)=x(2)*(dadc(1)+dcdc(1)-2.0d0*dbdc(1))
     &        -x(1)*(dadc(2)+dcdc(2)-2.0d0*dbdc(2))
      dbbdx(3)=2.0d0*x(2)*
     &        (dbdc(1)-dadc(1))-2.0d0*x(1)*(dbdc(2)-dadc(2))
      dccdx(3)=x(2)*dadc(1)-x(1)*dadc(2)
c
      do i=1,3
       dmudx(i)=0.5d0*daadx(i)*(bb+sqrt(dd))/(aa*aa)
     &         +0.5d0*(-dbbdx(i)-0.5d0/(sqrt(dd))*(2.0d0*bb*dbbdx(i)
     &           -4.0d0*daadx(i)*cc-4.0d0*aa*dccdx(i)))/aa
      enddo
c
c                                         ---- calc.  d2fdmu2(i) (i=1~2)
      do i=1,2
         d2fdmu2(i)=2.0d0*(a(i)+c(i)-2.0d0*b(i))
      enddo
c
c                                        ---- calc.  d2fdmudc(i) (i=1~2)
      do i=1,2
        d2fdmudc(i)=2.0d0*(dbdc(i)-dadc(i))+2.0d0*mu*(dadc(i)+dcdc(i)
     &                                             -2.0d0*dbdc(i))
      enddo
c
c                                        ---- calc.  d2fdcdmu(i) (i=1~2)
      do i=1,2
        d2fdcdmu(i)=2.0d0*(dbdc(i)-dadc(i))+2.0d0*mu*(dadc(i)+dcdc(i)
     &                                             -2.0d0*dbdc(i))
      enddo
c
c                                          ---- calc.  d2fdc2(i) (i=1~2)
      do i=1,2
        d2fdc2(i)=d2adc2(i)+2.0d0*mu*(d2bdc2(i)-d2adc2(i))
     &            +mu*mu*(d2adc2(i)+d2cdc2(i)-2.0d0*d2bdc2(i))
      enddo
c
c                                                 ---- for d2phidx2(k,l)
c
      vcommon=1.0d0/(f(1)*dfdmu(2)-f(2)*dfdmu(1))
      vtmp1=dfdc(1)*dfdmu(2)+f(1)*d2fdmudc(2)
     &                           -dfdc(2)*dfdmu(1)-f(2)*d2fdmudc(1)
      vtmp2=dfdmu(1)*dfdmu(2)+f(1)*d2fdmu2(2)
     &                           -dfdmu(2)*dfdmu(1)-f(2)*d2fdmu2(1)
      vtmp3=d2fdcdmu(2)*dfdmu(1)+dfdc(2)*d2fdmu2(1)
     &                     -d2fdcdmu(1)*dfdmu(2)-dfdc(1)*d2fdmu2(2)
      vtmp4=d2fdc2(2)*dfdmu(1)+dfdc(2)*d2fdmudc(1)
     &                      -d2fdc2(1)*dfdmu(2)-dfdc(1)*d2fdmudc(2)
c
      va=vcommon
      vc=dfdc(2)*dfdmu(1)-dfdc(1)*dfdmu(2)
c
      do i=1,2
        dvadx(i)=-vtmp2*vcommon*vcommon*dmudx(i)
      enddo
        dvadx(3)=-vtmp1*vcommon*vcommon
     &           -vtmp2*vcommon*vcommon*dmudx(3)
c
      do i=1,3
        dsedx(i)=dphidx(i)
      enddo
c
      do i=1,2
        dvcdx(i)=vtmp3*dmudx(i)
      enddo
        dvcdx(3)=vtmp4+vtmp3*dmudx(3)
c
c                                    ---- calc.  d2phidx2(i,j) (i,j=1~3)
      do j=1,2
        d2phidx2(1,j)=(-dfdmu(2)*vtmp2*vcommon*vcommon
     &                             +d2fdmu2(2)*vcommon)*dmudx(j)
      enddo
        d2phidx2(1,3)=-dfdmu(2)*vtmp1*vcommon*vcommon
     &                +d2fdmudc(2)*vcommon+(-dfdmu(2)*vtmp2
     &                  *vcommon*vcommon+d2fdmu2(2)*vcommon)*dmudx(3)
c
      do j=1,2
        d2phidx2(2,j)=(dfdmu(1)*vtmp2*vcommon*vcommon
     &                             -d2fdmu2(1)*vcommon)*dmudx(j)
      enddo
        d2phidx2(2,3)=dfdmu(1)*vtmp1*vcommon*vcommon
     &                -d2fdmudc(1)*vcommon+(dfdmu(1)*vtmp2
     &                  *vcommon*vcommon-d2fdmu2(1)*vcommon)*dmudx(3)
c
      do j=1,3
        d2phidx2(3,j)=dvadx(j)*se*vc+va*dsedx(j)*vc+va*se*dvcdx(j)
      enddo
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     calc. d2seds2(i,j) (i,j=1~3)
c
      subroutine jancae_vegter_calc_d2seds2 (d2seds2,d2phidx2,se,a,b,c,
     &      mu,x,vcos2t,vsin2t,iareaflag,dxds,dphidx,isflag,s,dseds,
     &      pryld,ndyld)
c-----------------------------------------------------------------------
      implicit none
c
      integer, intent(in) :: iareaflag,isflag,ndyld
      real*8, intent(in) :: d2phidx2(3,3),se,a(2),b(2),c(2),mu,x(4),
     &            vcos2t,vsin2t,dxds(3,3),dphidx(3),s(3),dseds(3),
     &            pryld(ndyld)
      real*8, intent(out) :: d2seds2(3,3)
      real*8, parameter :: TOL1 = 0.996515d0       ! =44.9deg
      real*8, parameter :: TOL2 = 1.003497d0       ! =45.1deg
      real*8, parameter :: TOL3a= 1.0d-6
      real*8, parameter :: TOL3b= 1.0d0-TOL3a
      real*8, parameter :: TOL4 = 1.0d-7
      real*8, parameter :: TOL4a= 1.0d0-TOL4       !thera<0.012812deg
      real*8, parameter :: TOL4b=-1.0d0+TOL4       !thera>89.98719deg
      real*8 vtan,d2xds2(3,3,3),vx1x2
      integer i,j,k,l,iflag
c
c                             ---- if condition to apply numerical diff.
c
      if(iareaflag.eq.3.or.iareaflag.eq.4) then
      vtan=x(2)/x(1)
      end if
c
      if(iareaflag.eq.4.and.vtan.ge.TOL1.and.vtan.le.TOL2) then
      iflag=1
c
      else if(iareaflag.eq.3.and.vtan.ge.TOL1.and.vtan.le.TOL2) then
      iflag=2
c
      else if(abs(mu).le.TOL3a.or.abs(mu).ge.TOL3b) then
      iflag=3
c
      else if(vcos2t.ge.TOL4a.or.vcos2t.le.TOL4b) then
      iflag=4
c
      else
      vx1x2=x(1)-x(2)
      iflag=0
      endif
c
c                                    ---- set d2xds2(k,i,j)  (k,i,j=1~3)
c
      if (iflag.eq.0) then
c
      d2xds2(1,1,1)=0.5d0*(1.0d0-vcos2t*vcos2t)/vx1x2
      d2xds2(1,1,2)=0.5d0*(vcos2t*vcos2t-1.0d0)/vx1x2
      d2xds2(1,1,3)=-vcos2t*vsin2t/vx1x2
c
      d2xds2(1,2,1)=0.5d0*(vcos2t*vcos2t-1.0d0)/vx1x2
      d2xds2(1,2,2)=0.5d0*(1-vcos2t*vcos2t)/vx1x2
      d2xds2(1,2,3)=vcos2t*vsin2t/vx1x2
c
      d2xds2(1,3,1)=-vcos2t*vsin2t/vx1x2
      d2xds2(1,3,2)=vcos2t*vsin2t/vx1x2
      d2xds2(1,3,3)=-2.0d0*(vsin2t*vsin2t-1.0d0)/vx1x2
c
      d2xds2(2,1,1)=0.5d0*(vcos2t*vcos2t-1.0d0)/vx1x2
      d2xds2(2,1,2)=0.5d0*(1.0d0-vcos2t*vcos2t)/vx1x2
      d2xds2(2,1,3)=vcos2t*vsin2t/vx1x2
c
      d2xds2(2,2,1)=0.5d0*(1.0d0-vcos2t*vcos2t)/vx1x2
      d2xds2(2,2,2)=0.5d0*(vcos2t*vcos2t-1.0d0)/vx1x2
      d2xds2(2,2,3)=-vcos2t*vsin2t/vx1x2
c
      d2xds2(2,3,1)=vcos2t*vsin2t/vx1x2
      d2xds2(2,3,2)=-vcos2t*vsin2t/vx1x2
      d2xds2(2,3,3)=2.0d0*(vsin2t*vsin2t-1.0d0)/vx1x2
c
      d2xds2(3,1,1)=-3.0d0*vcos2t*vsin2t*vsin2t/(vx1x2*vx1x2)
      d2xds2(3,1,2)=3.0d0*vcos2t*vsin2t*vsin2t/(vx1x2*vx1x2)
      d2xds2(3,1,3)=2.0d0*vsin2t*(2.0d0-3.0d0*vsin2t*vsin2t)/
     &                                           (vx1x2*vx1x2)
c
      d2xds2(3,2,1)=3.0d0*vcos2t*vsin2t*vsin2t/(vx1x2*vx1x2)
      d2xds2(3,2,2)=-3.0d0*vcos2t*vsin2t*vsin2t/(vx1x2*vx1x2)
      d2xds2(3,2,3)=-2.0d0*vsin2t*(2.0d0-3.0d0*vsin2t*vsin2t)/
     &                                           (vx1x2*vx1x2)
c
      d2xds2(3,3,1)=2.0d0*vsin2t*(3.0d0*vcos2t*vcos2t-1.0d0)/
     &                                           (vx1x2*vx1x2)
      d2xds2(3,3,2)=-2.0d0*vsin2t*(3.0d0*vcos2t*vcos2t-1.0d0)/
     &                                           (vx1x2*vx1x2)
      d2xds2(3,3,3)=4.0d0*vcos2t*(3.0d0*vsin2t*vsin2t-1.0d0)/
     &                                           (vx1x2*vx1x2)
      end if
c
c                                      ---- calc. d2seds2(i,j) (i,j=1~3)
c
      call jancae_clear2 ( d2seds2,3,3 )
c
      if (iflag.ne.0) then
      call jancae_vegter_d2seds2n(d2seds2,s,dseds,pryld,ndyld,se)
c
      else
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              d2seds2(i,j)=d2seds2(i,j)
     &                    +d2phidx2(k,l)*dxds(l,j)*dxds(k,i)
            enddo
              d2seds2(i,j)=d2seds2(i,j)
     &                    +dphidx(k)*d2xds2(k,i,j)
          enddo
        enddo
      enddo
      endif
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     numerical differential for 2nd order differentials
c
      subroutine jancae_vegter_d2seds2n (d2seds2,s,dseds,
     &                                   pryld,ndyld,se)
c-----------------------------------------------------------------------
      implicit none
c
      integer, intent(in)  :: ndyld
      real*8, intent(out) :: d2seds2(3,3)
      real*8, intent(in)  :: dseds(3),pryld(ndyld),se,s(3)
      real*8, parameter :: delta=1.0d-3
      real*8 s0(3),sea,seb,a,b,seba,seaa,sebb,seab,se0,ss(3)
      integer j,k
c
      s0(:) = s(:)
      ss(:) = s(:)
        do j=1,3
          do k=1,3
           if ( j.eq.k ) then
             se0=se
             ss(j)=s0(j)-delta
             call jancae_vegter_yieldfunc(3,ss,sea,dseds,d2seds2,0,
     &                                    pryld,ndyld)
             ss(j)=s0(j)+delta
             call jancae_vegter_yieldfunc(3,ss,seb,dseds,d2seds2,0,
     &                                    pryld,ndyld)
             ss(j)=s0(j)
             a=(se0-sea)/delta
             b=(seb-se0)/delta
             d2seds2(j,k)=(b-a)/delta
           else
             ss(j)=s0(j)-delta
             ss(k)=s0(k)-delta
             call jancae_vegter_yieldfunc(3,ss,seaa,dseds,d2seds2,0,
     &                                    pryld,ndyld)
             ss(j)=s0(j)+delta
             ss(k)=s0(k)-delta
             call jancae_vegter_yieldfunc(3,ss,seba,dseds,d2seds2,0,
     &                                    pryld,ndyld)
             ss(j)=s0(j)-delta
             ss(k)=s0(k)+delta
             call jancae_vegter_yieldfunc(3,ss,seab,dseds,d2seds2,0,
     &                                    pryld,ndyld)
             ss(j)=s0(j)+delta
             ss(k)=s0(k)+delta
             call jancae_vegter_yieldfunc(3,ss,sebb,dseds,d2seds2,0,
     &                                    pryld,ndyld)
             ss(j)=s0(j)
             ss(k)=s0(k)
             a=(seba-seaa)/(2.0d0*delta)
             b=(sebb-seab)/(2.0d0*delta)
             d2seds2(j,k)=(b-a)/(2.0d0*delta)
           endif
          end do
         end do
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     calc. equivalent stress for d2seds2n
c
      subroutine jancae_vegter_yieldfunc (nttl,s,se,dseds,d2seds2,
     &                                    nreq,pryld,ndyld)
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: nttl,nreq,ndyld
      real*8, intent(in) :: pryld(ndyld),s(3)
      real*8, intent(out)  :: se,d2seds2(3,3),dseds(3)
c
      call jancae_vegter ( s,se,dseds,d2seds2,nreq,pryld,ndyld )
c
      return
      end
c
c
cc--------------------------------------------------------------(yld2000)
c     Barlat YLD2000 yield function and its dfferentials
c     ( IJP v.19(203) p1297-1319. )
c

      subroutine jancae_yld2000 ( s,se,dseds,d2seds2,nreq,
     &                            pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension s(3),dseds(3),d2seds2(3,3),pryld(ndyld)
      dimension a(8),am(2,3,3),x(2,2),y(2,3),phi(2),
     &          dsedphi(2),dphidx(2,2),
     &          dxdy(2,2,3),dyds(2,3,3),
     &          d2sedphi2(2,2),d2phidx2(2,2,2),
     &          d2xdy2(2,2,3,3)
c
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
      do i=1,8
        a(i)=pryld(i+1)
      enddo
      em=pryld(9+1)
c                                         ---- set linear transf. matrix
      call jancae_yld2000_2d_am ( a,am )
c                                                 ---- equivalent stress
      call jancae_yld2000_2d_xyphi ( s,em,am,x,y,phi )
      q=phi(1)+phi(2)
      if ( q.le.0.0 ) q=0.0
      se=(0.5d0*q)**(1.0d0/em)
c                                            ---- 1st order differential
      if ( nreq.ge.1 ) then
        call jancae_yld2000_2d_ds1 ( em,am,x,y,phi,
     &                               dsedphi,dphidx,
     &                               dxdy,dyds,se )
        call jancae_clear1 ( dseds,3 )
        do nd=1,2
          do m=1,2
            do k=1,3
              do i=1,3
                dseds(i)=dseds(i)+
     &                   dsedphi(      nd      )*
     &                      dphidx(    nd,m    )*
     &                          dxdy(  nd,m,k  )*
     &                            dyds(nd,  k,i)
              enddo
            enddo
          enddo
        enddo
      endif
c                                            ---- 2nd order differential
      if ( nreq.ge.2 ) then
        call jancae_yld2000_2d_ds2 ( phi,x,y,em,
     &                               d2sedphi2,d2phidx2,
     &                               d2xdy2,se )
        call jancae_clear2 ( d2seds2,3,3 )
        do i=1,3
        do j=1,3
          do nd1=1,2
          do nd2=1,2
            do k=1,2
            do l=1,2
              do m=1,3
              do n=1,3
                d2seds2(i,j)=d2seds2(i,j)+
     &                       d2sedphi2(     nd1,nd2  )*
     &                           dphidx(    nd1,k    )*
     &                               dxdy(  nd1,k,m  )*
     &                                 dyds(nd1,  m,i)*
     &                           dphidx(    nd2,l    )*
     &                               dxdy(  nd2,l,n  )*
     &                                 dyds(nd2,n  ,j)
              enddo
              enddo
            enddo
            enddo
          enddo
          enddo
          do nd=1,2
            do k=1,2
            do l=1,2
              do m=1,3
              do n=1,3
                d2seds2(i,j)=d2seds2(i,j)+
     &                        dsedphi(      nd      )*
     &                          d2phidx2(   nd,k,l  )*
     &                               dxdy(  nd,k,m  )*
     &                                 dyds(nd,  m,i)*
     &                               dxdy(  nd,l,n  )*
     &                                 dyds(nd,  n,j)
              enddo
              enddo
            enddo
            enddo
          enddo
          do nd=1,2
            do k=1,2
              do m=1,3
              do n=1,3
                d2seds2(i,j)=d2seds2(i,j)+
     &                       dsedphi(      nd      )*
     &                          dphidx(    nd,k    )*
     &                             d2xdy2( nd,k,m,n)*
     &                                dyds(nd,  m,i)*
     &                                dyds(nd,  n,j)
              enddo
              enddo
            enddo
          enddo
        enddo
        enddo
      endif
c
      return
      end
c
c
c
c--------------------------------------------------------------(yld2000)
c     set barlat-yld2k linear transformation matrix am
c
      subroutine jancae_yld2000_2d_am ( a,am )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(8),am(2,3,3)
c
c                                      ---- linear transformation matrix
      am(1,1,1)= 2.0d0*a(1)
      am(1,1,2)=-1.0d0*a(1)
      am(1,1,3)= 0.0
      am(1,2,1)=-1.0d0*a(2)
      am(1,2,2)= 2.0d0*a(2)
      am(1,2,3)= 0.0
      am(1,3,1)= 0.0
      am(1,3,2)= 0.0
      am(1,3,3)= 3.0d0*a(7)
c
      am(2,1,1)=-2.0d0*a(3)+2.0d0*a(4)+8.0d0*a(5)-2.0d0*a(6)
      am(2,1,2)=       a(3)-4.0d0*a(4)-4.0d0*a(5)+4.0d0*a(6)
      am(2,1,3)= 0.0
      am(2,2,1)= 4.0d0*a(3)-4.0d0*a(4)-4.0d0*a(5)+      a(6)
      am(2,2,2)=-2.0d0*a(3)+8.0d0*a(4)+2.0d0*a(5)-2.0d0*a(6)
      am(2,2,3)= 0.0
      am(2,3,1)= 0.0
      am(2,3,2)= 0.0
      am(2,3,3)= 9.0d0*a(8)
c
      do i=1,3
        do j=1,3
          am(1,i,j)=am(1,i,j)/3.0d0
          am(2,i,j)=am(2,i,j)/9.0d0
        enddo
      enddo
c
      return
      end
c
c
c
c--------------------------------------------------------------(yld2000)
c     calc. barlat-yld2k function x,y,phi
c
      subroutine jancae_yld2000_2d_xyphi ( s,em,am,x,y,phi )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension s(3),am(2,3,3),x(2,2),y(2,3),phi(2)
      dimension p(2)
c
      p(1)= 1.0d0
      p(2)=-1.0d0
c                                                       ---- {y}=[am]{s}
      call jancae_clear2 ( y,2,3 )
      do nd=1,2
        do i=1,3
          do j=1,3
            y(nd,i)=y(nd,i)+am(nd,i,j)*s(j)
          enddo
        enddo
      enddo
c                                        ---- {x}=principle value of {y}
      do nd=1,2
        a=(y(nd,1)-y(nd,2))**2.d0 + 4.0d0*y(nd,3)**2.d0
        a=sqrt(a)
        do i=1,2
          x(nd,i)=0.5d0*(y(nd,1)+y(nd,2)+p(i)*a)
        enddo
      enddo
c                                                 ---- phi(1) and phi(2)
      nd=1
      phi(nd)=abs(x(nd,1)-x(nd,2))**em
      nd=2
      phi(nd)=abs(2.0d0*x(nd,2)+x(nd,1))**em +
     &        abs(2.0d0*x(nd,1)+x(nd,2))**em
c
      return
      end
c
c
c
c--------------------------------------------------------------(yld2000)
c     set 1st order differential of parameters
c
      subroutine jancae_yld2000_2d_ds1 ( em,am,x,y,phi,
     &                                   dsedphi,dphidx,
     &                                   dxdy,dyds,se )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension am(2,3,3),x(2,2),y(2,3),phi(2),
     &          dsedphi(2),dphidx(2,2),
     &          dxdy(2,2,3),dyds(2,3,3)
c
      dimension p(2)
c
      eps=1.0d-16
c
      p(1)= 1.0d0
      p(2)=-1.0d0
      emi=1.0d0/em
c                                                          ---- dse/dphi
      q=phi(1)+phi(2)
      if ( q.le.0.0 ) q=0.0
      do nd=1,2
        dsedphi(nd)=(0.5d0**emi)*
     &               emi*q**(emi-1.0d0)
      enddo
c                                                           ---- dphi/dx
      nd=1
      a0=x(nd,1)-x(nd,2)
      b0=abs(a0)
      sgn0=0
      if ( b0.ge.eps*se ) sgn0=a0/b0
      dphidx(nd,1)= em*b0**(em-1.0d0) * sgn0
      dphidx(nd,2)=-em*b0**(em-1.0d0) * sgn0
c
      nd=2
      a1=2.0d0*x(nd,1)+      x(nd,2)
      a2=      x(nd,1)+2.0d0*x(nd,2)
      b1=abs(a1)
      b2=abs(a2)
      sgn1=0.0
      sgn2=0.0
      if ( b1.ge.eps*se )  sgn1=a1/b1
      if ( b2.ge.eps*se )  sgn2=a2/b2
      dphidx(nd,1)=em*(2.0d0*b1**(em-1.0d0) *sgn1 +
     &                       b2**(em-1.0d0) *sgn2  )
      dphidx(nd,2)=em*(      b1**(em-1.0d0) *sgn1 +
     &                 2.0d0*b2**(em-1.0d0) *sgn2  )
c
      do nd=1,2
        a=(y(nd,1)-y(nd,2))*(y(nd,1)-y(nd,2))+
     &    4.0d0*   y(nd,3) *         y(nd,3)
        a=sqrt(a)
        if ( a.gt.eps*se ) then
          do j=1,2
            dxdy(nd,j,1)=0.5d0*(1.0d0+p(j)*(y(nd,1)-y(nd,2))/a)
            dxdy(nd,j,2)=0.5d0*(1.0d0-p(j)*(y(nd,1)-y(nd,2))/a)
            dxdy(nd,j,3)=2.0d0*       p(j)* y(nd,3)         /a
          enddo
        else
          do j=1,2
            dxdy(nd,j,1)=0.5d0*(1.0d0+0.0)
            dxdy(nd,j,2)=0.5d0*(1.0d0-0.0)
            dxdy(nd,j,3)=2.0d0*       0.0
          enddo
        endif
      enddo
c
      do nd=1,2
        do i=1,3
          do j=1,3
            dyds(nd,i,j)=am(nd,i,j)
          enddo
        enddo
      enddo
c
      return
      end
c
c
c
c--------------------------------------------------------------(yld2000)
c     set 2nd order differential of parameters
c
      subroutine jancae_yld2000_2d_ds2 ( phi,x,y,em,
     &                                   d2sedphi2,d2phidx2,
     &                                   d2xdy2,se )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension phi(2),x(2,2),y(2,3),
     &          d2sedphi2(2,2),d2phidx2(2,2,2),
     &          d2xdy2(2,2,3,3)
      dimension p(2)
c
      eps=1.0d-16
c
      p(1)= 1.0d0
      p(2)=-1.0d0
      emi=1.0d0/em
c                                                        ---- d2se/dphi2
      q=phi(1)+phi(2)
      if ( q.le.0.0 ) q=0.0
      do nd1=1,2
        do nd2=1,2
          a=0.5d0**emi * emi * (emi-1.0d0) *
     &      q**(emi-2.0d0)
          d2sedphi2(nd1,nd2)=a
        enddo
      enddo
c                                                         ---- d2phi/dx2
      nd=1
      do i=1,2
        do j=1,2
          a=(em-1.0d0)*em*(abs(x(nd,1)-x(nd,2)))**(em-2.0d0)
          if ( i.ne.j ) a=-a
          d2phidx2(nd,i,j)=a
        enddo
      enddo
      nd=2
      do i=1,2
        do j=1,2
          if ( i.eq.j ) then
            if ( i.eq.1 ) then
              a=(em-1.0d0)*em*
     &        (4.0d0*(abs(2.0d0*x(nd,1)+      x(nd,2)))**(em-2.0d0)+
     &               (abs(      x(nd,1)+2.0d0*x(nd,2)))**(em-2.0d0) )
            else
              a=(em-1.0d0)*em*
     &        (      (abs(2.0d0*x(nd,1)+      x(nd,2)))**(em-2.0d0)+
     &         4.0d0*(abs(      x(nd,1)+2.0d0*x(nd,2)))**(em-2.0d0) )
            endif
          else
            a=  (em-1.0d0)*em*
     &        (2.0d0*(abs(2.0d0*x(nd,1)+      x(nd,2)))**(em-2.0d0)+
     &         2.0d0*(abs(      x(nd,1)+2.0d0*x(nd,2)))**(em-2.0d0) )
          endif
          d2phidx2(nd,i,j)=a
        enddo
      enddo
c                                                           ---- d2x/dy2
      do nd=1,2
        a= (y(nd,1)-y(nd,2))*(y(nd,1)-y(nd,2))+
     &     4.0d0*   y(nd,3) *         y(nd,3)
        if ( a.gt.eps*se ) then
          a=1.0d0/sqrt(a**3)
          do m=1,2
            do i=1,3
              do j=1,3
                ij=i*10+j
                if ( (ij.eq.11).or.(ij.eq.22) ) then
                  q= y(nd,3)         * y(nd,3)
                else if ( ij.eq.33 ) then
                  q=(y(nd,1)-y(nd,2))*(y(nd,1)-y(nd,2))
                else if ( (ij.eq.12).or.(ij.eq.21) ) then
                  q=-y(nd,3)         * y(nd,3)
                else if ( (ij.eq.23).or.(ij.eq.32) ) then
                  q= y(nd,3)         *(y(nd,1)-y(nd,2))
                else
                  q=-y(nd,3)         *(y(nd,1)-y(nd,2))
                endif
                d2xdy2(nd,m,i,j)=2.0d0*a*p(m)*q
              enddo
            enddo
          enddo
        else
          do m=1,2
            do i=1,3
              do j=1,3
                d2xdy2(nd,m,i,j)=0.0
              enddo
            enddo
          enddo
        endif
      enddo
c
      return
      end
c
c
c
c----------------------------------------------------------(yld2004-18p)
c     Barlat Yld2004-18p yield function and its dfferentials
c     ( IJP v.21(2005) p1009-1039. )
c
      subroutine jancae_yld2004_18p ( s,se,dseds,d2seds2,nreq,
     &                                pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension s(6),dseds(6),d2seds2(6,6),pryld(ndyld)
c
      dimension sp1(6),sp2(6),cp1(6,6),cp2(6,6),cl(6,6),ctp1(6,6)
      dimension ctp2(6,6),psp1(3),psp2(3),hp1(3),hp2(3)
      dimension dfadpsp1(3),dfadpsp2(3),dpsdhp1(3,3),dpsdhp2(3,3)
      dimension dfadhp1(3),dfadhp2(3)
      dimension dhdsp1(3,6),dhdsp2(3,6),dsdsp1(6,6), dsdsp2(6,6)
      dimension dfads(6)
      dimension d2fadpsp11(3,3),d2fadpsp22(3,3)
      dimension d2fadpsp12(3,3),d2fadpsp21(3,3)
      dimension d2fadhp11(3,3),d2fadhp22(3,3)
      dimension d2fadhp12(3,3),d2fadhp21(3,3)
      dimension d2psdhp11(3,3,3),d2psdhp22(3,3,3)
      dimension d2hdsp11(3,6,6),d2hdsp22(3,6,6)
      dimension d2fads2(6,6)
      dimension delta(3,3)
      dimension xx1(3,6),xx2(3,6)
c
c                                                         ---- variables
c
c     sp1(6),sp2(6) : linear-transformed stress
c     cp1(6,6),cp2(6,6) : matrix for anisotropic parameters
c     cl(6,6) : matrix for transforming Cauchy stress to deviator
c     ctp1(6,6) : matrix for transforming Cauchy stress to sp1
c     ctp2(6,6) : matrix for transforming Cauchy stress to sp2
c     dc : coefficient of equivalent stress(se=(fai/dc)**(1/am))
c     aap1(3),aap2(3) : for dc
c     ppp1,ppp2 : for dc
c     qqp1,qqp2 : for dc
c     ttp1,ttp2 : for dc
c     bbp1(3),bbp2(3) : for dc
c     fai : yield fuction
c     psp1(3) : principal values of sp1
c     psp2(3) : principal values of sp2
c     hp1(3) : invariants of sp1
c     hp2(3) : invariants of sp2
c     cep1,cep2 : coefficient of characteristic equation p
c     ceq1,ceq2 : coefficient of characteristic equation q
c     cet1,cet2 : arccos ( q/p^(3/2) )
c
c     dfadpsp1(3),dfadpsp2(3) : d(fai)/d(psp)
c     dfadhp1(3),dfadhp2(3) : d(fai)/d(hp)
c     dpsdhp1(3,3),dpsdhp2(3,3) : d(psp)/d(hp)
c     dhdsp1(3,6),dhdsp2(3,6) : d(hp)/d(sp)
c     dsdsp1(6,6), dsdsp2(6,6) : d(sp)/d(s)
c     dfads(6) : d(fai)/d(s)
c
c     d2fadpsp11(3,3),d2fadpsp22(3,3) : d2(fai)/d(psp)2
c     d2fadpsp12(3,3),d2fadpsp21(3,3) : d2(fai)/d(psp)2
c     d2psdhp11(3,3,3),d2psdhp22(3,3,3) :d2(psp)/d(hp)2
c     d2hdsp11(3,6,6),d2hdsp22(3,6,6) : d2(hp)/d(sp)2
c     d2fads2(6,6) : d2(fai)/d(s)2
c
c     eps2,eps3 : values for calculation of d(psp)/d(hp),d2(psp)/d(hp)2
c
      pi=acos(-1.0d0)
      eps2 = 1.0d-15
      eps3 = 1.0d-8
      del = 1.0d-4
c                                                   ---- Kronecker Delta
      call jancae_clear2( delta,3,3 )
      do i=1,3
        delta(i,i)=1.0d0
      end do
c                                        ---- set anisotropic parameters
      call jancae_clear2( cp1,6,6 )
      call jancae_clear2( cp2,6,6 )
      cp1(1,2) = -pryld(1+1)
      cp1(1,3) = -pryld(1+2)
      cp1(2,1) = -pryld(1+3)
      cp1(2,3) = -pryld(1+4)
      cp1(3,1) = -pryld(1+5)
      cp1(3,2) = -pryld(1+6)
      cp1(4,4) =  pryld(1+9)  ! for tau_xy (c'66 in original paper)
      cp1(5,5) =  pryld(1+8)  ! for tau_xz (c'55 in original paper)
      cp1(6,6) =  pryld(1+7)  ! for tau_zx (c'44 in original paper)
      cp2(1,2) = -pryld(1+10)
      cp2(1,3) = -pryld(1+11)
      cp2(2,1) = -pryld(1+12)
      cp2(2,3) = -pryld(1+13)
      cp2(3,1) = -pryld(1+14)
      cp2(3,2) = -pryld(1+15)
      cp2(4,4) =  pryld(1+18) ! for tau_xy (c"66 in original paper)
      cp2(5,5) =  pryld(1+17) ! for tau_xz (c"55 in original paper)
      cp2(6,6) =  pryld(1+16) ! for tau_zx (c"44 in original paper)
      am       =  pryld(1+19)
c      dc       =  4
      ami=1.0d0/am
c
c             ---- set matrix for transforming Cauchy stress to deviator
      call jancae_clear2( cl,6,6 )
      do i=1,3
        do j=1,3
          if(i==j) then
            cl(i,j)=2.0d0
          else
            cl(i,j)=-1.0d0
          end if
        end do
      end do
      do i=4,6
        cl(i,i) = 3.0d0
      end do
      do i=1,6
        do j=1,6
          cl(i,j)=cl(i,j)/3.0d0
        end do
      end do
c
c                  ---- matrix for transforming Cauchy stress to sp1,sp2
      call jancae_mm (ctp1,cp1,cl,6,6,6)
      call jancae_mm (ctp2,cp2,cl,6,6,6)
c                               ---- coefficient of equivalent stress dc
      call jancae_yld2004_18p_coef (cp1,cp2,pi,am,dc)
c                                  ---- calculation of equivalent stress
      call jancae_yld2004_18p_yf (ctp1,ctp2,s,am,ami,dc,pi,
     1       sp1,sp2,psp1,psp2,hp1,hp2,cetpq1,cetpq2,fai,se)
c
c                                            ---- 1st order differential
      if ( nreq.ge.1 ) then
c                                                     ---- d(fai)/d(psp)
        do i = 1,3
          dfadpsp1(i) = 0.0d0
          dfadpsp2(i) = 0.0d0
          do j= 1,3
            dfadpsp1(i)=dfadpsp1(i)+(psp1(i)-psp2(j))
     1                *abs(psp1(i)-psp2(j))**(am-2.0d0)
            dfadpsp2(i)=dfadpsp2(i)+(psp1(j)-psp2(i))
     1                *abs(psp1(j)-psp2(i))**(am-2.0d0)
          end do
          dfadpsp1(i)=dfadpsp1(i)*am
          dfadpsp2(i)=dfadpsp2(i)*(-am)
        end do
c                                         ---- d(psp)/d(hp)&d(fai)/d(hp)
        call jancae_clear2( dpsdhp1,3,3 )
        call jancae_clear2( dpsdhp2,3,3 )
        call jancae_clear1( dfadhp1,3 )
        call jancae_clear1( dfadhp2,3 )
c                                                  ---- theta'<>0 & <>pi
        if(abs(cetpq1-1.0d0)>=eps2 .and. abs(cetpq1+1.0d0)>=eps2) then
          do i=1,3
            call jancae_yld2004_18p_dpsdhp(i,psp1,hp1,dpsdhp1)
          end do
c                                                          ---- theta'=0
        else if(abs(cetpq1-1.0d0)<eps2) then
          i=1
          call jancae_yld2004_18p_dpsdhp(i,psp1,hp1,dpsdhp1)
          do i=2,3
            do j=1,3
              dpsdhp1(i,j)=-0.5d0*(dpsdhp1(1,j)-3.0d0*delta(1,j))
            end do
           end do
c                                                         ---- theta'=pi
        else
          i=3
          call jancae_yld2004_18p_dpsdhp(i,psp1,hp1,dpsdhp1)
          do i=1,2
            do j=1,3
              dpsdhp1(i,j)=-0.5d0*(dpsdhp1(3,j)-3.0d0*delta(1,j))
            end do
           end do
        end if
c                                                 ---- theta''<>0 & <>pi
        if(abs(cetpq2-1.0d0)>=eps2 .and. abs(cetpq2+1.0d0)>=eps2) then
          do i=1,3
            call jancae_yld2004_18p_dpsdhp(i,psp2,hp2,dpsdhp2)
          end do
c                                                         ---- theta''=0
        else if(abs(cetpq2-1.0d0)<eps2)then
          i=1
          call jancae_yld2004_18p_dpsdhp(i,psp2,hp2,dpsdhp2)
          do i=2,3
            do j=1,3
              dpsdhp2(i,j)=-0.5d0*(dpsdhp2(1,j)-3.0d0*delta(1,j))
            end do
           end do
c                                                        ---- theta''=pi
        else
          i=3
          call jancae_yld2004_18p_dpsdhp(i,psp2,hp2,dpsdhp2)
          do i=1,2
            do j=1,3
              dpsdhp2(i,j)=-0.5d0*(dpsdhp2(3,j)-3.0d0*delta(1,j))
            end do
           end do
        end if
c
        do i=1,3
          do j=1,3
            dfadhp1(i)=dfadhp1(i)+dfadpsp1(j)*dpsdhp1(j,i)
            dfadhp2(i)=dfadhp2(i)+dfadpsp2(j)*dpsdhp2(j,i)
          end do
        end do
c                                                       ---- d(hp)/d(sp)
        call jancae_clear2( dhdsp1,3,6 )
        call jancae_clear2( dhdsp2,3,6 )
        do i = 1,3
          j=mod(i,3)+1
          k=mod(i+1,3)+1
          l=mod(i,3)+4
          if ( i == 1 ) l=6
          if ( i == 2 ) l=5
          dhdsp1(1,i)=1.0d0/3.0d0
          dhdsp2(1,i)=1.0d0/3.0d0
          dhdsp1(2,i)=-1.0d0/3.0d0*(sp1(j)+sp1(k))
          dhdsp2(2,i)=-1.0d0/3.0d0*(sp2(j)+sp2(k))
          dhdsp1(3,i)=1.0d0/2.0d0*(sp1(j)*sp1(k)-sp1(l)**2)
          dhdsp2(3,i)=1.0d0/2.0d0*(sp2(j)*sp2(k)-sp2(l)**2)
        end do
        do i=4,6
          k=mod(i+1,3)+1
          l=mod(i,3)+4
          m=mod(i+1,3)+4
          if ( i == 5 ) k=2
          if ( i == 6 ) k=1
          dhdsp1(2,i)=2.0d0/3.0d0*sp1(i)
          dhdsp2(2,i)=2.0d0/3.0d0*sp2(i)
          dhdsp1(3,i)=sp1(l)*sp1(m)-sp1(k)*sp1(i)
          dhdsp2(3,i)=sp2(l)*sp2(m)-sp2(k)*sp2(i)
        end do
c                                                        ---- d(sp)/d(s)
        do i=1,6
          do j=1,6
            dsdsp1(i,j)=ctp1(i,j)
            dsdsp2(i,j)=ctp2(i,j)
          end do
        end do
c                                                       ---- d(fai)/d(s)
        call jancae_clear1( dfads,6 )
        call jancae_clear2( xx1,3,6 )
        call jancae_clear2( xx2,3,6 )
        do l=1,6
          do j=1,3
            do k=1,6
              xx1(j,l)=xx1(j,l)+dhdsp1(j,k)*dsdsp1(k,l)
              xx2(j,l)=xx2(j,l)+dhdsp2(j,k)*dsdsp2(k,l)
            end do
            dfads(l)=dfads(l)
     1               +dfadhp1(j)*xx1(j,l)
     2               +dfadhp2(j)*xx2(j,l)
          end do
        end do
c                                                        ---- d(se)/d(s)
        dsedfa=fai**(ami-1.0d0)/am/dc**ami
        do i=1,6
          dseds(i)=dsedfa*dfads(i)
        end do
c
      endif
c                                            ---- 2nd order differential
      if ( nreq.ge.2 ) then
c                                                   ---- d2(fai)/d(psp)2
        call jancae_clear2( d2fadpsp11,3,3 )
        call jancae_clear2( d2fadpsp22,3,3 )
        call jancae_clear2( d2fadpsp12,3,3 )
        call jancae_clear2( d2fadpsp21,3,3 )
        do i=1,3
          d2fadpsp11(i,i)=am*(am-1.0d0)
     1                   *(abs(psp1(i)-psp2(1))**(am-2.0d0)
     2                   + abs(psp1(i)-psp2(2))**(am-2.0d0)
     3                   + abs(psp1(i)-psp2(3))**(am-2.0d0))
          d2fadpsp22(i,i)=am*(am-1.0d0)
     1                   *(abs(psp1(1)-psp2(i))**(am-2.0d0)
     2                   + abs(psp1(2)-psp2(i))**(am-2.0d0)
     3                   + abs(psp1(3)-psp2(i))**(am-2.0d0))
          do j=1,3
            d2fadpsp12(i,j)=-am*(am-1.0d0)
     1                     * abs(psp1(i)-psp2(j))**(am-2.0d0)
            d2fadpsp21(i,j)=-am*(am-1.0d0)
     1                     * abs(psp1(j)-psp2(i))**(am-2.0d0)
          end do
        end do
c                                                    ---- d2(psp)/d(hp)2
        call jancae_clear3( d2psdhp11,3,3,3 )
        call jancae_clear3( d2psdhp22,3,3,3 )
c
        if(abs(cetpq1-1.0d0)>=eps3 .and. abs(cetpq2-1.0d0)>=eps3 .and.
     1     abs(cetpq1+1.0d0)>=eps3 .and. abs(cetpq2+1.0d0)>=eps3) then
c
          do i=1,3
            call jancae_yld2004_18p_d2psdhp(i,psp1,hp1,d2psdhp11)
            call jancae_yld2004_18p_d2psdhp(i,psp2,hp2,d2psdhp22)
          end do
c                                                    ---- d2(fai)/d(hp)2
          call jancae_clear2(d2fadhp11,3,3)
          call jancae_clear2(d2fadhp12,3,3)
          call jancae_clear2(d2fadhp21,3,3)
          call jancae_clear2(d2fadhp22,3,3)
c                                                 -- d2(fai)/d(hd)d(hd)
          do iq=1,3
            do m=1,3
              do ip=1,3
                do l=1,3
                  d2fadhp11(iq,m)=d2fadhp11(iq,m)
     1                   +d2fadpsp11(ip,l)*dpsdhp1(l,m)*dpsdhp1(ip,iq)
                end do
                d2fadhp11(iq,m)=d2fadhp11(iq,m)
     1                   +dfadpsp1(ip)*d2psdhp11(ip,iq,m)
              end do
            end do
          end do
c                                              ---- d2(fai)/d(hdd)d(hdd)
          do iq=1,3
            do m=1,3
              do ip=1,3
                do l=1,3
                  d2fadhp22(iq,m)=d2fadhp22(iq,m)
     1                   +d2fadpsp22(ip,l)*dpsdhp2(l,m)*dpsdhp2(ip,iq)
                end do
                d2fadhp22(iq,m)=d2fadhp22(iq,m)
     1                   +dfadpsp2(ip)*d2psdhp22(ip,iq,m)
              end do
            end do
          end do
c                         ---- d2(fai)/d(hdd)d(hd) & d2(fai)/d(hd)d(hdd)
          do iq=1,3
            do m=1,3
              do ip=1,3
                do l=1,3
                  d2fadhp12(iq,m)=d2fadhp12(iq,m)
     1                   +d2fadpsp12(ip,l)*dpsdhp2(l,m)*dpsdhp1(ip,iq)
                  d2fadhp21(iq,m)=d2fadhp21(iq,m)
     1                   +d2fadpsp21(ip,l)*dpsdhp1(l,m)*dpsdhp2(ip,iq)
                end do
               end do
            end do
          end do
c                                                     ---- d2(hp)/d(sp)2
          call jancae_clear3( d2hdsp11,3,6,6)
          call jancae_clear3( d2hdsp22,3,6,6)
          do i=1,3
            j=mod(i,3)+1
            k=mod(i+1,3)+1
            dummy=-1.0d0/3.0d0
            d2hdsp11(2,i,j)=dummy
            d2hdsp11(2,j,i)=dummy
            d2hdsp22(2,i,j)=dummy
            d2hdsp22(2,j,i)=dummy
            d2hdsp11(3,i,j)=sp1(k)/2.0d0
            d2hdsp11(3,j,i)=sp1(k)/2.0d0
            d2hdsp22(3,i,j)=sp2(k)/2.0d0
            d2hdsp22(3,j,i)=sp2(k)/2.0d0
          end do
          do i=4,6
            j=mod(i,3)+4
            k=mod(i+1,3)+1
            m=mod(i+1,3)+4
            if ( i == 5 ) k = 2
            if ( i == 6 ) k = 1
            dummy=2.0d0/3.0d0
            d2hdsp11(2,i,i)=dummy
            d2hdsp22(2,i,i)=dummy
            d2hdsp11(3,i,i)=-sp1(k)
            d2hdsp22(3,i,i)=-sp2(k)
            d2hdsp11(3,i,j)=sp1(m)
            d2hdsp11(3,j,i)=sp1(m)
            d2hdsp22(3,i,j)=sp2(m)
            d2hdsp22(3,j,i)=sp2(m)
          end do
          do i=1,3
            j=mod(i,3)+4
            if ( i == 1 ) j = 6
            if ( i == 2 ) j = 5
            d2hdsp11(3,i,j)=-sp1(j)
            d2hdsp11(3,j,i)=-sp1(j)
            d2hdsp22(3,i,j)=-sp2(j)
            d2hdsp22(3,j,i)=-sp2(j)
          end do
c                                                     ---- d2(fai)/d(s)2
          call jancae_clear2( d2fads2,6,6)
          call jancae_clear2( xx1,3,6 )
          call jancae_clear2( xx2,3,6 )
          do i=1,3
            do j=1,6
              do ip=1,6
                xx1(i,j)=xx1(i,j)+dhdsp1(i,ip)*dsdsp1(ip,j)
                xx2(i,j)=xx2(i,j)+dhdsp2(i,ip)*dsdsp2(ip,j)
              end do
            end do
          end do
          do i=1,6
            do j=1,6
              do iq=1,3
                do m=1,3
                  d2fads2(i,j)=d2fads2(i,j)
     1                       +d2fadhp11(iq,m)*xx1(iq,i)*xx1(m,j)
     2                       +d2fadhp12(iq,m)*xx1(iq,i)*xx2(m,j)
     3                       +d2fadhp21(iq,m)*xx2(iq,i)*xx1(m,j)
     4                       +d2fadhp22(iq,m)*xx2(iq,i)*xx2(m,j)
                end do
                do n=1,6
                  do ir=1,6
                    d2fads2(i,j)=d2fads2(i,j)
     1                       +d2hdsp11(iq,ir,n)*dfadhp1(iq)
     2                                  *dsdsp1(ir,i)*dsdsp1(n,j)
     3                       +d2hdsp22(iq,ir,n)*dfadhp2(iq)
     4                                  *dsdsp2(ir,i)*dsdsp2(n,j)
                  end do
                end do
              end do
            end do
          end do
c                                                      ---- d2(se)/d(s)2
          d2sedfa2=ami*(ami-1.0d0)*fai**(ami-2.0d0)/dc**ami
          do i=1,6
            do j=1,6
              d2seds2(i,j)=d2sedfa2*dfads(i)*dfads(j)
     1                  +dsedfa*d2fads2(i,j)
            end do
          end do
        else
c                                            ---- numerical differential
          call jancae_yld2004_18p_nu2 ( ctp1,ctp2,s,se,am,ami,
     1                                  dc,pi,del,d2seds2)
        end if
      endif
c
      return
      end
c
c
c
c----------------------------------------------------------(yld2004-18p)
c      calculate  coefficient of equivalent stress dc
c
      subroutine jancae_yld2004_18p_coef (cp1,cp2,pi,am,dc)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension cp1(6,6),cp2(6,6),bbp1(3),bbp2(3)
c
      call jancae_yld2004_18p_coef_sub(cp1,pi,bbp1)
      call jancae_yld2004_18p_coef_sub(cp2,pi,bbp2)
      dc=0.0d0
      do i=1,3
        do j=1,3
          dc=dc+(abs(bbp1(i)-bbp2(j)))**am
        end do
      end do
      return
      end
c
c
c
c----------------------------------------------------------(yld2004-18p)
c     calculate  coefficient of equivalent stress dc 2
c
      subroutine jancae_yld2004_18p_coef_sub (cp,pi,bbp)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension cp(6,6),aap(3),bbp(3)
c
c                                                  ---- coefficients aap
      aap(1)=(cp(1,2)+cp(1,3)-2.0d0*cp(2,1)
     1        +cp(2,3)-2.0d0*cp(3,1)+cp(3,2))/9.0d0
      aap(2)=((2.0d0*cp(2,1)-cp(2,3))*(cp(3,2)-2.0d0*cp(3,1))
     2         +(2.0d0*cp(3,1)-cp(3,2))*(cp(1,2)+cp(1,3))
     3         +(cp(1,2)+cp(1,3))*(2.0d0*cp(2,1)-cp(2,3)))/2.7d1
      aap(3)=(cp(1,2)+cp(1,3))*(cp(2,3)-2.0d0*cp(2,1))
     1        *(cp(3,2)-2.0d0*cp(3,1))/5.4d1
c                                          ---- coefficients ppp,qqp,ttp
      ppp=aap(1)**2+aap(2)
      qqp=(2.0d0*aap(1)**3+3.0d0*aap(1)*aap(2)+2.0d0*aap(3))/2.0d0
      ttp=acos(qqp/ppp**(3.0d0/2.0d0))
c                                                  ---- coefficients bbp
      bbp(1)=2.0d0*sqrt(ppp)*cos(ttp/3.0d0)+aap(1)
      bbp(2)=2.0d0*sqrt(ppp)*cos((ttp+4.0d0*pi)/3.0d0)+aap(1)
      bbp(3)=2.0d0*sqrt(ppp)*cos((ttp+2.0d0*pi)/3.0d0)+aap(1)
c
      return
      end
c
c
c
c----------------------------------------------------------(yld2004-18p)
c     calculate yield function
c
      subroutine jancae_yld2004_18p_yf (ctp1,ctp2,s,am,ami,dc,pi,
     1         sp1,sp2,psp1,psp2,hp1,hp2,cetpq1,cetpq2,fai,se)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension ctp1(6,6),ctp2(6,6),s(6)
      dimension sp1(6),sp2(6),psp1(3),psp2(3),hp1(3),hp2(3)
c
      call jancae_yld2004_18p_yfsub(ctp1,s,pi,
     1                              sp1,psp1,hp1,cetpq1)
      call jancae_yld2004_18p_yfsub(ctp2,s,pi,
     1                              sp2,psp2,hp2,cetpq2)
c                                                    ---- yield function
      fai=0.0d0
      do i=1,3
        do j=1,3
          fai=fai+(abs(psp1(i)-psp2(j)))**am
        end do
      end do
c                                                 ---- equivalent stress
      se=(fai/dc)**ami
c
      return
      end
c
c
c
c----------------------------------------------------------(yld2004-18p)
c     calculate yield function2
c
      subroutine jancae_yld2004_18p_yfsub(ctp,s,pi,
     1                                    sp,psp,hp,cetpq)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension ctp(6,6),s(6)
      dimension sp(6),psp(3),hp(3)
c
c                                         ---- linear-transformed stress
      call jancae_mv (sp,ctp,s,6,6)
c                                                  ---- invariants of sp
      hp(1)=(sp(1)+sp(2)+sp(3))/3.0d0
      hp(2)=(sp(5)**2+sp(6)**2+sp(4)**2
     1       -sp(2)*sp(3)-sp(3)*sp(1)-sp(1)*sp(2))/3.0d0
      hp(3)=(2.0d0*sp(5)*sp(6)*sp(4)+sp(1)*sp(2)*sp(3)
     1       -sp(1)*sp(6)**2-sp(2)*sp(5)**2-sp(3)*sp(4)**2)/2.0d0
c                           ---- coefficients of characteristic equation
      hpq=sqrt(hp(1)**2+hp(2)**2+hp(3)**2)
      if ( hpq.gt.1.0e-16 ) then
        cep=hp(1)**2+hp(2)
        ceq=(2.0d0*hp(1)**3+3.0d0*hp(1)*hp(2)+2.0d0*hp(3))/2.0d0
        cetpq=ceq/cep**(3.0d0/2.0d0)
        if ( cetpq >  1.0d0 ) cetpq= 1.0d0
        if ( cetpq < -1.0d0 ) cetpq=-1.0d0
        cet=acos(cetpq)
c                                           ---- principal values of sp1
        psp(1)=2.0d0*sqrt(cep)*cos(cet/3.0d0)+hp(1)
        psp(2)=2.0d0*sqrt(cep)*cos((cet+4.0d0*pi)/3.0d0)+hp(1)
        psp(3)=2.0d0*sqrt(cep)*cos((cet+2.0d0*pi)/3.0d0)+hp(1)
      else
        cetpq=0.0
        do i=1,3
           psp(i)=0.0
        enddo
      endif
c
      return
      end
c
c
c
c----------------------------------------------------------(yld2004-18p)
c     numerical differential for 2nd order differentials
c
      subroutine jancae_yld2004_18p_nu2 (ctp1,ctp2,s,se,am,ami,
     1                                   dc,pi,del,d2seds2)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension s(6),ctp1(6,6),ctp2(6,6),d2seds2(6,6)
      dimension sp1(6),sp2(6),psp1(3),psp2(3),hp1(3),hp2(3),s0(6)
c
      s0(:) = s(:)
      do i=1, 6
        do j=1, 6
          if(i == j) then
            s0(i)=s(i)-del
            call jancae_yld2004_18p_yf (ctp1,ctp2,s0,am,ami,dc,pi,
     1         sp1,sp2,psp1,psp2,hp1,hp2,cetpq1,cetpq2,fai,sea)
            s0(i)=s(i)+del
            call jancae_yld2004_18p_yf (ctp1,ctp2,s0,am,ami,dc,pi,
     1         sp1,sp2,psp1,psp2,hp1,hp2,cetpq1,cetpq2,fai,seb)
            s0(i)=s(i)
            abc1=(se-sea)/del
            abc2=(seb-se)/del
            d2seds2(i,j) = (abc2-abc1)/del
          else
            s0(i)=s(i)-del
            s0(j)=s(j)-del
            call jancae_yld2004_18p_yf (ctp1,ctp2,s0,am,ami,dc,pi,
     1         sp1,sp2,psp1,psp2,hp1,hp2,cetpq1,cetpq2,fai,seaa)
            s0(i)=s(i)+del
            s0(j)=s(j)-del
            call jancae_yld2004_18p_yf (ctp1,ctp2,s0,am,ami,dc,pi,
     1         sp1,sp2,psp1,psp2,hp1,hp2,cetpq1,cetpq2,fai,seba)
            s0(i)=s(i)-del
            s0(j)=s(j)+del
            call jancae_yld2004_18p_yf (ctp1,ctp2,s0,am,ami,dc,pi,
     1         sp1,sp2,psp1,psp2,hp1,hp2,cetpq1,cetpq2,fai,seab)
            s0(i)=s(i)+del
            s0(j)=s(j)+del
            call jancae_yld2004_18p_yf (ctp1,ctp2,s0,am,ami,dc,pi,
     1         sp1,sp2,psp1,psp2,hp1,hp2,cetpq1,cetpq2,fai,sebb)
            s0(i)=s(i)
            s0(j)=s(j)
            abc1 = (seba - seaa)/(2.0d0*del)
            abc2 = (sebb - seab)/(2.0d0*del)
            d2seds2(i,j) = (abc2 - abc1)/(2.0d0*del)
          end if
        end do
      end do
c
      return
      end
c
c
c
c----------------------------------------------------------(yld2004-18p)
c     calculate d(psp)/d(hp)
c
      subroutine jancae_yld2004_18p_dpsdhp(i,psp,hp,dpsdhp)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension psp(3),hp(3),dpsdhp(3,3)
c
      dummy=psp(i)**2-2.0d0*hp(1)*psp(i)-hp(2)
      dpsdhp(i,1)=psp(i)**2/dummy
      dpsdhp(i,2)=psp(i)/dummy
      dpsdhp(i,3)=2.0d0/3.0d0/dummy
c
      return
      end
c
c
c
c----------------------------------------------------------(yld2004-18p)
c     calculate d2(psp)/d(hp)2
c
      subroutine jancae_yld2004_18p_d2psdhp(i,psp,hp,d2psdhp)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension psp(3),hp(3),d2psdhp(3,3,3)
c
      dummy=(psp(i)**2-2.0d0*hp(1)*psp(i)-hp(2))**3
      d2psdhp(i,1,1)=2.0d0*psp(i)**3
     1        *(psp(i)**2-3.0d0*hp(1)*psp(i)-2.0d0*hp(2))/dummy
      d2psdhp(i,2,2)=-2.0d0*psp(i)*(hp(1)*psp(i)+hp(2))/dummy
      d2psdhp(i,3,3)=-8.0d0/9.0d0*(psp(i)-hp(1))/dummy
      d2psdhp(i,1,2)=psp(i)**2
     1        *(psp(i)**2-4.0d0*hp(1)*psp(i)-3.0d0*hp(2))/dummy
      d2psdhp(i,2,1)=d2psdhp(i,1,2)
      d2psdhp(i,2,3)=-2.0d0*(psp(i)**2+hp(2))/3.0d0/dummy
      d2psdhp(i,3,2)=d2psdhp(i,2,3)
      d2psdhp(i,3,1)=-4.0d0*psp(i)*(hp(1)*psp(i)+hp(2))/3.0d0/dummy
      d2psdhp(i,1,3)=d2psdhp(i,3,1)
c
      return
      end
c
c
c
c----------------------------------------------------------------(yld89)
c     akiyama YLD89
c
      subroutine jancae_yld89( s,se,dseds,d2seds2,nreq,
     &                         pryld,ndyld)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension s(3),dseds(3),d2seds2(3,3),pryld(ndyld)
      dimension s0(3),dsedsz(3),d2seds2z(3,3)
c
c
      call jancae_yld89_branch(s,se,dseds,d2seds2,nreq,
     &                         pryld,ndyld)
      if ( nreq.le.1 ) return
c
c                                         ---- Numerical differentiation
c                                            (mod. by H.Takizawa 150228)
      delta=1.0d-6*se
c
      do j=1,3
        s0(j)=s(j)
      enddo
c
      do i=1,3
c
        s0(i)=s(i)+delta
        call jancae_yld89_branch(s0,se0,dsedsz,d2seds2z,1,
     &                           pryld,ndyld)
        ta1=dsedsz(1)
        ta2=dsedsz(2)
        ta3=dsedsz(3)
c
        s0(i)=s(i)-delta
        call jancae_yld89_branch(s0,se0,dsedsz,d2seds2z,1,
     &                           pryld,ndyld)
        tb1=dsedsz(1)
        tb2=dsedsz(2)
        tb3=dsedsz(3)
c
        s0(i)=s(i)   ! re-set s0(i) for next loop
c
        d2seds2(1,i)=(ta1-tb1)/(2.0d0*delta)
        d2seds2(2,i)=(ta2-tb2)/(2.0d0*delta)
        d2seds2(3,i)=(ta3-tb3)/(2.0d0*delta)
      end do
c
c                                         ---- d2seds2 must be symmetric
      do i=1,3
        do j=i+1,3
          asy=0.5d0*(d2seds2(i,j)+d2seds2(j,i))
          d2seds2(i,j)=asy
          d2seds2(j,i)=asy
        enddo
      enddo
c
      return
      end
c
c
c----------------------------------------------------------------(yld89)
c
      subroutine jancae_yld89_branch(s,se,dseds,d2seds2,nreq,
     &                               pryld,ndyld)
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      dimension s(3),dseds(3),d2seds2(3,3),pryld(ndyld)
c
c     variable
c     f
c     pK1
c     pK2
c
c     First derivative variable
c     dsedf
      dimension dfdK(2),dKds(2,3)
c
c     Second derivative variable
c     d2sedf2
      dimension d2fdK2(2,2),dfdKdfdK(2,2),d2K1ds2(3,3),d2K2ds2(3,3)
c
c     f1,f2,f3
c     f=f1+f2+f3
      dimension df1dK(2),df2dK(2),df3dK(2)
c     dfdK=df1dK+df2dK+df3dK
      dimension d2f1dK2(2,2),d2f2dK2(2,2),d2f3dK2(2,2)
c     d2fdK2=d2f1dK2+d2f2dK2+d2f3dK2
      dimension u(2),v(2),w(2),x(2)
c      sc    //Scalar

c                                                        ---- parameters
      pM=pryld(1+1)
      a =pryld(1+2)
      h =pryld(1+3)
      p =pryld(1+4)
c
c                                                       ---- clear start
c
      call jancae_clear1 ( dfdK,2 )
      call jancae_clear2 ( dKds,2,3 )
c
      call jancae_clear2 ( d2fdK2,2,2 )
      call jancae_clear2 ( dfdKdfdK,2,2 )
      call jancae_clear2 ( d2K1ds2,3,3 )
      call jancae_clear2 ( d2K2ds2,3,3 )
c
      call jancae_clear1 ( df1dK,2 )
      call jancae_clear1 ( df2dK,2 )
      call jancae_clear1 ( df3dK,2 )
c
      call jancae_clear2 ( d2f1dK2,2,2 )
      call jancae_clear2 ( d2f2dK2,2,2 )
      call jancae_clear2 ( d2f3dK2,2,2 )
c                                                         ---- clear end
c
c     K1
c
      pK1=(s(1)+h*s(2))*0.5d0

c     K2
c     K22=K2^2
      pK3= (s(1)-h*s(2))*0.5d0
      pK22=pK3*pK3+p*p*s(3)*s(3)
      pK2=pK22**0.5d0
c
      f1=a*abs(pK1+pK2)**pM
      f2=a*abs(pK1-pK2)**pM
      f3=(2.0d0-a)*abs(2.0d0*pK2)**pM
c
      f=f1+f2+f3
c
      se=(0.5d0*f)**(1.0d0/pM)
      if ( nreq.eq.0 ) return
c
c
c                                                         ---- dKds(2,3)
c
      dKds(1,1)=0.5d0
      dKds(1,2)=h*0.5d0
      dKds(1,3)=0.0d0
c
      if(pK2.eq.0.0d0) then
      dKds(2,1)=0.0d0
      dKds(2,2)=0.0d0
      dKds(2,3)=0.0d0
      if(s(3).eq.0.0d0) dKds(2,3)=p*p
      else
      dKds(2,1)=pK3/(2.0d0*pK2)
      dKds(2,2)=-h*pK3/(2.0d0*pK2)
      dKds(2,3)=p*p*s(3)/pK2
      end if
c
c                                                  ---- d2K1ds2(3,3)=0.0
c
      call jancae_clear2 ( d2K1ds2,3,3 )
c
c                                                      ---- d2K2ds2(3,3)
c
c
      if(pK2.eq.0.0d0) then

      DpK22=1.0d-32
c
      d2K2ds2(1,1)=(DpK22**(-0.5d0)-pK3*DpK22**(-1.5d0))/4.0d0
      d2K2ds2(1,2)=-h*d2K2ds2(1,1)
      d2K2ds2(1,3)=-p*p*s(3)*pK3/2.0d0*DpK22**(-1.5d0)

      d2K2ds2(2,1)=d2K2ds2(1,2)
      d2K2ds2(2,2)=h*h*d2K2ds2(1,1)
      d2K2ds2(2,3)=-h*d2K2ds2(1,3)

      d2K2ds2(3,1)=d2K2ds2(1,3)
      d2K2ds2(3,2)=d2K2ds2(2,3)
      d2K2ds2(3,3)=p*p*(DpK22**(-0.5d0)-p*p*s(3)*s(3)*DpK22**(-1.5d0))
c
      else
c
      d2K2ds2(1,1)=(pK22**(-0.5d0)-pK3*pK22**(-1.5d0))/4.0d0
      d2K2ds2(1,2)=-h*d2K2ds2(1,1)
      d2K2ds2(1,3)=-p*p*s(3)*pK3/2.0d0*pK22**(-1.5d0)
c
      d2K2ds2(2,1)=d2K2ds2(1,2)
      d2K2ds2(2,2)=h*h*d2K2ds2(1,1)
      d2K2ds2(2,3)=-h*d2K2ds2(1,3)
c
      d2K2ds2(3,1)=d2K2ds2(1,3)
      d2K2ds2(3,2)=d2K2ds2(2,3)
      d2K2ds2(3,3)=p*p*(pK22**(-0.5d0)-p*p*s(3)*s(3)*pK22**(-1.5d0))

      end if
c
c
c                                                              ---- dfdK
c
c                                                             ---- df1dK
c
      do i=1,2
      df1dK(i)=a*pM*(pK1+pK2)*abs(pK1+pK2)**(pM-2.0d0)
      end do
c
c                                                             ---- df2dK
c
      df2dK(1)=a*pM*(pK1-pK2)*abs(pK1-pK2)**(pM-2.0d0)
      df2dK(2)=-a*pM*(pK1-pK2)*abs(pK1-pK2)**(pM-2.0d0)
c
c                                                             ---- df3dK
c
      df3dK(1)=0.0d0
      df3dK(2)=2.0d0*pM*(2.0d0-a)*(2.0d0*pK2)*abs(2.0d0*pK2)**(pM-2.0d0)
c
c                                            ---- dfdK=df1dK+df2dK+df3dK
c
      do i=1,2
        dfdK(i)=df1dK(i)+df2dK(i)+df3dK(i)
      end do
c
c                                                            ---- d2fdK2
c
c                                                           ---- d2f1dK2
c
      do i=1,2
        do j=1,2
      d2f1dK2(i,j)=a*pM*(pM-1.0d0)*abs(pK1+pK2)**(pM-2.0d0)
        end do
      end do
c
c                                                           ---- d2f2dK2
c
      d2f2dK2(1,1)=a*pM*(pM-1.0d0)*abs(pK1-pK2)**(pM-2.0d0)
      d2f2dK2(1,2)=-a*pM*(pM-1.0d0)*abs(pK1-pK2)**(pM-2.0d0)
      d2f2dK2(2,1)=d2f2dK2(1,2)
      d2f2dK2(2,2)=d2f2dK2(1,1)
c
c                                                           ---- d2f3dK2
c
      d2f3dK2(1,1)=0.0d0
      d2f3dK2(1,2)=0.0d0
      d2f3dK2(2,1)=0.0d0
      d2f3dK2(2,2)=4.0d0*pM*(pM-1.0d0)*(2.0d0-a)
     &                                 *abs(2.0d0*pK2)**(pM-2.0d0)
c
c                                    ---- d2fdK2=d2f1dK2+d2f2dK2+d2f3dK2
c
      do i=1,2
        do j=1,2
           d2fdK2(i,j)=d2f1dK2(i,j)+d2f2dK2(i,j)+d2f3dK2(i,j)
        end do
      end do
c
c                                                     ---- dfdKdfdK(2,2)
c
      do i=1,2
        do j=1,2
          dfdKdfdK(i,j)=dfdK(i)*dfdK(j)
        end do
      end do
c
c                                                             ---- dsedf
c
      dsedf=0.0d0
      dsedf=(0.5d0*f)**((1.0d0-pM)/pM)
      dsedf=dsedf/(2.0d0*pM)
c
c                                                           ---- d2sedf2
c
      d2sedf2=0.0d0
      d2sedf2=(0.5d0*f)**((1.0d0-2.0d0*pM)/pM)
      d2sedf2=d2sedf2*(1.0d0-pM)/(4.0d0*pM*pM)
c
c                                                             ---- dseds
      do i=1,3
        dseds(i)=0.0d0
        do j=1,2
          dseds(i)=dseds(i)+dfdK(j)*dKds(j,i)
        end do
        dseds(i)=dseds(i)*dsedf
      end do
      if ( nreq.eq.1 ) return
c
c
      if ( nreq.ge.2 ) return
c
c                                                           ---- d2seds2
c
      do i=1,3
        do j=1,3
c
          call jancae_clear1 (w,2)
          call jancae_clear1 (x,2)
c
                do k=1,2
                  do l=1,2
                w(k)=w(k)+dfdKdfdK(k,l)*dKds(l,j)
                x(k)=x(k)+d2fdK2(k,l)*dKds(l,j)
                  end do
                end do
c
                sc1=0.0d0
                sc2=0.0d0
                sc3=0.0d0
c
                do k=1,2
                  sc1=sc1+dKds(k,i)*w(k)
                  sc2=sc2+dKds(k,i)*x(k)
                end do
                  sc3=dfdK(1)*d2K1ds2(j,i)+dfdK(2)*d2K2ds2(j,i)
c
          d2seds2(i,j)=d2sedf2*sc1+dsedf*sc2+d2sedf2*sc3
c
        end do
      end do
c
c
      return
      end
c
c
c
c----------------------------------------------------------(yoshida2011)
c     F.Yoshida (2011,2013) yield function and its dfferentials
c
c     NUMISHEET 2011 Proceedings, AIP Conf. Proc.1383 (2011), pp.807-814
c
c     "A user-friendly 3D yield function to describe anisotropy of
c       steel sheets ",IJP,v.45(2013), pp.1119-139. )
c
      subroutine jancae_yoshida_2011 ( s,se,dseds,d2seds2,nreq,
     &                                 pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (maxa=100)
      dimension s(6),dseds(6),d2seds2(6,6),pryld(ndyld)
      dimension a(maxa),ipow(maxa,3)
c
c
      nd0=3
c
      n=0
      do it=0,nd0
        n=n+(nd0-it)*2+1
      enddo
      nterms=n
      if ( maxa.lt.nterms ) then
        write (6,*) 'increase maxa :',maxa,nterms
        call jancae_exit ( 9000 )
      endif
c
      n=0
      ipow=0
      do it=0,nd0
        ndmax=nd0*2-it*2
        do jy=0,ndmax
          jx=ndmax-jy
          n=n+1
          ipow(n,1)=jx
          ipow(n,2)=jy
          ipow(n,3)=it
        enddo
      enddo
c
      a     = 0.0
      a( 1)= 1.0d0        *pryld(1+ 1)    !       c1
      a( 2)=-3.0d0        *pryld(1+ 2)    !    -3*c2
      a( 3)= 6.0d0        *pryld(1+ 3)    !     6*c3
      a( 4)=-7.0d0        *pryld(1+ 4)    !    -7*c4
      a( 5)= 6.0d0        *pryld(1+ 5)    !     6*c5
      a( 6)=-3.0d0        *pryld(1+ 6)    !    -3*c6
      a( 7)= 1.0d0        *pryld(1+ 7)    !       c7
      a( 8)= 1.0d0* 9.0d0 *pryld(1+ 8)    !  9   *c8
      a( 9)=-2.0d0* 9.0d0 *pryld(1+ 9)    !  9*-2*c9
      a(10)= 3.0d0* 9.0d0 *pryld(1+10)    !  9* 3*c10
      a(11)=-2.0d0* 9.0d0 *pryld(1+11)    !  9*-2*c11
      a(12)= 1.0d0* 9.0d0 *pryld(1+12)    !  9   *c12
      a(13)= 1.0d0*27.0d0 *pryld(1+13)    ! 27   *c13
      a(14)=-1.0d0*27.0d0 *pryld(1+14)    ! 27*-1*c14
      a(15)= 1.0d0*27.0d0 *pryld(1+15)    ! 27   *c15
      a(16)= 1.0d0*27.0d0 *pryld(1+16)    ! 27   *c16
c
      call jancae_hy_polytype ( s,se,dseds,d2seds2,nreq,
     &                          nd0,a,ipow,maxa,nterms )
c
      return
      end
c
c
c
c----------------------------------------------------------(yoshida2011)
c     Weilong Hu & F.Yoshida style polynominal type yield function
c
      subroutine jancae_hy_polytype ( s,se,dseds,d2seds2,nreq,
     &                                nd0,a,ipow,maxa,nterms )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension s(6),dseds(6),d2seds2(6,6)
      dimension a(maxa),ipow(maxa,3)
      dimension sterm(3),ii(3),v(6)
c
c     nd        : order of polynominal nd=2*nd0
c     a(n)      : constants of function
c     ipow(n,i) : power of terms
c
      nd  =nd0*2
      dinv=1.0d0/float(nd)
c
      sterm(1)=s(1)-s(3)                ! sx-sz
      sterm(2)=s(2)-s(3)                ! sy-sz
      sterm(3)=s(4)**2+s(5)**2+s(6)**2  ! txy^2+tyz^2+tzx^2
c
      fai=0.0
      do n=1,nterms
        q=a(n)
        do k=1,3
          if ( ipow(n,k).gt.0 ) then
            q=q*sterm(k)**ipow(n,k)
          endif
        enddo
        fai=fai+q
      enddo
      se=fai**dinv
      if ( nreq.eq.0 ) return
c
      v=0.0
      do i=1,6
        idmax=1
        if ( i.eq.3 ) idmax=2
        do id=1,idmax
          do n=1,nterms
            do k=1,3
              ii(k)=ipow(n,k)
            enddo
            select case ( i )
            case ( 1,2 )
              dd=float(ii(i))
              ii(i)=ii(i)-1
            case ( 3 )
              dd=-1.0d0*float(ii(id))
              ii(id)=ii(id)-1
            case default
              dd=2.0d0*s(i)*float(ii(3))
              ii(3)=ii(3)-1
            end select
            q=dd*a(n)
            do k=1,3
              if ( ii(k).gt.0 ) then
                q=q*sterm(k)**ii(k)
              else if ( ii(k).lt.0 ) then
                q=0.0
              endif
            enddo
            v(i)=v(i)+q
          enddo
        enddo
      enddo
      ff=dinv*fai**(dinv-1.0d0)
      do i=1,6
        dseds(i)=ff*v(i)
      enddo
      if ( nreq.eq.1 ) return
c
      fff=dinv*(dinv-1.0d0)*fai**(dinv-2.0d0)
      do i=1,6
        do j=1,6
          d2seds2(i,j)=fff*v(i)*v(j)
        enddo
      enddo
      do i=1,6
        idmax=1
        if ( i.eq.3 ) idmax=2
        do id=1,idmax
          do j=1,6
            jdmax=1
            if (  j.eq.3               ) jdmax=2
            if ( (j.gt.3).and.(i.eq.j) ) jdmax=2
            do jd=1,jdmax
              do n=1,nterms
                do k=1,3
                  ii(k)=ipow(n,k)
                enddo
                select case ( i )
                case ( 1,2 )
                  ddi=float(ii(i))
                  ii(i)=ii(i)-1
                case ( 3 )
                  ddi=-1.0d0*float(ii(id))
                  ii(id)=ii(id)-1
                case default
                  ddi=2.0d0*s(i)*float(ii(3))
                  ii(3)=ii(3)-1
                end select
                select case ( j )
                case ( 1,2 )
                  ddj=float(ii(j))
                  ii(j)=ii(j)-1
                case ( 3 )
                  ddj=-1.0d0*float(ii(jd))
                  ii(jd)=ii(jd)-1
                case default
                  if ( jd.eq.1 ) then
                    ddj=2.0d0*s(j)*float(ii(3))
                    ii(3)=ii(3)-1
                  else
                    ddi=2.0d0*float(ipow(n,3))
                    ddj=1.0d0
                  endif
                end select
                q=a(n)*ddi*ddj
                do k=1,3
                  if ( ii(k).gt.0 ) then
                    q=q*sterm(k)**ii(k)
                  else if ( ii(k).lt.0 ) then
                    q=0.0
                  endif
                enddo
                d2seds2(i,j)=d2seds2(i,j)+ff*q
              enddo
            enddo
          enddo
        enddo
      enddo
      return
c
      end
c
c
c

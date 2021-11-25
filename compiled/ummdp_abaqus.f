************************************************************************
*                                                                      *
*                                  UMMDp                               *
*                                                                      *           
*                             <><><><><><><>                           *
*                                                                      *
*              UNIFIED MATERIAL MODEL DRIVER FOR PLASTICITY            *
*                                                                      *
*                         < PLUG-IN FOR ABAQUS >                       *
*                                                                      *
************************************************************************
*                                                                      *
*     > Copyright (c) 2018 JANCAE                                      *
*       . This software includes code originally developed by the      *
*       Material Modeling Working group of JANCAE.                     *
*
*     > Extended by M.G. Oliveira from University of Aveiro, Portugal  *
*       . Added additional isotropic hardening laws                    *
*       . Corrected Voigt notation for Yld2004-18p with Abaqus         *
*       . Linked kinematic hardening laws to the core of UMMDp         *
*       . Added Chaboche kinematic hardening law as used by Abaqus     *
*     	. Implemented uncoupled rupture criteria                       *
*       . Modified code to use only explicit variables                 *
*                                                                      *
************************************************************************
c
      SUBROUTINE UMAT ( STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     &    RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     &    TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     &    NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     &    DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC )
c
      INCLUDE 'ABA_PARAM.INC'
c
      CHARACTER*80 MATERL
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     &          DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     &          STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     &          PROPS(NPROPS),COORDS(3),DROT(3,3),
     &          DFGRD0(3,3),DFGRD1(3,3)
c
c***********************************************************************
c-----------------------------------------------------------------------
      common /jancae1/ne,ip,lay
      common /jancae3/prop
      common /jancaea/nsdv
      common /jancaeb/propdim
c
      parameter (mxpbs=10)
c
      dimension s2(ntens),dpe(ntens),x1(mxpbs,ntens),x2(mxpbs,ntens),
     &          pe(ntens)
c
      dimension ustatev(6)
c
      parameter (mxprop=100)
      dimension prop(mxprop)
c-----------------------------------------------------------------------
c
c                        ne  : element no.
c                        ip  : integration point no.
c                        lay : layer no. of shell
      ne = noel
      ip = npt
      lay = kspt
      if ( lay == 0 ) lay = 1
      nsdv = nstatv
      nprop = mxprop
      propdim = nprops - 1
c
c                                        ---- set debug and verbose mode
      nvbs0 = props(1)
      call ummdp_debugmode ( nvbs,nvbs0 )
c                                       ---- output detailed information
      if ( nvbs >= 4 ) then
        call ummdp_print_info  ( kinc,ndi,nshr )
        call ummdp_print_inout ( 0,stress,dstran,ddsdde,ntens,statev,
     1                           nstatv )
      end if
c
c                                           ---- set material properties
      do i = 2,nprops
        prop(i-1) = props(i)
      end do
c
      call ummdp_prop_dim ( prop,nprop,propdim,ndela,ndyld,ndihd,ndkin,
     1                      npbs,ndrup )                       
      if ( npbs > mxpbs ) then
        write (6,*) 'npbs > mxpbs error in umat'
        write (6,*) 'npbs =',npbs
        write (6,*) 'mxpbs=',mxpbs
        call ummdp_exit ( 9000 )
      end if
c                                                      ---- check nstatv
      call ummdp_check_nisv ( nstatv,ntens,npbs )
c                             ---- copy current internal state variables
      call ummdp_isvprof ( isvrsvd,isvsclr )
      call ummdp_isv2pex ( isvrsvd,isvsclr,statev,nstatv,p,pe,x1,ntens,
     1                     mxpbs,npbs )                    
c
c                             ---- update stress and set tangent modulus
      mjac = 1
      call ummdp_plasticity ( stress,s2,dstran,p,dp,dpe,de33,x1,x2,
     1                        mxpbs,ddsdde,ndi,nshr,ntens,nvbs,mjac,
     2                        prop,nprop,propdim )                      
c                                                     ---- update stress
      do i = 1,ntens
        stress(i) = s2(i)
      end do
c                                            ---- update eq.plast,strain
      statev(isvrsvd+1) = p + dp
c                                         ---- update plast.strain comp.
      call rotsig ( statev(isvrsvd+2),drot,ustatev,2,ndi,nshr )
c
      do i = 1,ntens
        is = isvrsvd + isvsclr + i
        statev(is) = ustatev(i) + dpe(i)
      end do
c                                       ---- update of back stress comp.
      if ( npbs /= 0 ) then
        do n = 1,npbs
          do i = 1,ntens
            is = isvrsvd + isvsclr + ntens*n + i
            statev(is) = x2(n,i)
          end do
        end do
      end if
c                           ----  if debug mode, output return arguments
      if ( nvbs >= 4 ) then
        call ummdp_print_inout ( 1,stress,dstran,ddsdde,ntens,statev,
     1                           nstatv )
      end if
c
      return
      end subroutine umat
c
c
c
c***********************************************************************
      SUBROUTINE SDVINI ( STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,
     1                    LAYER,KSPT )
c-----------------------------------------------------------------------
      INCLUDE 'ABA_PARAM.INC'
c
      DIMENSION STATEV(NSTATV),COORDS(NCRDS)
c-----------------------------------------------------------------------
c
      ne = noel
      ip = npt
      lay = kspt
      if ( lay == 0 ) lay = 1
c
      if ( ne*ip*lay == 1 ) then
        write (6,*) 'SDVINI is called. '
      end if
c
      do n = 1,nstatv
        statev(n) = 0.0
      end do
c
      return
      end subroutine sdvini
c
c
c
c***********************************************************************
      SUBROUTINE UVARM ( UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     1    NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     2    JMAC,JMATYP,MATLAYO,LACCFLA )
c-----------------------------------------------------------------------
      INCLUDE 'ABA_PARAM.INC'
c
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  flgray(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION array(15),jarray(15),JMAC(*),JMATYP(*),COORD(*)
c-----------------------------------------------------------------------
c     The dimensions of the variables flgray, array and jarray
c     must be set equal to or greater than 15.
c
      parameter (maxsdv=50)
      parameter (mxpbs=10)
      parameter (mxprop=100)
c
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
      dimension   sdv(maxsdv),uvar1(nuvarm)
c
      dimension prop(mxprop)
      real*8,allocatable,dimension(:) :: prela,pryld,prihd,prkin,prrup
c-----------------------------------------------------------------------
c
      nprop = mxprop
c
c     variables list :
c        uvar codes and state variables arrays
c        nt : ntens
c
c     statev(1                    ) : equivalent plastic strain
c     statev(2        ~ 1+nt      ) : plastic strain components
c     statev(1+nt*i+1 ~ 1+nt*(i+1)) : partial back stress component xi
c
c     uvar(1     ) : equivalent stress
c     uvar(2     ) : flow stress
c     uvar(3~2+nt) : total back stress components
c     uvar(3~2+nt) : rupture criterion
c
      ne = noel
      ip = npt
      lay = kspt
      if ( lay == 0 ) lay = 1
      ntens = ndi + nshr

c                                            ---- get uvar before update
      do i = 1,nuvarm
        uvar1(i) = uvar(i)
      end do
c                                                        ---- get stress
      call getvrm ( 'S',array,jarray,flgray,jrcd,jmac,jmatyp,matlayo,
     1              laccfla )
      if ( jrcd /= 0 ) then
        write (6,*) 'request error in uvarm for s'
        write (6,*) 'stop in uvrm.'
        call ummdp_exit ( 9000 )
      end if
c
      do i = 1,ndi
        s(i) = array(i)
      end do
      do i = 1,nshr
        i1 = ndi + i
        i2 = 3 + i
        s(i1) = array(i2)
      end do
c                                               ---- get state variables
      if ( nsdv > maxsdv ) then
        write (6,*) 'increase dimension of ARRAY2 and JARRAY2'
        write (6,*) 'stop in uvrm.'
        call ummdp_exit ( 9000 )
      end if
      call getvrm ( 'SDV',ARRAY2,JARRAY2,FLGRAY2,jrcd,jmac,jmatyp,
     1              matlayo,laccfla)
      if ( jrcd /= 0 ) then
        write (6,*) 'request error in uvarm for sdv'
        write (6,*) 'stop in uvrm.'
        call ummdp_exit ( 9000 )
      end if
      do i = 1,nsdv
        sdv(i) = array2(i)
      end do
c                                           ---- set material properties
      call ummdp_prop_dim ( prop,nprop,propdim,ndela,ndyld,ndihd,ndkin,
     1                      npbs,ndrup )                        
      allocate( prela(ndela) )
      allocate( pryld(ndyld) )
      allocate( prihd(ndihd) )
      allocate( prkin(ndkin) )
      allocate( prrup(ndrup) )
      k = 0
      do i = 1,ndela
        k = k + 1
        prela(i) = prop(k)
      end do
      do i = 1,ndyld
        k = k + 1
        pryld(i) = prop(k)
      end do
      do i = 1,ndihd
        k = k + 1
        prihd(i) = prop(k)
      end do
      do i = 1,ndkin
        k = k + 1
        prkin(i) = prop(k)
      end do
      do i = 1,ndrup
        k = k + 1
        prrup(i) = prop(k)
      end do
c
c                                                  ---- calc back stress
      call ummdp_isvprof ( isvrsvd,isvsclr )
      call ummdp_isv2pex ( isvrsvd,isvsclr,sdv,maxsdv,p,pe,x,ntens,
     1                     mxpbs,npbs )                      
      do i = 1,ntens
         xsum(i) = 0.0
      end do
      if ( npbs /= 0 ) then
        do i = 1,ntens
          do nb = 1,npbs
            xsum(i) = xsum(i) + x(nb,i)
          end do
        end do
      end if
c                                                 ---- equivalent stress
      if ( nuvarm >= 1 ) then
        do i = 1,ntens
          eta(i) = s(i) - xsum(i)
        end do
        call ummdp_yield ( se,dseds,d2seds2,0,eta,ntens,ndi,nshr,pryld,
     1                     ndyld )                        
        uvar(1) = se
      end if
c                                                       ---- flow stress
      if ( nuvarm >= 2 ) then
        call ummdp_isotropic ( sy,dsydp,d2sydp2,0,p,prihd,ndihd )
        uvar(2) = sy
      end if
c                                                       ---- back stress
      if ( npbs /= 0 ) then
        if ( nuvarm >= 3 ) then
          do i = 1,ntens
            uvar(2+i) = xsum(i)
          end do
        end if
      end if
c                                                 ---- rupture criterion
      if ( prrup(1) /= 0) then
        nt = ntens
        if ( npbs == 0 ) nt = 0
        if ( nuvarm >= (3+nt) ) then
          call ummdp_rupture ( ntens,sdv,nsdv,uvar,uvar1,nuvarm,jrcd,
     1                         jmac,jmatyp,matlayo,laccfla,nt,ndrup,
     2                         prrup) 
        end if
      end if
c
      return
      end subroutine uvarm
c
c
c
************************************************************************
c     SET INTERNAL STATE VARIABLES PROFILE
c
      subroutine ummdp_isvprof ( isvrsvd,isvsclr )
c-----------------------------------------------------------------------
      INCLUDE 'ABA_PARAM.INC'
c
      isvrsvd = 0             ! no reserved variables
c
      isvsclr = 1             ! statev(1) is for eq.plast.strain
c
      return
      end subroutine ummdp_isvprof
c
c
c
************************************************************************
c     EXIT PROGRAM BY ERROR
c
      subroutine ummdp_exit (nexit)
c-----------------------------------------------------------------------
      INCLUDE 'ABA_PARAM.INC'
      common /jancae1/ne,ip,lay
c-----------------------------------------------------------------------
      write (6,*) 'error code :',nexit
      write (6,*) 'element no.           :',ne
      write (6,*) 'integration point no. :',ip
      write (6,*) 'layer no.             :',lay
c
      call xit
c
      return
      end subroutine ummdp_exit
c
c
cc
c                               OVER THIS LINE IS CODE DEPENDENT
c<<-->><<-->><<-->><<->><<-->><<-->><<-->><<-->><<-->><<-->><<-->><<-->>
c                            UNDER THIS LINE IS DOCE INDEPENDENT
c
c     UMMDp: UNIFIED MATERIAL MODEL DRIVER FOR PLASTICITY
c
c
************************************************************************
c     PLASTICITY DUMMY
c
      subroutine ummdp_plasticity ( s1,s2,de,p,dp,dpe,de33,x1,x2,mxpbs,
     1                              ddsdde,nnrm,nshr,nttl,nvbs,mjac,
     2                              prop,nprop,propdim )                               
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
      character*32 text
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
          write (6,*) 'prop(',i,')=',prop(i)
        end do
        call ummdp_exit ( 9000 )
      end if
      if ( nvbs >= 4 ) then
        do i = 1,n
          write (6,*) 'prop(',i,')=',prop(i)
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
      common /jancae1/ne,ip,lay
      common /jancae2/n1234
      integer ne,ip,lay
c
      integer,intent(in) :: mxpbs,nnrm,nshr,nttl,nvbs,mjac,nprop,npbs,
     1                      ndela,ndyld,ndihd,ndkin,ndrup,nnn
      real*8 ,intent(in) :: p
      real*8 ,intent(in) :: s1(nttl),de(nttl),prop(nprop)
c
      real*8,intent(out) :: de33,dp
      real*8,intent(out) :: s2(nttl),dpe(nttl),ddsdde(nttl,nttl)
     1                      
      real*8 ,intent(inout) :: x1(mxpbs,nttl),x2(mxpbs,nttl)
c
      integer i,j,k,n,m,maxnr,ndiv,maxnest,nout,n1234,i1,i2,j1,j2,k1,k2,
     1        nest,newmstg,nite,nstg,mstg,knr,ip1,ip2,nsym
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
      logical debug
      character*32 text
c
c-----------------------------------------------------------------------
c     >>> Arguments List
c
c       ne      | index number of element                           (in)
c       ip      | index number of integration point                 (in)
c       lay     | index number of layer (shell and membrane)        (in)
c
c       nnrm    | number of normal components                       (in)
c       nshr    | number of shear  components                       (in)
c       nttl    | total number of components                        (in)
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
c       de33    | strain increment in thickness direction          (out)  
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
c       dseds   | 1st order derivative of equivalent stress wrt stress            
c       d2seds2 | 2nd order derivative of equivalent stress wrt stress
c       stry    | trial stress predicted elastically
c       sy      | flow stress (function of equivalent plastic strain)
c       dsydp   | 1st order derivative of flow stress wrt 
c                  equivalent plastic strain
c       g1      | error of stress point to yield surface
c       g2      | error of direction of plastic strain increment to
c                  normal of yield surface (error vector)
c       g2n     | norm of g2 vector
c       eta     | stress for yield function {s}-{xt}
c       xt1     | total back stress before update
c       xt2     | total back stress after update
c       vk      | evolution for kinematic hardening dx=dp*vk
c       dvkdp   | derivative of v wrt equivalent plastic strain
c       dvkds   | derivative of v wrt stress
c       dvkdx   | derivative of v wrt partial back stress
c       dvkdxt  | derivative of v wrt total back stress
c       g3      | error of evolution for back stress (error vector)
c       g3n     | norm of g3 vectors
c       sgap    | stress gap to be eliminated in multistage steps
c       tol     | covergence tolerance
c       maxnr   | maximum iterations of Newton-Raphson
c       maxnest | maximum trial times of multistage gap reduction
c       ndiv    | division number of multistage
c
c-----------------------------------------------------------------------
c
      debug = .true.
      debug = .false.
      tol = 1.0d-5
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
      if ( (nvbs >= 1) .or. (nout /= 0) ) then
        write (6,*)
        write (6,*) '******************************'
        write (6,*) '******* START OF UMMDp *******'
        write (6,*) '******************************'
      end if
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
        write (6,*) 'MATERIAL DATA LIST'
        call ummdp_print_elastic   ( prela,ndela )
        call ummdp_print_yield     ( pryld,ndyld )
        call ummdp_print_isotropic ( prihd,ndihd )
        call ummdp_print_kinematic ( prkin,ndkin,npbs )
        call ummdp_print_rupture   ( prrup,ndrup )
      end if
c                                                           ---- set [U]
      call ummdp_utility_clear2( um,nttl,nnn )
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
            x2(n,i) = 0.0d0
            x1(n,i) = 0.0d0
          end do
        end do
      end if
c
      de33 = 0.0d0
      dp = 0.0d0
      do i = 1,nttl
        dpe(i) = 0.0d0
      end do
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
        text = 'current stress (input)'
        call ummdp_utility_print1 ( text,s1,nttl )
        text = 'strain inc. (input)'
        call ummdp_utility_print1 ( text,de,nttl )
        if ( npbs /= 0 ) then
          text = 'part. back stess (input)'
          call ummdp_backprint ( text,npbs,x1,nttl,mxpbs )
          text = 'total  back stess (input)'
          call ummdp_utility_print1 ( text,xt1,nttl )
        end if
      end if
c                                            ---- set elastic [D] matrix
      call ummdp_setdelast ( delast,prela,ndela,nttl,nnrm,nshr,d33d )                        
c                                             ---- copy delast to ddsdde
      do i = 1,nttl
        do j = 1,nttl
          ddsdde(i,j) = delast(i,j)
        end do
      end do
      if ( nvbs >= 5 ) then
        text = 'elastic matrix'
        call ummdp_utility_print2 ( text,ddsdde,nttl,nttl )
      end if
c                                                ---- elastic prediction
      call ummdp_utility_mv ( vv,ddsdde,de,nttl,nttl )
      do i = 1,nttl
        s2(i) = s1(i) + vv(i)
      end do
      if ( nvbs >= 5 ) then
        text = 'elastic predicted stress'
        call ummdp_utility_print1 ( text,s2,nttl )
      end if
c                                                       ---- back stress
      do i = 1,nttl
        eta(i) = s2(i) - xt2(i)
      end do
c                                                       ---- check yield
      call ummdp_yield  ( se,dseds,d2seds2,0,eta,nttl,nnrm,nshr,pryld,
     1                    ndyld )                      
      call ummdp_isotropic ( sy,dsydp,d2sydp2,0,p,prihd,ndihd )                        
c
      if ( nvbs >= 3 ) then
        write (6,*) 'plastic strain p=',p
        write (6,*) 'flow stress   sy=',sy
        write (6,*) 'equiv.stress  se=',se
        if ( npbs /= 0 ) then
          call ummdp_yield  ( xe,dseds,d2seds2,0,xt1,nttl,nnrm,nshr,
     1                        pryld,ndyld )                          
          write (6,*) 'equiv.back.s  xe=',xe
        end if
      end if
      if ( se <= sy ) then
        if ( nvbs >= 3 ) write (6,*) 'judge : elastic'
        if ( (nttl == 3) .or. (nttl == 5) ) then
          de33 = 0.0d0
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
        call ummdp_yield ( se,dseds,d2seds2,2,eta,nttl,nnrm,nshr,pryld,
     1                     ndyld )                       
c
        if ( nvbs >= 5 ) then
          text = 's2'
          call ummdp_utility_print1 ( text,s2,nttl )
          if ( npbs /= 0 ) then
            text = 'xt2'
            call ummdp_utility_print1 ( text,xt2,nttl )
            text = 'eta'
            call ummdp_utility_print1 ( text,eta,nttl )
          end if
          text = 'dse/ds'
          call ummdp_utility_print1 ( text,dseds,nttl )
          text = 'd2se/ds2'
          call ummdp_utility_print2 ( text,d2seds2,nttl,nttl )
        end if
c                                        ---- calc. sy and differentials
        call ummdp_isotropic ( sy,dsydp,d2sydp2,1,pt,prihd,ndihd )                        
c
        if ( nvbs >= 5 ) then
          write (6,*) 'plastic strain p=',pt
          write (6,*) 'flow stress   sy=',sy
          write (6,*) 'hardening dsy/dp=',dsydp
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
          write (6,*) 'g1 (yield surf) =',g1
          write (6,*) 'g2n (normality) =',g2n
          if ( nvbs >= 5 ) then
            text = 'g2 vector'
            call ummdp_utility_print1 ( text,g2,nttl )
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
                  call ummdp_utility_print1 ( text,uv,nttl )
                end if
              end do
            end if
          end if
        end if
c                      ---- calc. dependencies common for NR and Dds/Dde
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
        call ummdp_utility_clear1( wv,nnn )
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
        call ummdp_utility_minv ( ami,am,nnn,det )
c                                                     ---- [C]=[U][A]^-1
        call ummdp_utility_mm ( cm,um,ami,nttl,nnn,nnn )
c
c
c                                                 ---- check convergence
        if ( (abs(g1  /sy) <= tol) .and.
     1       (abs(g2n /sy) <= tol) .and.
     2       (abs(g3nn/sy) <= tol) ) then
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
        call ummdp_utility_mv ( vv,cm,gv,nttl,nnn )
        call ummdp_utility_vvs ( top0,dseds,vv,nttl )
        top = g1 - top0
        call ummdp_utility_mv ( vv,cm,wv,nttl,nnn )
        call ummdp_utility_vvs ( bot0,dseds,vv,nttl )
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
          call ummdp_utility_clear1( vv,nttl )
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
          text = 'current stress (input)'
          call ummdp_utility_print1 ( text,s1,nttl )
          text = 'strain inc. (input)'
          call ummdp_utility_print1 ( text,de,nttl )
          write (6,*) 'eq.plast.strain (input)'
          write (6,*) p
          write (6,*) 'the proposals to fix this error'
          write (6,*) ' reduce the amount of strain per inc.'
          write (6,*) ' increase maxnest in program',maxnest
          write (6,*) ' increase ndiv    in program',ndiv
          write (6,*) ' increase maxnr   in program',maxnr
          write (6,*) ' increase tol     in program',tol
          call ummdp_exit ( 9000 )
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
        call ummdp_utility_print1 ( text,s2,nttl )
        text = 'plastic strain inc'
        call ummdp_utility_print1 ( text,dpe,nttl )
        if ( npbs /= 0 ) then
          text = 'updated part. back stess'
          call ummdp_backprint ( text,npbs,x2,nttl,mxpbs )
          text = 'updated total back stess'
          call ummdp_utility_print1 ( text,xt2,nttl )
        end if
      end if
c                                    ---- calc. strain inc. in thickness
      if ( (nttl == 3) .or. (nttl == 5) ) then
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
            ddsdde(i,j) = 0.0d0
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
      call ummdp_utility_clear2( bm,nnn,nttl )
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
      call ummdp_utility_clear2 ( em2,nttl,nnn )
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
          call ummdp_utility_print2 ( text,ddsdde,nttl,nttl )
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
        call ummdp_utility_print2 ( text,ddsdde,nttl,nttl )
      end if
c
      return
      end subroutine ummdp_plasticity_core
c
c
c
c***********************************************************************
c     SET DEBUG AND PRINT VERBOSE MODE
c
      subroutine ummdp_debugmode ( nvbs,nvbs0 )
c-----------------------------------------------------------------------
      implicit none
c
      common /jancae1/ne,ip,lay
      integer ne,ip,lay
c
      integer,intent(in) :: nvbs0
c
      integer,intent(out) :: nvbs
c
      integer nechk,ipchk,laychk,nchk
c-----------------------------------------------------------------------
c                             specify verbose level and point
c      nvbs0 = 0   ! verbose mode
c
c           0  error message only
c           1  summary of MsRM
c           2  detail of MsRM and summary of NR
c           3  detail of NR
c           4  input/output
c           5  all status for debug
c
c       MsRM : Multistage Return Mapping
c       NR   : Newton-Raphson
c-----------------------------------------------------------------------
      nechk = 1     ! element no. to be checked
      ipchk = 1     ! integration point no. to checked
      laychk = 1    ! layer no. to be checked
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
c     SET ELASTIC MATERIAL JACOBIAN MATRIX
c
      subroutine ummdp_setdelast ( delast,prela,ndela,nttl,nnrm,nshr,
     1                             d33d )
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
        call ummdp_utility_clear2( delast3d,6,6 )
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
        write (6,*) 'elasticity code error in ummdp_setelast'
        write (6,*) 'ntela=',ntela
        call ummdp_exit ( 9000 )
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
     1                      - delast3d(i3,3)*delast3d(3,j3)/d33
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
      end subroutine ummdp_setdelast
c
c
c
************************************************************************
c     CHECK DIMENSIONS OF INTERNAL STATE VARIABLES
c
      subroutine ummdp_check_nisv ( nisv,nttl,npbs )
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
        call ummdp_exit ( 9000 )
      end if
c
      return
      end subroutine ummdp_check_nisv
c
c
c
************************************************************************
c     SET VARIABLES FROM STATE VARIABLES
c
      subroutine ummdp_isv2pex ( isvrsvd,isvsclr,stv,nstv,p,pe,x,nttl,
     1                           mxpbs,npbs )                           
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
c                                                 ---- eq.plastic strain
      p = stv(isvrsvd+1)
c                                              ---- plastic strain comp.
      do i = 1,nttl
        pe(i) = stv(isvrsvd + isvsclr + i)
      end do
c                                         ---- partial back stress comp.
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
c     SUM PARTIAL BACK STRESS FOR TOTAL BACK STREESS
c
      subroutine ummdp_backsum ( npbs,xt,x,nttl,mxpbs )
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
      integer,intent(in) :: npbs,nttl,mxpbs
      real*8 ,intent(in) :: x(mxpbs,nttl)
      character*32,intent(in) :: text
c
      integer i,j
      real*8 xx(npbs,nttl)
c-----------------------------------------------------------------------
c
      if ( npbs == 0 ) return
c
      do i = 1,nttl
        do j = 1,npbs
          xx(j,i) = x(j,i)
        end do
      end do
      call ummdp_utility_print2 ( text,xx,npbs,nttl )
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
c
      case (0) ; nd = 2
        case (1) ; nd = 2
c
        case default
          write (6,*) 'error elastic property id :',nela
          call ummdp_exit ( 9000 )
c
      end select
      ndela = nd + 1
c
      n = ndela
      nyld = nint(prop(n+1))
      select case (nyld)
c
        case ( 0 ) ! von Mises
          nd = 0                           
        case ( 1 ) ! Hill 1948
          nd = 6                           
        case ( 2 ) ! Yld2004-18p
          nd = 19                          
        case ( 3 ) ! CPB2005
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
        case ( -4 ) ! BBC2005
          nd = 9                          
        case ( -5 ) ! Yld89
          nd = 4                          
        case ( -6 ) ! BBC2008
          nd = 2 + 8*nint(prop(n+2))      
        case ( -7 ) ! Hill1990
          nd = 0.5d0                      
c
        case default
          write (6,*) 'error yield function id :',nyld
          call ummdp_exit ( 9000 )
c
      end select
      ndyld = nd + 1
c
      n = ndela + ndyld
      nihd = nint(prop(n+1))
      select case ( nihd )
c
        case ( 0 ) ! Perfecty Plastic
          nd = 1
        case ( 1 ) ! Linear
          nd = 2
        case ( 2 ) ! Swift
          nd = 3 
        case ( 3 ) ! Ludwick
          nd = 3
        case ( 4 ) ! Voce
          nd = 3
        case ( 5 ) ! Voce + Linear
          nd = 4
        case ( 6 ) ! Voce + Swift
          nd = 7
c
        case default
          write (6,*) 'error work hardening curve id :',nihd
          call ummdp_exit ( 9000 )
c
      end select
      ndihd = nd + 1
c
      n = ndela + ndyld + ndihd
      nkin = nint(prop(n+1))
      select case ( nkin )
c
        case ( 0 ) ! No Kinematic Hardening
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
        case ( 4 ) ! Chaboche
          npbs = nint(prop(n+2))
          nd = 2*npbs + 1
        case ( 5 ) ! Chaboche - Ziegler
          npbs = nint(prop(n+2))
          nd = 2*npbs + 1
        case ( 6 ) ! Yoshida-Uemori
          nd = 5
          npbs = 2      
c
        case default
          write (6,*) 'error kinematic hardening id :',nkin
          call ummdp_exit ( 9000 )
c
      end select
      ndkin = nd + 1
c
      n = ndela + ndyld + ndihd + ndkin
      nrup = nint(prop(n+1))
      select case (nrup)
c
        case ( 0 ) ! No Uncoupled Rupture Criterion
          nd = 0             
        case ( 1 ) ! Equivalent Plastic Strain
          nd = 2
        case ( 2 ) ! Cockroft and Latham
          nd = 2
        case ( 3 ) ! Rice and Tracey
          nd = 2
        case ( 4 ) ! Ayada
          nd = 2
        case ( 5 ) ! Brozzo
          nd = 2
        case ( 6 ) ! Forming Limit Diagram
          nd = 2
c
        case default
          write (6,*) 'error uncoupled rupture criterion id :',nrup
          call ummdp_exit ( 9000 )
c
      end select
      ndrup = nd + 1
c
      return
      end subroutine ummdp_prop_dim
c
c
c
************************************************************************
*
*     ISOTROPIC HARDENING LAWS
*
************************************************************************
c
c      0 : Perfectly Plastic
c      1 : Linear
c      2 : Swift
c      3 : Ludwick
c      4 : Voce
c      5 : Voce + Linear
c      6 : Voce + Swift
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     CALCULATE ISOTROPIC HARDENING LAW 
c
      subroutine ummdp_isotropic ( sy,dsydp,d2sydp2,nreq,p,prihd,ndihd )              
c-----------------------------------------------------------------------
      implicit none
c
      integer ndihd,nreq
      real*8 sy,dsydp,d2sydp2,p
      real*8 prihd(ndihd)
c
      integer ntihd
      real*8 sy0,hard,c,e0,en,q,b,a
c-----------------------------------------------------------------------
c
      ntihd = nint(prihd(1))
      select case ( ntihd )
c
      case ( 0 )                                     ! Perfectly Plastic
        sy = prihd(1+1)
        if ( nreq . ge.1 ) then
          dsydp = 0.0
          if ( nreq >= 2 ) then
            d2sydp2 = 0.0
          end if
        end if
c
      case ( 1 )                                                ! Linear
        sy0  = prihd(1+1)
        hard = prihd(1+2)
        sy = sy0 + hard*p
        if ( nreq >= 1 ) then
          dsydp = hard
          if ( nreq >= 2 ) then
            d2sydp2 = 0.0
          end if
        end if
c
      case ( 2 )                                                 ! Swift
        c  = prihd(1+1)
        e0 = prihd(1+2)
        en = prihd(1+3)
        sy = c*(e0+p)**en
        if ( nreq >= 1 ) then
          dsydp = en*c*(e0+p)**(en-1.0d0)
          if ( nreq >= 2 ) then
            d2sydp2 = en*c*(en-1.0d0)*(e0+p)**(en-2.0d0)
          end if
        end if
c
      case ( 3 )                                               ! Ludwick
        sy0 = prihd(1+1)
        c   = prihd(1+2)
        en  = prihd(1+3)
        sy = sy0+c*p**en
        if ( nreq >= 1 ) then
          dsydp = en*c*p**(en-1.0d0)
          if ( nreq >= 2 ) then
            d2sydp2 = en*c*(en-1.0d0)*p**(en-2.0d0)
          end if
        end if
c
      case ( 4 )                                                  ! Voce
        sy0 = prihd(1+1)
        q   = prihd(1+2)
        b   = prihd(1+3)
        sy = sy0+q*(1.0d0-exp(-b*p))
        if ( nreq >= 1 ) then
          dsydp = q*b*exp(-b*p)
          if ( nreq >= 2 ) then
            d2sydp2 = -q*b*b*exp(-b*p)
          end if
        end if
c
      case ( 5 )                                         ! Voce + Linear
        sy0 = prihd(1+1)
        q   = prihd(1+2)
        b   = prihd(1+3)
        c   = prihd(1+4)
        sy = sy0+q*(1.0d0-exp(-b*p))+c*p
        if ( nreq >= 1 ) then
          dsydp = q*b*exp(-b*p)+c
          if ( nreq >= 2 ) then
            d2sydp2 = -q*b*b*exp(-b*p)
          end if
        end if
c
      case ( 6 )                                          ! Voce + Swift
        a   = prihd(1+1)
        sy0 = prihd(1+2)
        q   = prihd(1+3)
        b   = prihd(1+4)
        c   = prihd(1+5)
        e0  = prihd(1+6)
        en  = prihd(1+7)
        sy = a*(sy0+q*(1.0d0-exp(-b*p))) + (1.0d0-a)*(c*(e0+p)**en)
        if ( nreq >= 1 ) then
          dsydp = a*(q*b*exp(-b*p)) +(1.0d0-a)*(en*c*(e0+p)**(en-1.0d0))
          if ( nreq >= 2 ) then
            d2sydp2 = a*(-q*b*b*exp(-b*p)) +
     &                (1.0d0-a)*(en*c*(en-1.0d0)*(e0+p)**(en-2.0d0))
          end if
        end if
c
      case default
        write (6,*) 'hardening type error',ntihd
        call ummdp_exit (9000)
      end select
c
      return
      end subroutine ummdp_isotropic
c
c
c
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
c     CALCULATE KINEMATIC HARDENING LAW
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
      case ( 3 )                          ! Armstrong & Frederick (1966)
        call ummdp_kinematic_armstrong ( vk,dvkdp,dvkds,dvkdx,dvkdxt,p,
     1                                   s,x,xt,nttl,nnrm,nshr,mxpbs,
     2                                   npbs,prkin,ndkin,pryld,ndyld )
c
      case ( 4 )                                       ! Chaboche (1979)
        call ummdp_kinematic_chaboche ( vk,dvkdp,dvkds,dvkdx,dvkdxt,p,
     1                                  s,x,xt, nttl,nnrm,nshr,mxpbs,
     2                                  npbs,prkin,ndkin,pryld,ndyld )
c
      case ( 5 )                             ! Chaboche (1979) - Ziegler
        call ummdp_kinematic_chaboche_ziegler ( vk,dvkdp,dvkds,dvkdx,
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
        call ummdp_exit ( 9000 )
      end select
c
      return
      end subroutine ummdp_kinematic
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CALCULATE 1ST AND 2ND DERIVATIVES FOR KINEMATIC HARDENING LAWS
c
      subroutine ummdp_kinematic_dseds ( eta,seta,dseds,d2seds2,nttl,
     1                                   nnrm,nshr,pryld,ndyld )
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
c     PRAGER KINEMATIC HARDENING LAW
c
      subroutine ummdp_kinematic_prager ( vk,dvkdp,dvkds,dvkdx,dvkdxt,
     1                                    p,s,x,xt,nttl,nnrm,nshr,mxpbs,
     2                                    npbs,prkin,ndkin,pryld,ndyld )
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
c     ZIEGLER KINEMATIC HARDENING LAW
c
      subroutine ummdp_kinematic_ziegler ( vk,dvkdp,dvkds,dvkdx,dvkdxt,
     1                                     p,s,x,xt,nttl,nnrm,nshr,
     2                                     mxpbs,npbs,prkin,ndkin,
     3                                     pryld,ndyld )
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
c     ARMSTRONG-FREDERICK (1966) KINEMATIC HARDENING LAW
c
      subroutine ummdp_kinematic_armstrong ( vk,dvkdp,dvkds,dvkdx,
     1                                       dvkdxt,p,s,x,xt,nttl,nnrm,
     2                                       nshr,mxpbs,npbs,prkin,
     3                                       ndkin,pryld,ndyld )
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
c     CHABOCHE (1979) KINEMATIC HARDENING LAW
c
      subroutine ummdp_kinematic_chaboche ( vk,dvkdp,dvkds,dvkdx,dvkdxt,
     1                                      p,s,x,xt,nttl,nnrm,nshr,
     2                                      mxpbs,npbs,prkin,ndkin,
     3                                      pryld,ndyld )
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
      end subroutine ummdp_kinematic_chaboche
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CHABOCHE (1979) - ZIEGLER KINEMATIC HARDENING LAW
c
      subroutine ummdp_kinematic_chaboche_ziegler ( vk,dvkdp,dvkds,
     1                                              dvkdx,dvkdxt,p,s,x,
     2                                              xt,nttl,nnrm,nshr,
     3                                              mxpbs,npbs,prkin,
     4                                              ndkin,pryld,ndyld )
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
      end subroutine ummdp_kinematic_chaboche_ziegler
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     YOSHIDA-UEMORI KINEMATIC HARDENING LAW
c
      subroutine ummdp_kinematic_yoshida_uemori ( vk,dvkdp,dvkds,dvkdx,
     1                                            dvkdxt,p,s,x,xt,nttl,
     2                                            nnrm,nshr,mxpbs,npbs,
     3                                            prkin,ndkin,pryld,
     4                                            ndyld )
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
c************************************************************************
*
*     PRINT SUBROUTINES
*
************************************************************************
c
c     ummdp_print_elastic ( prela,ndela )
c       print elasticity parameters
c
c     ummdp_print_yield ( pryld,ndyld )
c       print yield criteria parameters
c
c     ummdp_print_isotropci ( prihd,ndihd )
c       print isotropic hardening law parameters
c
c     ummdp_print_kinematic ( prkin,ndkin,npbs )
c       print kinematic hardening law parameters
c
c     ummdp_print_rupture ( prrup,ndrup )
c       print uncoupled rupture criterion parameters
c     
c     ummdp_print_info
c       print informations for debug (info)
c
c     ummdp_print_inout
c       print informations for debug (input/output)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     PRINT ELASTICITY PARAMETERS
c
      subroutine ummdp_print_elastic ( prela,ndela )
c-----------------------------------------------------------------------
      implicit none
c
      integer ndela
      real*8 prela(ndela)
c
      integer ntela
c-----------------------------------------------------------------------
c
      ntela = nint(prela(1))
      write(6,*)
      write (6,'(4XA)') '>>> Elasticity'
      select case ( ntela )
      case ( 0 )
        write (6,'(8xA,I2)') '>> Young Modulus and Poisson Ratio'
        ! write (6,'(12xA,I2)') ' > ID:',ntela
        write (6,'(12xA,I2)') '. prela(1) =',prela(1)
        write (6,'(12xA,F10.1)') '. prela(2) =',prela(1+1)
        write (6,'(12xA,F5.2)') '. prela(3) =',prela(1+2)
      case ( 1 )
        write (6,*) 'Bulk Modulus and Modulus of Rigidity'
        write (6,*) 'Bulk modulus  =',prela(1+1)
        write (6,*) 'Shear modulus =',prela(1+2)
      end select
c
      return
      end subroutine ummdp_print_elastic
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     PRINT YIELD FUNCTION PARAMETERS
c
      subroutine ummdp_print_yield ( pryld,ndyld )
c-----------------------------------------------------------------------
      implicit none
c
      integer ndyld
      real*8 pryld(ndyld)
c
      integer i,j
      integer ntyld,n0,n
c-----------------------------------------------------------------------
c
      ntyld = pryld(1)
      write (6,*)
      write (6,*) '>> Yield Function',ntyld
      select case ( ntyld )
c
      case ( 0 )                                             ! von Mises
        write (6,*) 'von Mises'
c
      case ( 1 )                                             ! Hill 1948
        write (6,*) 'Hill 1948'
        write (6,*) 'F=',pryld(2)
        write (6,*) 'G=',pryld(3)
        write (6,*) 'H=',pryld(4)
        write (6,*) 'L=',pryld(5)
        write (6,*) 'M=',pryld(6)
        write (6,*) 'N=',pryld(7)
c
      case ( 2 )                                           ! Yld2004-18p
        write (6,*) 'Yld2004-18p'
        n0 = 1
        do i = 1,18
           n0 = n0 + 1
           write (6,*) 'a(',i,')=',pryld(n0)
        end do
        write (6,*) 'M=',pryld(1+18+1)
c
      case ( 3 )                                              ! CPB 2006
        write (6,*) 'CPB 2006'
        n0 = 1
        do i = 1,3
          do j = 1,3
            n0 = n0 + 1
            write (6,*) 'c(',i,',',j,')=',pryld(n0)
          end do
        end do
        do i = 4,6
          n0 = n0 + 1
          write (6,*) 'c(',i,',',i,')=',pryld(n0)
        end do
        n0 = n0 + 1
        write (6,*) 'a =',pryld(n0)
        n0 = n0 + 1
        write (6,*) 'ck=',pryld(n0)
c
      case ( 4 )                                 ! Karafillis-Boyce 1993
        write (6,*) 'Karafillis-Boyce 1993'
        n0 = 1
        do i = 1,6
          do j = i,6
            n0 = n0 + 1
            write (6,*) 'L(',i,',',j,') =',pryld(n0)
          end do
        end do
        n0 = n0 + 1
        write (6,*) 'k =',pryld(n0)
        n0 = n0 + 1
        write (6,*) 'c =',pryld(n0)
c
      case ( 5 )                                               ! Hu 2005
        write (6,*) 'Hu 2005'
        n0 = 1
        do i = 1,5
          n0 = n0 + 1
          write (6,*) 'X(',i,')=',pryld(n0)
        end do
        n0 = n0 + 1
        write (6,*) 'X(',7,')=',pryld(n0)
        do i = 1,3
          n0 = n0 + 1
          write (6,*) 'C(',i,')=',pryld(n0)
        end do
c
      case ( 6 )                                          ! Yoshida 2011
        write (6,*) 'Yoshida 2011'
        n0 = 1
        do i = 1,16
          n0 = n0+1
          write (6,*) 'c(',i,')=',pryld(n0)
        end do
c
      case ( -1 )                                    ! Gotoh Biquadratic
        write (6,*) 'Gotoh Biquadratic'
        do i = 1,9
          write (6,*) 'A(',i,')=',pryld(i+1)
        end do
c
      case ( -2 )                                           ! Yld2000-2d
        write (6,*) 'Yld2000-2d'
        do i = 1,8
          write (6,*) 'a(',i,')=',pryld(i+1)
        end do
        write (6,*) 'M=',pryld(9+1)
c
      case ( -3 )                                               ! Vegter
        write (6,*) 'Vegter '
        write (6,*) 'nf=',nint(pryld(2))
        write (6,*) 'f_bi0=',pryld(3)
        write (6,*) 'r_bi0=',pryld(4)
        do i = 0,nint(pryld(2))
          write (6,*) 'test angle=',90.0d0*float(i)/pryld(2)
          write (6,*) 'phi_un(',i,')=',pryld(4+i*4+1)
          write (6,*) 'phi_sh(',i,')=',pryld(4+i*4+2)
          write (6,*) 'phi_ps(',i,')=',pryld(4+i*4+3)
          write (6,*) 'omg(   ',i,')=',pryld(4+i*4+4)
        end do
c       do i = 1,7
c         write (6,*) 'phi_un(',i-1,')=',pryld(1+i   )
c         write (6,*) 'phi_sh(',i-1,')=',pryld(1+i+ 7)
c         write (6,*) 'phi_ps(',i-1,')=',pryld(1+i+14)
c         write (6,*) 'omg   (',i-1,')=',pryld(1+i+23)
c       end do
c       write (6,*)   'f_bi0=',pryld(1+22)
c       write (6,*)   'r_bi0=',pryld(1+23)
c       write (6,*)   'nf   =',nint(pryld(1+31))
c
      case ( -4 )                                             ! BBC 2005
        write (6,*) 'BBC 2005'
        write (6,*) 'k of order 2k',pryld(1+1)
        write (6,*) 'a=',pryld(1+2)
        write (6,*) 'b=',pryld(1+3)
        write (6,*) 'L=',pryld(1+4)
        write (6,*) 'M=',pryld(1+5)
        write (6,*) 'N=',pryld(1+6)
        write (6,*) 'P=',pryld(1+7)
        write (6,*) 'Q=',pryld(1+8)
        write (6,*) 'R=',pryld(1+9)
c
      case ( -5 )                                                ! Yld89
        write (6,*) 'Yld89'
        write (6,*) 'order M=',pryld(1+1)
        write (6,*) 'a      =',pryld(1+2)
        write (6,*) 'h      =',pryld(1+3)
        write (6,*) 'p      =',pryld(1+4)
c
      case ( -6 )                                             ! BBC 2008
        write (6,*) 'BBC 2008'
        write (6,*) 's      =',nint(pryld(1+1))
        write (6,*) 'k      =',nint(pryld(1+2))
        do i = 1,nint(pryld(1+1))
          write (6,*) 'i=',i
          n = 2 + (i-1)*8
          write (6,*) 'l_1=',pryld(n+1)
          write (6,*) 'l_2=',pryld(n+2)
          write (6,*) 'm_1=',pryld(n+3)
          write (6,*) 'm_2=',pryld(n+4)
          write (6,*) 'm_3=',pryld(n+5)
          write (6,*) 'n_1=',pryld(n+6)
          write (6,*) 'n_2=',pryld(n+7)
          write (6,*) 'n_3=',pryld(n+8)
        end do
c
      case ( -7 )                                            ! Hill 1990
        write (6,*) 'Hill 1990'
        write (6,*) 'a   =',pryld(1+1)
        write (6,*) 'b   =',pryld(1+2)
        write (6,*) 'tau =',pryld(1+3)
        write (6,*) 'sigb=',pryld(1+4)
        write (6,*) 'M   =',pryld(1+5)
      end select
c
      return
      end subroutine ummdp_print_yield
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     PRINT ISOTROPIC HARDENING LAW PARAMETERS
c
      subroutine ummdp_print_isotropic ( prihd,ndihd )
c-----------------------------------------------------------------------
      implicit none
c
      integer ndihd
      real*8 prihd(ndihd)
c
      integer ntihd
c-----------------------------------------------------------------------
c
      ntihd = nint(prihd(1))
      write (6,*)
      write (6,*) '>> Isotropic Hardening Law',ntihd
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
     &                    c*(e0+p)^en)'
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
      end subroutine ummdp_print_isotropic
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     PRINT KINEMATIC HARDENING LAW PARAMETERS
c
      subroutine ummdp_print_kinematic ( prkin,ndkin,npbs )
c-----------------------------------------------------------------------
      implicit none
c
      integer ndkin,npbs
      real*8 prkin(ndkin)
c
      integer i
      integer ntkin,n0
c-----------------------------------------------------------------------
c
      ntkin = nint(prkin(1))
      write (6,*)
      write (6,*) '>> Kinematic Hardening Law',ntkin
      select case ( ntkin )
      case ( 0 )                                ! No Kinematic Hardening
        write (6,*) 'No Kinematic Hardening'
c
      case ( 1 )                                                ! Prager
        write (6,*) 'Prager dX=(2/3)*c*{dpe}'
        write (6,*) 'c =',prkin(1+1)
c
      case ( 2 )                                               ! Ziegler
        write (6,*) 'Ziegler dX=dp*c*{{s}-{X}}'
        write (6,*) 'c =',prkin(1+1)
c
      case ( 3 )                          ! Armstrong & Frederick (1966)
        write (6,*) 'Armstrong-Frederick (1966)'
        write (6,*) 'dX=(2/3)*c*{dpe}-dp*g*{X}'
        write (6,*) 'c =',prkin(1+1)
        write (6,*) 'g =',prkin(1+2)
c
      case ( 4 )                                       ! Chaboche (1979)
        write (6,*) 'Chaboche (1979)'
        write (6,*) 'dx(j)=c(j)*(2/3)*{dpe}-dp*g(j)*{x(j)}'
        write (6,*) 'no. of x(j) =',npbs
        do i = 1,npbs
          n0 = 1+(i-1)*2
          write (6,*) 'c(',i,')=',prkin(1+n0+1)
          write (6,*) 'g(',i,')=',prkin(1+n0+2)
        end do
c
      case ( 5 )                       ! Chaboche (1979) - Ziegler Model
        write (6,*) 'Chaboche (1979) - Ziegler Model'
        write (6,*) 'dx(j)=((c(j)/se)*{{s}-{X}}-g(j)*{x(j)})*dp'
        write (6,*) 'no. of x(j) =',npbs
        do i = 1,npbs
          n0 = 1+(i-1)*2
          write (6,*) 'c(',i,')=',prkin(1+n0+1)
          write (6,*) 'g(',i,')=',prkin(1+n0+2)
        end do
c
      case ( 6 )                                        ! Yoshida-Uemori
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
      end subroutine ummdp_print_kinematic
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     PRINT UNCOUPLED RUPTURE CRITERION PARAMETERS
c
      subroutine ummdp_print_rupture ( prrup,ndrup )
c-----------------------------------------------------------------------
      implicit none
c
      integer ndrup
      real*8 prrup(ndrup)
c
      integer ntrup
c-----------------------------------------------------------------------
c
      ntrup = nint(prrup(1))
      write (6,*)
      write (6,*) '>> Uncoupled Rupture Criterion',ntrup
c
      select case ( ntrup )
c
      case ( 0 ) 										   	! No Uncoupled Rupture Criterion
        write (6,*) 'No Uncoupled Rupture Criterion'
c
      case ( 1 ) 														 ! Equivalent Plastic Strain
        write (6,*) 'Equivalent Plastic Strain'
        write (6,*) 'W=int[dp]'
        write (6,*) 'Wl=',prrup(3)
c
      case ( 2 )  																 ! Cockroft and Latham
        write (6,*) 'Cockroft and Latham'
        write (6,*) 'W=int[(sp1/se)*dp]'
        write (6,*) 'Wl=',prrup(3)
c
      case ( 3 ) 																	     ! Rice and Tracey
        write (6,*) 'Rice and Tracey'
        write (6,*) 'W=int[exp(1.5*sh/se)*dp]'
        write (6,*) 'Wl=',prrup(3)
c
      case ( 4 ) 																						     ! Ayada
        write (6,*) 'Ayada'
        write (6,*) 'W=int[(sh/se)*dp]'
        write (6,*) 'Wl=',prrup(3)
c
      case ( 5 ) 																							  ! Brozzo
        write (6,*) 'Brozzo'
        write (6,*) 'W=int[(2/3)*(sp1/(sp1-se))*dp]'
        write (6,*) 'Wl=',prrup(3)
c
      case ( 6 ) 	   														 ! Forming Limit Diagram
        write (6,*) 'Forming Limit Diagram'
        write (6,*) 'W=e1/e1(fld)'
        write (6,*) 'Wl=',prrup(3)
c
      end select
c
      return
      end subroutine ummdp_print_rupture
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     PRINT INFORMATIONS FOR DEBUG (INFO)
c
      subroutine ummdp_print_info ( inc,nnrm,nshr )
c-----------------------------------------------------------------------
      implicit none
c
			integer inc,nnrm,nshr
c
      integer nttl,nerr
c
			integer ne,ip,lay
      common /jancae1/ne,ip,lay
c-----------------------------------------------------------------------
c
      nttl = nnrm + nshr
c
      write (6,*) '----- JANCAE.UMMDp Debug Info -----'
      write (6,*) 'increment=',inc
      write (6,*) 'elem,ip,lay=',ne,ip,lay
      write (6,*) 'nttl,nnrm,nshr=',nttl,nnrm,nshr
      nerr = 0
      if ( nnrm == 3 ) then
        if ( nshr == 3 ) then
          write (6,*) '3d solid element'
        else if ( nshr == 1 ) then
          write (6,*) 'plane strain or axi-sym solid element'
        else
          nerr = nerr + 1
        end if
      else if ( nnrm == 2 ) then
        if ( nshr == 1 ) then
          write (6,*) 'plane stress or thin shell element'
        else if ( nshr == 3 ) then
          write (6,*) 'thick shell element'
        else
          nerr = nerr + 1
        end if
      else
        nerr = nerr + 1
      end if
      if ( nerr /= 0 ) then
        write (6,*) 'no supported element type',nnrm,nshr
        call ummdp_exit ( 9000 )
      end if
c
      return
      end subroutine ummdp_print_info
c
c
c
************************************************************************
c     PRINT INFORMATIONS FOR DEBUG (INPUT/OUTPUT)
c
      subroutine ummdp_print_inout ( io,s,de,d,nttl,stv,nstv )
c-----------------------------------------------------------------------
      implicit none
c
			integer io,nttl,nstv
			real*8 s(nttl),stv(nstv),de(nttl),d(nttl,nttl)
c
      character*32 text
c-----------------------------------------------------------------------
c
      if ( io == 0 ) then
        text = 'initial stresses'
      else
        text = 'updated stresses'
      end if
      call ummdp_utility_print1 ( text,s,nttl )
c
      if ( io == 0 ) then
        text = 'initial internal state var.'
      else
        text = 'updated internal state var.'
      end if
      call ummdp_utility_print1 ( text,stv,nstv )
c
      if ( io == 0 ) then
        text = 'driving strain increment'
        call ummdp_utility_print1 ( text,de,nttl )
      else
        text = 'tangent modulus matrix'
        call ummdp_utility_print2 ( text,d,nttl,nttl )
      end if
c
      return
      end subroutine ummdp_print_inout
c
c
c************************************************************************
*
*     UNCOUPLED RUPTURE CRITERIA
*
c***********************************************************************
c
c      0 : No Rupture Criterion
c
c      1 : Equivalent Plastic Strain
c      2 : Cockroft and Latham
c      3 : Rice and Tracey
c      4 : Ayada
c      5 : Brozzo
c      6 : Forming Limit Diagram (only plane-stress)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     CALCULATE UNCOUPLED RUPTURE CRITERIA
c
      subroutine ummdp_rupture ( ntens,sdv,nsdv,uvar2,uvar1,nuvarm,
     1                           jrcd,jmac,jmatyp,matlayo,laccfla,
     2                           nt,ndrup,prrup )
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension UVAR1(NUVARM),jmac(*),jmatyp(*)
      dimension sdv(nsdv),uvar2(nuvarm),prrup(ndrup)
      real*8 lim,wlimnorm
c-----------------------------------------------------------------------
c
c      prrup(1) : criteria id
c      prrup(2) : flag to terminate analysis if limit is reached
c      prrup(3) : rupture limit
c
c 																		       ---- rupture criteria limit
      lim = prrup(3)
c                                           ---- select rupture criteria
      ntrup = nint(prrup(1))
      select case ( ntrup )
c
      case ( 0 )                                  ! No Rupture Criterion
        return
c
      case ( 1 )                             ! Equivalent Plastic Strain
        call ummdp_rupture_eqvstrain ( sdv,nsdv,uvar2,uvar1,nuvarm,
     1                                 nt,lim,wlimnorm )
c
      case ( 2 )                                   ! Cockroft and Latham
        call ummdp_rupture_cockroft ( sdv,nsdv,uvar2,uvar1,nuvarm,
     1                                jrcd,jmac,jmatyp,matlayo,laccfla,
     2                                nt,lim,wlimnorm )
c
      case ( 3 )                                       ! Rice and Tracey
        call ummdp_rupture_rice ( sdv,nsdv,uvar2,uvar1,nuvarm,
     1                            jrcd,jmac,jmatyp,matlayo,laccfla,
     2                            nt,lim,wlimnorm )
c
      case ( 4 )                                                 ! Ayada
        call ummdp_rupture_ayada ( sdv,nsdv,uvar2,uvar1,nuvarm,
     1                             jrcd,jmac,jmatyp,matlayo,laccfla,
     2                             nt,lim,wlimnorm )
c
      case ( 5 )                                                ! Brozzo
        call ummdp_rupture_brozzo ( sdv,nsdv,uvar2,uvar1,nuvarm,
     1                              jrcd,jmac,jmatyp,matlayo,laccfla,
     2                              nt,lim,wlimnorm )
c
      case ( 6 )                                 ! Forming Limit Diagram
        call ummdp_rupture_fld ( ntens,uvar2,uvar1,nuvarm,
     1                           jrcd,jmac,jmatyp,matlayo,laccfla,
     2                           nt,lim,wlimnorm )
c
      case default
        write (6,*) 'error in ummdp_rupture'
        write (6,*) 'ntrup error :',ntrup
        call ummdp_exit ( 9000 )
      end select
c
c                    ---- terminate analysis if rupture limit is reached
      end = nint(prrup(2))
      if ( end == 1 ) then
        if ( wlimnorm >= 1.0d0 ) then 
          write (6,*) 'analysis terminated by rupture criterion'
          write (6,*) 'stop in uvrm.'
          call ummdp_exit( 10000 )
        end if
      end if
c
      return
      end
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     Equivalent Plastic Strain
c
      subroutine ummdp_rupture_eqvstrain ( sdv,nsdv,uvar2,uvar1,nuvarm,
     1                                     nt,lim,wlimnorm )
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension sdv(nsdv),uvar2(nuvarm),uvar1(nuvarm)
      real*8 lim,peeq
c-----------------------------------------------------------------------
c
c     nuvarm : 2
c
c     uvar(2+nt+1) : equivalent plastic strain
c     uvar(2+nt+2) : rupture parameter normalised by critical value
c
c                                            ---- get uvar before update
      wlimnorm  = uvar1(2+nt+2)
c
c                                              ---- get sdv after update
      peeq = sdv(1)
c
c                                                       ---- update uvar
      uvar2(2+nt+1) = peeq
      uvar2(2+nt+2) = peeq / lim
c
      return
      end subroutine ummdp_rupture_eqvstrain
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     Cockroft and Latham
c
      subroutine ummdp_rupture_cockroft ( sdv,nsdv,uvar2,uvar1,nuvarm,
     1                                    jrcd,jmac,jmatyp,matlayo,
     2                                    laccfla,nt,lim,wlimnorm )
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension UVAR1(NUVARM),jmac(*),jmatyp(*)
      dimension array(15),jarray(15)
      dimension sdv(nsdv),uvar2(nuvarm)
      character*3 flgray(15)
      real*8 lim,wlimnorm
      real*8 se1,peeq1,maxsp1,maxsp1se1,wlim1
      real*8 se2,peeq2,maxsp2,maxsp2se2,wlim2
c-----------------------------------------------------------------------
c
c     nuvarm : 4
c
c     uvar(1)      : equivalent stress
c     uvar(2+nt+1) : equivalent plastic strain
c     uvar(2+nt+2) : maximum principal stress
c     uvar(2+nt+3) : rupture parameter
c     uvar(2+nt+4) : rupture parameter normalised by critical value
c
c                                            ---- get uvar before update
      se1    = uvar1(1)
      peeq1  = uvar1(2+nt+1)
      maxsp1 = uvar1(2+nt+2)
      wlim1  = uvar1(2+nt+3)
c
      wlimnorm  = uvar1(2+nt+4)
c
c                                     ---- get sdv and uvar after update
      se2 = uvar2(1)
      peeq2 = sdv(1)
c
c                                 ---- get principal stress after update
      call getvrm ('SP',array,jarray,flgray,jrcd,jmac,jmatyp,
     1                  matlayo,laccfla )
      if ( jrcd /= 0 ) then
        write (6,*) 'request error in uvarm for sp'
        write (6,*) 'stop in uvrm.'
        call ummdp_exit ( 9000 )
      end if
      maxsp2 = array(3)
c
c                                                 ---- rupture criterion
      maxsp1se1 = 0.0d0
      maxsp2se2 = 0.0d0
      if ( se1 > 0.0d0 ) maxsp1se1 = maxsp1 / se1
      if ( se2 > 0.0d0 ) maxsp2se2 = maxsp2 / se2
c
      wlim2 = wlim1 + (maxsp2se2+maxsp1se1)*(peeq2-peeq1)/2.0d0
c
c                                                       ---- update uvar
      uvar2(2+nt+1) = peeq2
      uvar2(2+nt+2) = maxsp2
      uvar2(2+nt+3) = wlim2
      uvar2(2+nt+4) = wlim2/lim
c
      return
      end subroutine ummdp_rupture_cockroft
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     Rice and Tracey
c
      subroutine ummdp_rupture_rice ( sdv,nsdv,uvar2,uvar1,nuvarm,
     1                                jrcd,jmac,jmatyp,matlayo,laccfla,
     2                                nt,lim,wlimnorm )
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension UVAR1(NUVARM),jmac(*),jmatyp(*)
      dimension array(15),jarray(15)
      dimension sdv(nsdv),uvar2(nuvarm)
      character*3 flgray(15)
      real*8 lim,wlimnorm
      real*8 se1,peeq1,shyd1,shyd1se1,wlim1
      real*8 se2,peeq2,shyd2,shyd2se2,wlim2
c-----------------------------------------------------------------------
c
c     nuvarm : 4
c
c     uvar(1)      : equivalent stress
c     uvar(2+nt+1) : equivalent plastic strain
c     uvar(2+nt+2) : hydrostatic stress
c     uvar(2+nt+3) : rupture parameter
c     uvar(2+nt+4) : rupture parameter normalised by critical value
c
c                                            ---- get uvar before update
      se1   = uvar1(1)
      peeq1 = uvar1(2+nt+1)
      shyd1 = uvar1(2+nt+2)
      wlim1 = uvar1(2+nt+3)
c
      wlimnorm  = uvar1(2+nt+4)
c
c                                     ---- get sdv and uvar after update
      se2 = uvar2(1)
      peeq2 = sdv(1)
c
c                               ---- get hydrostatic stress after update
      call getvrm ( 'SINV',array,jarray,flgray,jrcd,jmac,jmatyp,
     1                     matlayo,laccfla )
      if ( jrcd /= 0 ) then
        write (6,*) 'request error in uvarm for sinv'
        write (6,*) 'stop in uvrm.'
        call ummdp_exit ( 9000 )
      end if
      shyd2 = -array(3)
c
c                                                 ---- rupture criterion
      shyd1se1 = 0.0d0
      shyd2se2 = 0.0d0
      if ( se1 > 0.0d0 ) shyd1se1 = exp(1.5d0*shyd1/se1)
      if ( se2 > 0.0d0 ) shyd2se2 = exp(1.5d0*shyd2/se2)
c
      wlim2 = wlim1 + (shyd1se1+shyd2se2)*(peeq2-peeq1)/2.0d0
c
c                                                       ---- update uvar
      uvar2(2+nt+1) = peeq2
      uvar2(2+nt+2) = shyd2
      uvar2(2+nt+3) = wlim2
      uvar2(2+nt+4) = wlim2/lim
c
      return
      end subroutine ummdp_rupture_rice
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     Ayada
c
      subroutine ummdp_rupture_ayada ( sdv,nsdv,uvar2,uvar1,nuvarm,
     1                                 jrcd,jmac,jmatyp,matlayo,laccfla,
     2                                 nt,lim,wlimnorm )
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension UVAR1(NUVARM),jmac(*),jmatyp(*)
      dimension array(15),jarray(15)
      dimension sdv(nsdv),uvar2(nuvarm)
      character*3 flgray(15)
      real*8 lim,wlimnorm
      real*8 se1,peeq1,shyd1,shyd1se1,wlim1
      real*8 se2,peeq2,shyd2,shyd2se2,wlim2
c-----------------------------------------------------------------------
c
c     nuvarm : 4
c
c     uvar(1)      : equivalent stress
c     uvar(2+nt+1) : equivalent plastic strain
c     uvar(2+nt+2) : hydrostatic stress
c     uvar(2+nt+3) : rupture parameter
c     uvar(2+nt+4) : rupture parameter normalised by critical value
c
c                                            ---- get uvar before update
      se1   = uvar1(1)
      peeq1 = uvar1(2+nt+1)
      shyd1 = uvar1(2+nt+2)
      wlim1 = uvar1(2+nt+3)
c
      wlimnorm  = uvar1(2+nt+4)
c
c                                     ---- get sdv and uvar after update
      se2 = uvar2(1)
      peeq2 = sdv(1)
c
c                               ---- get hydrostatic stress after update
      call getvrm ( 'SINV',array,jarray,flgray,jrcd,jmac,jmatyp,
     1                     matlayo,laccfla )
      if ( jrcd /= 0 ) then
        write (6,*) 'request error in uvarm for sinv'
        write (6,*) 'stop in uvrm.'
        call ummdp_exit ( 9000 )
      end if
      shyd2 = -array(3)
c
c                                                 ---- rupture criterion
      shyd1se1 = 0.0d0
      shyd2se2 = 0.0d0
      if ( se1 > 0.0d0 ) shyd1se1 = shyd1 / se1
      if ( se2 > 0.0d0 ) shyd2se2 = shyd2 / se2
c
      wlim2 = wlim1 + (shyd1se1+shyd2se2)*(peeq2-peeq1)/2.0d0
c
c                                                       ---- update uvar
      uvar2(2+nt+1) = peeq2
      uvar2(2+nt+2) = shyd2
      uvar2(2+nt+3) = wlim2
      uvar2(2+nt+4) = wlim2/lim
c
      return
      end subroutine ummdp_rupture_ayada
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     Brozzo
c
      subroutine ummdp_rupture_brozzo ( sdv,nsdv,uvar2,uvar1,nuvarm,
     1                                  jrcd,jmac,jmatyp,matlayo,
     2                                  laccfla,nt,lim,wlimnorm )
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension UVAR1(NUVARM),jmac(*),jmatyp(*)
      dimension array(15),jarray(15)
      dimension sdv(nsdv),uvar2(nuvarm)
      character*3 flgray(15)
      real*8 lim,wlimnorm
      real*8 se1,peeq1,shyd1,maxsp1,maxsp1shyd1,wlim1
      real*8 se2,peeq2,shyd2,maxsp2,maxsp2shyd2,wlim2
c-----------------------------------------------------------------------
c
c     nuvarm : 5
c
c     uvar(2+nt+1) : equivalent plastic strain
c     uvar(2+nt+2) : maximum principal stress
c     uvar(2+nt+3) : hydrostatic stress
c     uvar(2+nt+4) : rupture parameter
c     uvar(2+nt+5) : rupture parameter normalised by critical value
c
c                                            ---- get uvar before update
      peeq1  = uvar1(2+nt+1)
      maxsp1 = uvar1(2+nt+2)
      shyd1  = uvar1(2+nt+3)
      wlim1  = uvar1(2+nt+4)
c
      wlimnorm  = uvar1(2+nt+5)
c
c                                              ---- get sdv after update
      peeq2 = sdv(1)
c
c                                 ---- get principal stress after update
      call getvrm ('SP',array,jarray,flgray,jrcd,jmac,jmatyp,
     1                  matlayo,laccfla )
      if ( jrcd /= 0 ) then
        write (6,*) 'request error in uvarm for sp'
        write (6,*) 'stop in uvrm.'
        call ummdp_exit ( 9000 )
      end if
      maxsp2 = array(3)
c
c                               ---- get hydrostatic stress after update
      call getvrm ('SINV',array,jarray,flgray,jrcd,jmac,jmatyp,
     1                    matlayo,laccfla )
      if ( jrcd /= 0 ) then
        write (6,*) 'request error in uvarm for sinv'
        write (6,*) 'stop in uvrm.'
        call ummdp_exit ( 9000 )
      end if
      shyd2 = -array(3)
c
c                                                 ---- rupture criterion
      maxsp1shyd1 = 0.0d0
      maxsp2shyd2 = 0.0d0
      if ( shyd1 > 0.0d0 ) then
        maxsp1shyd1 = (2.0d0/3.0d0) * maxsp1 / (maxsp1-shyd1)
      end if
      if ( shyd2 > 0.0d0 ) then 
        maxsp2shyd2 = (2.0d0/3.0d0) * maxsp2 / (maxsp2-shyd2)
      end if
c
      wlim2 = wlim1 + (maxsp2shyd2+maxsp1shyd1)*(peeq2-peeq1)/2.0d0
c
c                                                       ---- update uvar
      uvar2(2+nt+1) = peeq2
      uvar2(2+nt+2) = maxsp2
      uvar2(2+nt+3) = shyd2
      uvar2(2+nt+4) = wlim2
      uvar2(2+nt+5) = wlim2/lim
c
      return
      end subroutine ummdp_rupture_brozzo
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     Forming Limit Diagram (FLD)
c       . only plane-stress formulation
c
      subroutine ummdp_rupture_fld ( ntens,uvar2,uvar1,nuvarm,
     1                               jrcd,jmac,jmatyp,matlayo,laccfla,
     2                               nt,lim,wlimnorm )
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      parameter (mxflc=50)
c
      dimension jmac(*),jmatyp(*)
      dimension array(15),jarray(15)
      dimension uvar2(nuvarm),uvar1(nuvarm)
      character*3 flgray(15)
      real*8 lim,wlimnorm
      real*8 e1,e2,e1fld,wlim
      real*8 le(3,3),es(3),ev(3,3)
      dimension dum1(0),dum2(mxflc,0)
      dimension fld1(mxflc),fld2(mxflc)
c-----------------------------------------------------------------------
c
c     nuvarm : 5
c
c     uvar(2+nt+1) : maximum principal strain
c     uvar(2+nt+2) : minimum principal strain
c     uvar(2+nt+3) : projection of major principal strain on FLD
c     uvar(2+nt+4) : rupture parameter
c     uvar(2+nt+5) : rupture parameter normalised by critical value
c
c                                            ---- get uvar before update
      wlimnorm  = uvar1(2+nt+6)
c      
c                                 ---- get principal strain after update
      if ( ntens == 3 ) then
        call getvrm ( 'LEP',array,jarray,flgray,jrcd,jmac,jmatyp,
     1                      matlayo,laccfla )
        if ( jrcd /= 0 ) then
          write (6,*) 'request error in uvarm for lep'
          write (6,*) 'stop in uvrm.'
          call ummdp_exit ( 9000 )
        end if
        e2 = array(1)
        e1 = array(2)
c                               ---- get logarithmic strain after update
      else
        write (6,*) 'request error in uvarm for fld, only plane-stress'
        write (6,*) 'stop in uvrm.'
        call ummdp_exit ( 9000 )
!         call getvrm ( 'LE',array,jarray,flgray,jrcd,jmac,jmatyp,
!      &                     matlayo,laccfla )
!         if ( jrcd /= 0 ) then
!           write (6,*) 'request error in uvarm for le'
!           write (6,*) 'stop in uvrm.'
!           call ummdp_exit ( 9000 )
!         end if
! c                                            ---- assemble strain tensor
!         call ummdp_clear2 ( le,3,3 )
!         le(1,1) = array(1)
!         le(2,2) = array(2)
!         le(3,3) = array(3)
!         le(1,2) = array(4)/2
!         le(1,3) = array(5)/2
!         le(2,3) = array(6)/2
!         le(2,1) = le(1,2)
!         le(3,1) = le(1,3)
!         le(3,2) = le(2,3)
! c                           ---- strain tensor eigen- values and vectors
!         call ummdp_clear1 ( es,3 )
!         call ummdp_clear2 ( ev,3,3 )
!         call ummdp_eigen_sym3 ( es,ev,le )
!         e2 = es(2)
!         e1 = es(1)
      end if
c
c                                     ---- activate fld table collection
      call settablecollection ( 'FLD',jerror )

      if ( jerror /= 0 ) then
        write (6,*) 'request error in uvarm for table collection fld'
        write (6,*) 'stop in uvrm.'
        call ummdp_exit ( 9000 )
      end if
c	                     		                        ---- get fld E1 values
      call getpropertytable ( 'FLD1',dum1,dum1,dum1,nfld,fld1,dum2,0,
     1                               jerror )
      if ( jerror /= 0 ) then
        write (6,*) 'request error in uvarm for property table fld1'
        write (6,*) 'stop in uvrm.'
        call ummdp_exit ( 9000 )
      end if
c                                                 ---- get fld E2 values
      call getpropertytable ( 'FLD2',dum1,dum1,dum1,nfld,fld2,dum2,0,
     1                               jerror )
      if ( jerror /= 0 ) then
        write (6,*) 'request error in uvarm for property table fld2'
        write (6,*) 'stop in uvrm.'
        call ummdp_exit ( 9000 )
      end if
c
c                         ---- linear extra/inter -polation of E1 on FLD
      n = nfld
c                              --- linear extrapolation on the left side
      if ( e2 < fld2(1) ) then
        e1fld = fld1(2) + ( (e2-fld2(2)) / (fld2(1)-fld2(2)) )
     1          * ( fld1(1) - fld1(2) )
c
c                             --- linear extrapolation on the right side
      else if ( e2 > fld2(n) ) then
        e1fld = fld1(n-1) + ( (e2-fld2(n-1)) / (fld2(n)-fld2(n-1)) ) 
     1          * ( fld1(n) - fld1(n-1) )
c
c                                  --- linear interpolation inside range
      else
        k = 0
        do i = 1,n-1
          if ( ( e2 >= fld2(i) ) .and. ( e2 <= fld2(i+1) ) ) then
            k = i
          end if
        end do
        e1fld = fld1(k) + ( fld1(k+1) - fld1(k) )
     1          * ( (e2-fld2(k)) / (fld2(k+1)-fld2(k)) )
      end if
c
c                                                 ---- rupture criterion
      wlim = e1/e1fld
c
c                                                       ---- update uvar
      uvar2(2+nt+1) = e1
      uvar2(2+nt+2) = e2
      uvar2(2+nt+3) = e1fld
      uvar2(2+nt+4) = wlim
      uvar2(2+nt+5) = wlim/lim
c
      return
      end subroutine ummdp_rupture_fld
c
c
c************************************************************************
*
*     UTILITY SUBROUTINES
*
************************************************************************
c
c     ummdp_utility_clear1 ( a,n )
c       clear 1st order vector
c
c     ummdp_utility_clear2 ( a,n,m )
c       clear 2nd order matrix
c
c     ummdp_utility_clear3 ( a,n,m,l )
c       clear 3rd order tensor
c
c     ummdp_utility_setunitm ( a,n )
c       set unit 2nd order matrix
c
c     ummdp_utility_print1 ( text,a,n )
c       print vector with text
c
c     ummdp_utility_print2 ( text,a,n,m )
c       print matrix with text
c
c     ummdp_utility_mv (v,a,u,nv,nu)
c       mutiply matrix and vector
c
c     ummdp_utility_mm (a,b,c,na1,na2,nbc)
c       mutiply matrix and matrix
c
c     ummdp_utility_vvs ( s,u,v,n )
c       calculate scalar product of vectors 
c
c     ummdp_utility_minv ( b,a,n,d )
c       calculate inverse matrix using lu decomposition
c
c         ummdp_utility_ludcmp( a,n,indx,d )
c           lu decomposition
c         ummdp_utility_lubksb(a,n,indx,b)
c           lu backward substitution
c         ummdp_utility_minv2 ( b,a,deta )
c           calculate inverse matrix 2x2
c         ummdp_utility_minv3 ( b,a,deta )
c           calculate inverse matrix 3x3
c
c     ummdp_utility_eigen_sym3 ( es,ev,a )
c       calculate eigenvalues and eigenvectors by jacobi method
c
c     ummdp_utility_file_exist ( flname )
c       checking existence of files
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     CLEAR 1st ORDER VECTOR A(N)
c
      subroutine ummdp_utility_clear1 ( a,n )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: n
c
      real*8,intent(inout) :: a(n)
c
      integer i
c-----------------------------------------------------------------------
c
      do i = 1,n
        a(i) = 0.0d0
      end do
c
      return
      end subroutine ummdp_utility_clear1
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     CLEAR 2ND ORDER MATRIX
c
      subroutine ummdp_utility_clear2 ( a,n,m )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: n,m
c
      real*8,intent(inout) :: a(n,m)
c
      integer i,j
c-----------------------------------------------------------------------
c
      do i = 1,n
        do j = 1,m
          a(i,j) = 0.0d0
        end do
      end do
c
      return
      end subroutine ummdp_utility_clear2
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     CLEAR 3RD ORDER MATRIX
c
      subroutine ummdp_utility_clear3 ( a,n,m,l )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: n,m,l
c
      real*8,intent(inout) ::  a(n,m,l)
c
      integer i,j,k
c-----------------------------------------------------------------------
c
      do i = 1,n
        do j = 1,m
          do k = 1,l
            a(i,j,k) = 0.0d0
          end do
        end do
      end do
c
      return
      end subroutine ummdp_utility_clear3
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     SET UNIT 2ND ORDER MATRIX
c
      subroutine ummdp_utility_setunitm ( a,n )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: n
c
      real*8,intent(out) :: a(n,n)
c
      integer i
c-----------------------------------------------------------------------
c
      call ummdp_utility_clear2 ( a,n,n )
      do i = 1,n
        a(i,i) = 1.0d0
      end do
c
      return
      end subroutine ummdp_utility_setunitm
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     PRINT VECTOR WITH TEXT
c
      subroutine ummdp_utility_print1 ( text,a,n )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer     ,intent(in) :: n
      real*8      ,intent(in) :: a(n)
      character*32,intent(in) :: text
c
      integer i
c-----------------------------------------------------------------------
c
      write (6,*) text
      write (6,9000) (a(i),i=1,n)
 9000 format (6e16.8)
c
      return
      end subroutine ummdp_utility_print1
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     PRINT MATRIX WITH TEXT
c
      subroutine ummdp_utility_print2 ( text,a,n,m )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer     ,intent(in) :: n,m
      real*8      ,intent(in) :: a(n,m)
      character*32,intent(in) :: text
c
      integer i,j
c-----------------------------------------------------------------------
      write (6,*) text
      do i = 1,n
        write (6,9000) (a(i,j),j=1,m)
      end do
 9000 format (6e16.8)
c
      return
      end subroutine ummdp_utility_print2
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     MULTIPLY MATRIX AND VECTOR
c
      subroutine ummdp_utility_mv ( v,a,u,nv,nu )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: nv,nu
      real*8 ,intent(in) :: u(nu)
      real*8 ,intent(in) :: a(nv,nu)
c
      real*8,intent(out) :: v(nv)
c
      integer i,j
c-----------------------------------------------------------------------
c
      call ummdp_utility_clear1 ( v,nv )
      do i = 1,nv
        do j = 1,nu
          v(i) = v(i) + a(i,j)*u(j)
        end do
      end do
c
      return
      end subroutine ummdp_utility_mv
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     MULTIPLY MATRIX AND MATRIX
c     
      subroutine ummdp_utility_mm ( a,b,c,na1,na2,nbc )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: na1,na2,nbc
      real*8 ,intent(in) :: b(na1,nbc),c(nbc,na2)
c
      real*8,intent(out) :: a(na1,na2)
c
      integer i,j,k
c-----------------------------------------------------------------------
c
      call ummdp_utility_clear2 ( a,na1,na2 )
      do i = 1,na1
        do j = 1,na2
          do k = 1,nbc
            a(i,j) = a(i,j) + b(i,k)*c(k,j)
          end do
        end do
      end do
c
      return
      end subroutine ummdp_utility_mm
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     CALCULATE SCALAR PRODUCT OF VECTORS
c
      subroutine ummdp_utility_vvs ( s,u,v,n )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: n
      real*8 ,intent(in) :: v(n),u(n)
c
      real*8,intent(out) :: s
c
      integer i
c-----------------------------------------------------------------------
c
      s = 0.0d0
      do i = 1,n
        s = s + u(i)*v(i)
      end do
c
      return
      end subroutine ummdp_utility_vvs
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     CALCULATE INVERSE MATRIX USING LU DECOMPOSITION
c
c       Ref.: http://astr-www.kj.yamagata-u.ac.jp/~shibata/kbg/
c
      subroutine ummdp_utility_minv ( b,a,n,d )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer n
      real*8 d
      real*8 a(n,n),b(n,n)
c
      integer i,j,k
      real*8 eps,anorm,ani
      real*8 indx(n),y(n),c(n,n),aorg(n,n)
      character*32 text
      logical check
c-----------------------------------------------------------------------
c
      check = .false.
      eps = 1.0d-36
c
      do i = 1,n
        do j = 1,n
          aorg(i,j) = a(i,j)
        end do
      end do
c
      anorm = 0.0
      do i = 1,n
        do j = 1,n
          if ( anorm < abs(a(i,j)) ) anorm = abs(a(i,j))
        end do
      end do
      do i = 1,n
        do j = 1,n
          a(i,j) = a(i,j) / anorm
        end do
      end do
c
      if ( n == 2 ) then
        call ummdp_utility_minv2 ( b,a,d,eps )
        goto 100
      else if ( n == 3 ) then
        call ummdp_utility_minv3 ( b,a,d,eps )
        goto 100
      end if
c
      call ummdp_utility_ludcmp ( a,n,indx,d,eps )
c                                                 ---- check determinant
      if ( abs(d) <= eps ) then
         write (6,*) 'determinant det[a] error',d
         write (6,*) 'stop in minv'
         call ummdp_exit ( 9000 ) 
      end if
c                                                            ---- B=A^-1
      do j = 1,n
        call ummdp_utility_clear1 ( y,n )
        y(j) = 1.0d0
        call ummdp_utility_lubksb ( a,n,indx,y,eps )
        do i = 1,n
          b(i,j) = y(i)
        end do
      end do
c
  100 continue
      ani = 1.0d0/anorm
      do i = 1,n
        do j = 1,n
          a(i,j) = aorg(i,j)
          b(i,j) = b(i,j) * ani
        end do
      end do
c                                                             ---- check
      if ( check ) then
        write (6,*) 'check inverse matrix',n
        text = 'original matrix [A]'
        call ummdp_utility_print2 ( text,a,n,n )
        text = 'inversed matrix [A]^-1'
        call ummdp_utility_print2 ( text,b,n,n )
        call ummdp_utility_clear2 ( c,n,n )
        do i = 1,n
          do j = 1,n
            do k = 1,n
              c(i,j) = c(i,j) + b(i,k)*a(k,j)
            end do
          end do
        end do
        text = '[A]^-1*[A]=[I] ?'
        call ummdp_utility_print2 ( text,c,n,n )
      end if
c
      return
      end subroutine ummdp_utility_minv
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     LU DECOMPOSITION
c
      subroutine ummdp_utility_ludcmp ( a,n,indx,d,eps )
c
c-----------------------------------------------------------------------
      implicit none
c
			integer n
			integer indx(n)
			real*8 d,eps
      real*8 a(n,n)
c
			integer i,j,k
			integer imax
			real*8 aamax,sum,dum,ajj
			real*8 vtemp(n)
      character*32 text
c-----------------------------------------------------------------------
c
      d = 1.0d0
      do i = 1,n
        aamax = 0.0d0
        do j = 1,n
          if ( abs(a(i,j)) > aamax ) aamax = abs(a(i,j))
        end do
        if ( aamax <= eps ) then
          write (6,*) 'singular matrix in ummdp_ludcmp'
          text = 'matrix detail'
          call ummdp_utility_print2 ( text,a,n,n )
          call ummdp_exit ( 9000 )
        end if
        vtemp(i) = 1.0d0 / aamax
      end do
c
      do j = 1,n
        do i = 1,j-1
          sum = a(i,j)
          do k = 1,i-1
            sum = sum - a(i,k)*a(k,j)
          end do
          a(i,j) = sum
        end do
        aamax = 0.0d0
        do i = j,n
          sum = a(i,j)
          do k = 1,j-1
            sum = sum - a(i,k)*a(k,j)
          end do
          a(i,j) = sum
          dum = vtemp(i)*abs(sum)
          if ( dum >= aamax ) then
            imax = i
            aamax = dum
          end if
        end do
        if ( j /= imax ) then
          do k = 1,n
            dum = a(imax,k)
            a(imax,k) = a(j,k)
            a(j,k) = dum
          end do
          d = -d
          vtemp(imax) = vtemp(j)
        end if
        indx(j) = imax
c       if ( abs(a(i,j)) <= eps ) a(i,j) = eps     !2010.07.02 c.out
        if ( j /= n ) then
          ajj = a(j,j)                               !2010.07.02 add
          if ( abs(ajj) <= eps ) ajj = eps         !2010.07.02 add
          dum = 1.0d0 / ajj                          !2010.07.02 mod
          do i = j+1,n
            a(i,j) = a(i,j) * dum
          end do
        end if
      end do
c                                                 ---- get the det. of A
      do j = 1,n
        d = d*a(j,j)
      end do
c
      return
      end subroutine ummdp_utility_ludcmp
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     LU BACKWARD SUBSTITUTION
c
      subroutine ummdp_utility_lubksb ( a,n,indx,b,eps )
c-----------------------------------------------------------------------
      implicit none
c
			integer n
			integer indx(n)
			real*8 eps
			real*8 a(n,n),b(n)
c
			integer i,j
			integer ii,ll
			real*8 sum
c-----------------------------------------------------------------------
c
      ii = 0
      do i = 1,n
        ll = indx(i)
        sum = b(ll)
        b(ll) = b(i)
        if ( ii /= 0 ) then
          do j = ii,i-1
            sum = sum - a(i,j)*b(j)
          end do
        else if ( abs(sum) >= eps ) then
          ii = i
        end if
        b(i) = sum
      end do
      do i = n,1,-1
        sum = b(i)
        do j = i+1,n
          sum = sum - a(i,j)*b(j)
        end do
        b(i) = sum / a(i,i)
      end do
c
      return
      end subroutine ummdp_utility_lubksb
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     CALCULATE INVERSE MATRIX 2x2 
c
      subroutine ummdp_utility_minv2 ( b,a,deta,eps )
c
c-----------------------------------------------------------------------
      implicit none
c
			real*8,intent(in) :: eps
			real*8,intent(in) :: a(2,2)
c
      real*8,intent(out) :: deta
      real*8,intent(out) :: b(2,2)
c
			real*8 detai
c-----------------------------------------------------------------------
c
      deta = a(1,1)*a(2,2) - a(1,2)*a(2,1)
      if ( abs(deta) <= eps ) then
         write (6,*) 'determinant det[a] error',deta
         write (6,*) 'stop in minv2'
         call ummdp_exit ( 9000 )
      end if
c
      detai = 1.0d0 / deta
      b(1,1) =  a(2,2) * detai
      b(1,2) = -a(1,2) * detai
      b(2,1) = -a(2,1) * detai
      b(2,2) =  a(1,1) * detai
c
      return
      end subroutine ummdp_utility_minv2
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CALCULATE INVERSE MATRIX 3x3
c
      subroutine ummdp_utility_minv3 ( b,a,deta,eps )
c-----------------------------------------------------------------------
      implicit none
c
			real*8,intent(in) :: eps
      real*8,intent(in) :: a(3,3)
c
      real*8,intent(out) :: deta 
      real*8,intent(out) :: b(3,3)
c
			real*8 detai
c-----------------------------------------------------------------------
c
      deta = a(1,1) * (a(2,2)*a(3,3) - a(2,3)*a(3,2))
     1       + a(1,2) * (a(2,3)*a(3,1) - a(2,1)*a(3,3))
     2       + a(1,3) * (a(2,1)*a(3,2) - a(2,2)*a(3,1))
      if ( abs(deta) <= eps ) then
         write (6,*) 'determinant det[a] error',deta
         write (6,*) 'stop in minv3'
         call ummdp_exit ( 9000 )
      end if
c
      detai = 1.0d0 / deta
      b(1,1) = (a(2,2)*a(3,3) - a(2,3)*a(3,2)) * detai
      b(1,2) = (a(1,3)*a(3,2) - a(1,2)*a(3,3)) * detai
      b(1,3) = (a(1,2)*a(2,3) - a(1,3)*a(2,2)) * detai
      b(2,1) = (a(2,3)*a(3,1) - a(2,1)*a(3,3)) * detai
      b(2,2) = (a(1,1)*a(3,3) - a(1,3)*a(3,1)) * detai
      b(2,3) = (a(1,3)*a(2,1) - a(1,1)*a(2,3)) * detai
      b(3,1) = (a(2,1)*a(3,2) - a(2,2)*a(3,1)) * detai
      b(3,2) = (a(1,2)*a(3,1) - a(1,1)*a(3,2)) * detai
      b(3,3) = (a(1,1)*a(2,2) - a(1,2)*a(2,1)) * detai
c
      return
      end subroutine ummdp_utility_minv3
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c     CALCULATE EIGENVALUES AND EIGENVECTORS BY JACOBI METHOD
c
c     Ref.: http://www.flagshyp.com/
c
c     input
c       a(3,3)  : symmetric matrix to be analyzed
c
c     output
c       es(i)   : i-th eigenvalue
c       ev(i,3) : normalized eigenvector for i-th eigenvalue
c
      subroutine ummdp_utility_eigen_sym3 ( es,ev,a )
c-----------------------------------------------------------------------
      implicit none
c
      real*8 es(3)
			real*8 ev(3,3),a(3,3)
c
			integer i,j,is,ip,iq,ir
			integer msweep
			real*8 eps,ax,er,sum,od,hd,theta,t,c,s,tau,h,g
      real*8 w(3,3),prc(3,3)
c-----------------------------------------------------------------------
c
      msweep = 100
      eps = 1.0d-8
c
c                                                       ---- preparation
      ax = 0.0d0
      er = 0.0d0
      do i = 1,3
        do j = 1,3
          if ( abs(a(i,j)) > ax ) ax = abs(a(i,j))
          er = er + abs(a(i,j)-a(j,i))
        end do
      end do
      if ( er/ax > eps ) then
        write (6,*) 'a is not symmetric'
        write (6,*) 'stop in ummdp_eigen_sym3'
        call ummdp_exit ( 9000 )
      end if
      do i = 1,3
        do j = 1,3
          w(i,j) = a(i,j) / ax
        end do
      end do
c                                    ---- initialise prc to the identity
      do i = 1,3
        do j = 1,3
          prc(i,j) = 0.0d0
        end do
        prc(i,i) = 1.0d0
        es(i) = w(i,i)
      end do
c                                                   ---- starts sweeping
      do is = 1,msweep
c
        sum = 0.0d0
        do ip = 1,2
          do iq = ip+1,3
            sum = sum + abs( w(ip,iq) )
          end do
        end do
c       write (6,*) 'ite',is,sum,eps
c            ---- if the sum of off-diagonal terms is zero evaluates the
c                                                     esches and returns
        if ( abs(sum) < eps ) then
          do i = 1,3
            do j = 1,3
              ev(i,j) = prc(j,i)
            end do
            es(i) = es(i)*ax
          end do
          return
        end if
c                             ---- performs the sweep in three rotations
c                                         ---- one per off diagonal term
        do ip = 1,2
          do iq = ip+1,3
            od = 100.0d0 * abs( w(ip,iq) )
            if ( abs(od) > eps ) then
              hd = es(iq) - es(ip)
c                                      ---- evaluates the rotation angle
              theta = 0.5d0 * hd / w(ip,iq)
              t = 1.0d0/(abs(theta) + sqrt(1.0d0+theta**2))
              if ( theta < 0.0d0 ) t = -t
c                                   ---- re-evaluates the diagonal terms
              c = 1.0d0 / sqrt(1.0d0+t**2)
              s = t * c
              tau = s / (1.0d0+c)
              h = t * w(ip,iq)
              es(ip) = es(ip) - h
              es(iq) = es(iq) + h
c                     ---- re-evaluates the remaining off-diagonal terms
              ir = 6 - ip - iq
              g = w( min(ir,ip),max(ir,ip) )
              h = w( min(ir,iq),max(ir,iq) )
              w( min(ir,ip),max(ir,ip) ) = g - s*(h+g*tau)
              w( min(ir,iq),max(ir,iq) ) = h + s*(g-h*tau)
c                                          ---- rotates the eigenvectors
              do ir = 1,3
                g = prc(ir,ip)
                h = prc(ir,iq)
                prc(ir,ip) = g - s*(h+g*tau)
                prc(ir,iq) = h + s*(g-h*tau)
              end do
            end if
            w(ip,iq) = 0.0d0
          end do
        end do
      end do
c                              ---- if convergence is not achieved stops
      write (6,*) 'did not converge in eigen calculation.'
      write (6,*) 'msweep=',msweep
      write (6,*) 'eps=',eps
      write (6,*) 'sum=',sum
      write (6,*) 'stop in ummdp_eigen_sym3'
      call ummdp_exit ( 9000 )
c
      return
      end subroutine ummdp_utility_eigen_sym3
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     CHECKING EXISTENCE OF FILE NAMES 'FLNAME'
c
      logical function ummdp_utility_file_exist ( flname )
c
c-----------------------------------------------------------------------
      implicit none
c
      character*16,intent(in) :: flname
c
			integer nio
c-----------------------------------------------------------------------
c
      nio = 616
      open  ( nio,file=flname,status='old',err=10 )
c
      close ( nio,status='keep' )
      ummdp_utility_file_exist = .true.
      return
c
   10 ummdp_utility_file_exist = .false.
      return
c
      end function ummdp_utility_file_exist
c
c
c
************************************************************************
*
*     YIELD FUNCTIONS
*
************************************************************************
c
c      0 : von Mises (1913)
c
c     3D
c       1 : Hill 1948
c       2 : Yld2004-18p
c       3 : CPB 2006
c       4 : Karafillis-Boyce 1993
c       5 : Hu 2005
c       6 : Yoshida 2011
c
c     2D
c      -1 : Gotoh Biquadratic
c      -2 : Yld2000-2d
c      -3 : Vegter
c      -4 : BBC 2005
c      -5 : Yld89
c      -6 : BBC 2008
c      -7 : Hill 1990
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     CALCULATE YIELD FUNCTION AND DERIVATIVES
c
      subroutine ummdp_yield ( se,cdseds,cd2seds2,nreq,cs,nttl,nnrm,
     1                         nshr,pryld,ndyld )                 
c
c-----------------------------------------------------------------------
      implicit none
c
      integer nreq,nttl,nnrm,nshr,ndyld
      real*8 se
			real*8 cdseds(nttl),cs(nttl),pryld(ndyld)
			real*8 cd2seds2(nttl,nttl)
c       
      integer i,j
      integer ntyld
      integer indx(6)
			real*8 ss
	    real*8 s(6),dseds(6)
			real*8 d2seds2(6,6)
c-----------------------------------------------------------------------
c
      ntyld = nint(pryld(1))
c
      if ( ntyld < 0 ) then
        if ( ( nnrm /= 2 ) .or. ( nshr /= 1 ) ) then
          write (6,*) 'error in ummdp_yield'
          write (6,*) 'ntyld<0 for plane stress'
          write (6,*) 'nnrm,nshr,ntyld:',nnrm,nshr,ntyld
          call ummdp_exit (9000)
        end if
        goto 100
      end if
c
      ss = 0.0
      do i = 1,nttl
        ss = ss + cs(i)**2
      end do
      if ( (ss <= 0.0) .and. (nreq == 0) ) then
        se = 0.0
        return
      end if
c
c                                                ---- 3D yield functions
c
c                                        ---- set index to s(i) to cs(i)
      do i = 1,6
        indx(i) = 0
      end do
      if ( nnrm == 3 ) then
        do i = 1,nttl
          indx(i) = i
        end do
      else if ( nnrm == 2 ) then
        indx(1) = 1
        indx(2) = 2
        indx(3) = 0
        do i = 1,nshr
          indx(3+i) = 2 + i
        end do
      end if
c                                                          ---- set s(i)
      call ummdp_utility_clear1 ( s,6 )
      do i = 1,6
        if ( indx(i) /= 0 ) then
          s(i) = cs(indx(i))
        end if
      end do
c
      select case ( ntyld )
      case ( 0 )                                             ! von Mises
        call ummdp_mises ( s,se,dseds,d2seds2,nreq )
c
      case ( 1 )                                             ! Hill 1948
        call ummdp_hill1948 ( s,se,dseds,d2seds2,nreq,pryld,ndyld ) 
c
      case ( 2 )                                           ! Yld2004-18p
        call ummdp_yld2004_18p ( s,se,dseds,d2seds2,nreq,pryld,ndyld )                         
c
      case ( 3 )                                              ! CPB 2006
        call ummdp_cpb2006 ( s,se,dseds,d2seds2,nreq,pryld,ndyld )                      
c
      case ( 4 )                                 ! Karafillis-Boyce 1993
        call ummdp_karafillis_boyce ( s,se,dseds,d2seds2,nreq,
     1                                pryld,ndyld )                         
c
      case ( 5 )                                               ! Hu 2005
        call ummdp_hu2005 ( s,se,dseds,d2seds2,nreq,pryld,ndyld )        
c
      case ( 6 )                                          ! Yoshida 2011
        call ummdp_yoshida2011 ( s,se,dseds,d2seds2,nreq,pryld,ndyld )            
c
      case default
        write (6,*) 'error in ummdp_yield'
        write (6,*) 'ntyld error :',ntyld
        call ummdp_exit (9000)
      end select
c
c                                                        ---- set dse/ds
      if ( nreq >= 1 ) then
        do i = 1,6
          if ( indx(i) /= 0 ) cdseds(indx(i)) = dseds(i)
        end do
      end if
c                                                      ---- set d2se/ds2
      if ( nreq >= 2 ) then
        do i = 1,6
          if ( indx(i) /= 0 ) then
            do j = 1,6
              if ( indx(j) /= 0 ) then
                cd2seds2(indx(i),indx(j)) = d2seds2(i,j)
              end if
            end do
          end if
        end do
      end if
c
      return
c
c
  100 continue
c                                       ---- plane stress yield criteria
c
      select case ( ntyld )
      case ( -1 )                                    ! Gotoh Biquadratic
        call ummdp_gotoh ( cs,se,cdseds,cd2seds2,nreq,pryld,ndyld )
c
      case ( -2 )                                           ! Yld2000-2d
        call ummdp_yld2000 ( cs,se,cdseds,cd2seds2,nreq,pryld,ndyld )
c
      case ( -3 )                                               ! Vegter
        call ummdp_vegter ( cs,se,cdseds,cd2seds2,nreq,pryld,ndyld )
c
      case ( -4 )                                             ! BBC 2005
        call ummdp_bbc2005 ( cs,se,cdseds,cd2seds2,nreq,pryld,ndyld )  
c
      case ( -5 )                                                ! Yld89
        call ummdp_yld89 ( cs,se,cdseds,cd2seds2,nreq, pryld,ndyld )                 
c
      case ( -6 )                                             ! BBC 2008
        call ummdp_bbc2008 ( cs,se,cdseds,cd2seds2,nreq,pryld,ndyld )               
c
      case ( -7 )                                            ! Hill 1990
        call ummdp_hill1990 ( cs,se,cdseds,cd2seds2,nreq,pryld,ndyld )                 
c
      case default
        write (6,*) 'error in ummdp_yield'
        write (6,*) 'ntyld error :',ntyld
        call ummdp_exit (9000)
      end select
c
      return
      end subroutine ummdp_yield
c
c
c************************************************************************
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
      ! implicit none
c
      integer,intent(in) :: n
c
      ! integer i
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
c************************************************************************
c     BBC2008 YIELD FUNCTION AND DERIVATIVES
c
c       doi: 
c
      subroutine ummdp_bbc2008 ( s,se,dseds,d2seds2,nreq,pryld,ndyld )
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
      call ummdp_bbc2008_core ( s,se,dseds,d2seds2,nreq,
     1                          pryld,ndyld,nds,ndk )
c
      return
      end subroutine ummdp_bbc2008
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      subroutine ummdp_bbc2008_core ( s,se,dseds,d2seds2,nreq,
     1                                pryld,ndyld,sp,kp )
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
     1       ummdp_bbc2008_get_se 
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
      call ummdp_bbc2008_setup ( pryld,ndyld,sp,kp,wp,Lp,Mp,Np,kCm,se1 )
c ----------------
c  se section
c ----------------
c
c                             ---- The unit of this se is (stress)^(2kp)
      se = ummdp_bbc2008_get_se ( sp,kp,wp,s,Lp,Mp,Np,kCm )
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
      call ummdp_utility_clear1 ( dFds,3 )
      call ummdp_utility_clear2 ( d2Fds2,3,3 )

c                              ---- long-long-long do loops starts here.
      do csp = 1,sp
c
        call ummdp_bbc2008_get_w_phi ( wpi,phiL,phiM,phiN,
     1                                 csp,sp,wp,Lp,Mp,Np,s )
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
          call ummdp_bbc2008_get_dphiXds ( dphiLds,Lp,s,csp,m,sp )
          call ummdp_bbc2008_get_dphiXds ( dphiMds,Mp,s,csp,(kp-m),sp )
          call ummdp_bbc2008_get_dphiXds ( dphiNds,Np,s,csp,m,sp )
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
            call ummdp_bbc2008_get_d2phiXds2 (d2phiLds2,Lp,s,csp,m,sp)
            call ummdp_bbc2008_get_d2phiXds2 (d2phiMds2,Mp,s,csp,
     1                                                       (kp-m),sp)
            call ummdp_bbc2008_get_d2phiXds2 (d2phiNds2,Np,s,csp,m,sp)
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
      end subroutine ummdp_bbc2008_core
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     ummdp_bbc2008_get_w_phi ()
c     A subroutine to get w^(i-1), w^(s-i) and phiX variables
c
      subroutine ummdp_bbc2008_get_w_phi ( wpi,phiL,phiM,phiN,csp,sp,
     1                                     wp,Lp,Mp,Np,s )
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
      real*8 ummdp_bbc2008_get_phiX
c-----------------------------------------------------------------------
c
      wpi(1) = wp**(csp-1)
      wpi(2) = wp**(sp-csp)
      phiL = ummdp_bbc2008_get_phiX (Lp, s, csp , sp)
      phiM = ummdp_bbc2008_get_phiX (Mp, s, csp , sp)
      phiN = ummdp_bbc2008_get_phiX (Np, s, csp , sp)
c
      return
      end subroutine ummdp_bbc2008_get_w_phi
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CALCULATE EQUIVALENT STRESS
c
      real*8 function ummdp_bbc2008_get_se ( sp,kp,wp,s,Lp,Mp,Np,kCm )
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
      ummdp_bbc2008_get_se = 0.0d0
c
      do csp = 1,sp
c
        call ummdp_bbc2008_get_w_phi ( wpi,phiL,phiM,phiN,csp,sp,wp,
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
          ummdp_bbc2008_get_se = 
     1    ummdp_bbc2008_get_se 
     2     + kCm(m) * phiM_kp_m * ( wpi(1) * phiL_m + wpi(2) * phiN_m )
c
        end do
c
      end do
c
      return
      end function ummdp_bbc2008_get_se
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     ummdp_bbc2008_get_phiX (Xp, s, csp , sp)
c     A function to calculate s(a)*X(a,b)*s(b) (summation convention)
c
      real*8 function ummdp_bbc2008_get_phiX ( Xp,s,csp,sp )
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
      call ummdp_utility_mv ( v,XXp,s,nc,nc)
      call ummdp_utility_vvs (ummdp_bbc2008_get_phiX,v,s,nc)
c
      return
      end function ummdp_bbc2008_get_phiX
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     ummdp_bbc2008_get_dphiXds (dphiXds, Xp, s, csp, lambda,sp)
c     A subroutine to calculate d(phiX^(lambda))/ds.
c     It returns dphiXds(nc).
c
      subroutine ummdp_bbc2008_get_dphiXds ( dphiXds,Xp,s,csp,lambda,sp)
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
c
      call ummdp_utility_clear1(dphiXds,nc)
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
c       from ummdp_mv().
c
      call ummdp_utility_mv (v, XXp, s, nc, nc)
c
      if ( lambda == 1 ) then
        dphiXds(1:nc) = 2.0d0 * v(1:nc)
      else
        call ummdp_utility_vvs (phi, v, s, nc)
        dphiXds(1:nc) = 2.0d0 * lambda * phi**(lambda-1) * v(1:nc)
      end if
c
      return
      end subroutine ummdp_bbc2008_get_dphiXds
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     ummdp_bbc2008_get_d2phiXds2 (d2phiXds2, Xp, s, csp, lambda,sp)
c     A subroutine to calculate d2(phiX^(lambda))/(dsds').
c     It returns d2phiXdsds(nc,nc).
c
      subroutine ummdp_bbc2008_get_d2phiXds2 ( d2phiXds2,Xp,s,csp,
     1                                         lambda,sp )
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
      call ummdp_utility_mv (v, XXp, s, nc, nc)
      call ummdp_utility_vvs (phi, v, s, nc)
c
      phi_lambda2 = 1.0d0
      if ( lambda /= 2 ) then
        phi_lambda2 = phi**(lambda-2)
      end if
c
      call ummdp_utility_clear2 ( d2phiXds2,nc,nc )
c
c                                            ---- d2phiX/(ds(i)ds(1:nc))
      do i = 1,nc
        d2phiXds2(i, 1:nc) = 2.0d0 * lambda * phi_lambda2 * 
     1   ( 2.0d0 * (lambda - 1) * v(1:nc) * v(i)
     2    + phi * XXp(1:nc, i) )
      end do
c
      return
      end subroutine ummdp_bbc2008_get_d2phiXds2
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     setup_bbc2008_parameters()
c     A routine to setup local variables.
c
      subroutine ummdp_bbc2008_setup ( pryld,ndyld,sp,kp,wp,
     1                                 Lp,Mp,Np,kCm,se1 )
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
      real*8 ummdp_bbc2008_get_se
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
        call ummdp_bbc2008_setup_MN_tensors (m,csp,pryld,ndyld,Mp,sp)
c                                                      ---- N^(i) tensor
        n = m + 3
        call ummdp_bbc2008_setup_MN_tensors (n,csp,pryld,ndyld,Np,sp)
      end do
c
c
c     equiv. stress in uniaxial stress state.
c     dummy_s = (1.0d0, 0.0d0, 0.0d0)
c     ** The unit of this se1 is (stress)^(2kp)
c
      call ummdp_utility_clear1 (dummy_s, 3)
      dummy_s(1) = 1.0d0
      se1 = ummdp_bbc2008_get_se (sp, kp, wp, dummy_s, Lp, Mp, Np, kCm)
c
      return
      end subroutine ummdp_bbc2008_setup
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
      subroutine ummdp_bbc2008_setup_MN_tensors ( ic,csp,
     1                                            pryld,ndyld,Xp,sp )
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
      end subroutine ummdp_bbc2008_setup_MN_tensors
c
c
c************************************************************************
c     CPB2006 YIELD FUNCTION AND DERIVATIVES
c
c       doi: 
c
c     !!! CAUTION !!!
c     Plane stress condition is NOT implemented in this code.
c
      subroutine ummdp_cpb2006 ( s,se,dseds,d2seds2,nreq,pryld,ndyld )
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: nreq,ndyld
      real*8 ,intent(in) :: s(6),pryld(ndyld)
c
      real*8,intent(out) :: se
      real*8,intent(out) :: dseds(6)
			real*8,intent(out) :: d2seds2(6,6)
c
      integer i,j,k,m,n,l,iq,ip,ir
      real*8 pi,eps,a,ck,ai,H1,H2,H3,p,q,theta,F,D,DseDF,denom,D2seDF2,
     1       del,sea,seb,abc1,abc2,seaa,seba,seab,sebb,
     2       ummdp_cpb2006_seND
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
      call ummdp_utility_clear2 ( c,6,6 )
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
      call ummdp_utility_clear2 ( ct,6,6 )
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
      call ummdp_utility_mv ( sigma,ct,s,6,6 )
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
        call ummdp_utility_clear1 ( DFDH,3 )
        call ummdp_utility_clear2 ( DpsigmaDH,3,3 )
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
        call ummdp_utility_clear2 ( DHDsigma,3,6 )
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
        call ummdp_utility_clear1 ( DFDs,6 )
        call ummdp_utility_clear2 ( dummat,3,6 )
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
        call ummdp_utility_clear2 ( D2FDpsigma2,3,3 )
c
        do i = 1,3
          D2FDpsigma2(i,i) = a*(psigma(i)/abs(psigma(i))-ck)**2.0d0 *
     1                       (abs(psigma(i))-ck*psigma(i))**(a-2.0d0)
        end do
c
c                              ---- D2(psigma)/D(H)2 -> 3 x 3 x 3 Matrix
        call ummdp_utility_clear3 ( D2psigmaDH2,3,3,3 )
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
          call ummdp_utility_clear2 ( D2FDH2,3,3 )
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
          call ummdp_utility_clear3 ( D2HDsigma2,3,6,6 )
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
          call ummdp_utility_clear2 ( D2FDs2,6,6 )
          call ummdp_utility_clear2 ( dummat,3,6 )
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
                sea = ummdp_cpb2006_seND ( s0,ct,phi,ck,a,ai )
                s0(i) = s(i) + del
                seb = ummdp_cpb2006_seND ( s0,ct,phi,ck,a,ai )
c
                s0(i) = s(i)
                abc1 = (se-sea) / del
                abc2 = (seb-se) / del
                d2seds2(i,j) = (abc2-abc1) / del
              else
                s0(i) = s(i) - del
                s0(j) = s(j) - del
                seaa = ummdp_cpb2006_seND ( s0,ct,phi,ck,a,ai )
c
                s0(i) = s(i) + del
                s0(j) = s(j) - del
                seba = ummdp_cpb2006_seND ( s0,ct,phi,ck,a,ai )
c
                s0(i) = s(i) - del
                s0(j) = s(j) + del
                seab = ummdp_cpb2006_seND ( s0,ct,phi,ck,a,ai )
c
                s0(i) = s(i) + del
                s0(j) = s(j) + del
                sebb = ummdp_cpb2006_seND ( s0,ct,phi,ck,a,ai )
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
      end subroutine ummdp_cpb2006
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     NUMERICAL DIFFERENTIATION FOR EQUIVALENT STRESS
c
      real*8 function ummdp_cpb2006_seND ( s,ct,phi,ck,a,ai )
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
      call ummdp_utility_mv ( sigma,ct,s,6,6 )
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
      ummdp_cpb2006_seND = (F/D) ** ai
c
      return
      end function ummdp_cpb2006_seND
c
c
c
************************************************************************
c     GOTOH BIQUADRATIC YIELD FUNCTION AND DERIVATIVES
c
c       doi:
c
      subroutine ummdp_gotoh ( s,se,dseds,d2seds2,nreq,pryld,ndyld )          
c------------------------------------------------------------- variables
      implicit none
c
      integer,intent(in) :: nreq,ndyld
      real*8 ,intent(in) :: s(3),pryld(ndyld)
c
      real*8,intent(out) :: se
      real*8,intent(out) :: dseds(3)
			real*8,intent(out) :: d2seds2(3,3)
c
      integer i,j,k,m,n
      real*8 phi,q
      real*8 a(9),v(4),t(4)
      real*8 c(4,4),dtds(4,3)
      real*8 d2tds2(4,3,3)
c-----------------------------------------------------------------------
c     a(i)      : coef.s of Gotoh's 4th order function
c     t(i)      : stress^2 vector
c                 ={ sx^2, sx*sy, sy^2, txy^2 }
c     c(i,j)    : matrix to calc. se
c              Gotoh's function se = ( {t}^T*[c]*{t} )^(1/4)
c
c     dtds(i,j) : diff. of t(i) with respect to s(j)
c     d2tds2(i,j,k)
c               : 2nd order diff. of t(i) w.r.t s(j) & s(k)
c-----------------------------------------------------------------------
c
c                                            ---- anisotropic parameters
      do i = 1,9
        a(i) = pryld(1+i)
      end do
c                                                      ---- coef. matrix
      c(1,1) = a(1)
      c(1,2) = a(2)*0.5d0
      c(1,3) = 0.0
      c(1,4) = a(6)*0.5d0
      c(2,2) = a(3)
      c(2,3) = a(4)*0.5d0
      c(2,4) = a(7)*0.5d0
      c(3,3) = a(5)
      c(3,4) = a(8)*0.5d0
      c(4,4) = a(9)
      do i = 2,4
        do j = 1,i-1
          c(i,j) = c(j,i)
        end do
      end do
c                                                    ---- t-vector (s^2)
      t(1) = s(1) * s(1)
      t(2) = s(1) * s(2)
      t(3) = s(2) * s(2)
      t(4) = s(3) * s(3)
c                                                 ---- equivalent stress
      call ummdp_utility_mv  ( v,c,t,4,4 )
      call ummdp_utility_vvs ( phi,t,v,4 )
c
      if ( phi <= 0.0d0 ) phi = 0.0d0
      se = sqrt(sqrt(phi))
c                                            ---- 1st order differential
      if ( nreq >= 1 ) then
        call ummdp_utility_clear2 ( dtds,4,3 )
        dtds(1,1) = s(1) * 2.0d0
        dtds(2,1) = s(2)
        dtds(2,2) = s(1)
        dtds(3,2) = s(2) * 2.0d0
        dtds(4,3) = s(3) * 2.0d0
        call ummdp_utility_clear1 ( v,4 )
        do i = 1,3
          do j = 1,4
            do k = 1,4
              v(i) = v(i) + 2.0d0*t(j)*c(j,k)*dtds(k,i)
            end do
          end do
        end do
        q = 0.25d0 * phi**(-0.75d0)
        do i = 1,3
          dseds(i) = q * v(i)
        end do
      end if
c                                            ---- 2nd order differential
      if ( nreq >= 2 ) then
        call ummdp_utility_clear3 ( d2tds2,4,3,3 )
        d2tds2(1,1,1) = 2.0d0
        d2tds2(2,1,2) = 1.0d0
        d2tds2(2,2,1) = 1.0d0
        d2tds2(3,2,2) = 2.0d0
        d2tds2(4,3,3) = 2.0d0
        call ummdp_utility_clear2 ( d2seds2,3,3 )
        do i = 1,3
          do j = 1,3
            do m = 1,4
              do n = 1,4
                d2seds2(i,j) = d2seds2(i,j)+
     1                       2.0d0*c(m          ,n       )*
     2                        ( dtds(m,i)*dtds(  n  ,j)+
     3                           t(  m)  *d2tds2(n,i,j)  )
              end do
            end do
          end do
        end do
        do i = 1,3
          do j = 1,3
            d2seds2(i,j) = q*(d2seds2(  i,   j)
     1                      -0.75d0*v(i)*v(j)/phi)
          end do
        end do
      end if
c
      return
      end subroutine ummdp_gotoh
c
c
c
c************************************************************************
c     Hill 1948 YIELD FUNCTION AND DERIVATIVES
c
c       doi:
c
      subroutine ummdp_hill1948 ( s,se,dseds,d2seds2,nreq,pryld,ndyld )
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: nreq,ndyld
      real*8 ,intent(in) :: s(6),pryld(ndyld)
c
      real*8,intent(out) :: se
      real*8,intent(out) :: dseds(6)
			real*8,intent(out) :: d2seds2(6,6)
c
      integer i,j
			real*8 pf,pg,ph,pl,pm,pn,phi
      real*8 v(6)
			real*8 c(6,6)
c-----------------------------------------------------------------------
c
c                                            ---- anisotropic parameters
      pf = pryld(1+1)
      pg = pryld(1+2)
      ph = pryld(1+3)
      pl = pryld(1+4)
      pm = pryld(1+5)
      pn = pryld(1+6)
c                                               ---- coefficients matrix
      call ummdp_utility_clear2 ( c,6,6 )
      c(1,1) = pg + ph
      c(1,2) = -ph
      c(1,3) = -pg
      c(2,1) = -ph
      c(2,2) = pf + ph
      c(2,3) = -pf
      c(3,1) = -pg
      c(3,2) = -pf
      c(3,3) = pf + pg
      c(4,4) = 2.0d0 * pn
      c(5,5) = 2.0d0 * pm
      c(6,6) = 2.0d0 * pl
      do i = 1,6
        do j = 1,6
          c(i,j) = c(i,j) / (pg+ph)
        end do
      end do
c
      call ummdp_utility_mv  ( v,c,s,6,6 )
      call ummdp_utility_vvs ( phi,s,v,6 )
c                                                 ---- equivalent stress
      if ( phi <= 0.0 ) phi = 0.0
      se = sqrt(phi)
c                                              ---- 1st order derivative
      if ( nreq >= 1 ) then
        do i = 1,6
          dseds(i) = v(i) / se
        end do
      end if
c                                              ---- 2nd order derivative
      if ( nreq >= 2 ) then
        do i = 1,6
          do j = 1,6
            d2seds2(i,j) = (-v(i)*v(j)/phi+c(i,j)) / se
          end do
        end do
      end if
c
      return
      end subroutine ummdp_hill1948
c
c
c
************************************************************************
c     HILL 1990 YIELD FUNCTION AND DERIVATIVES
c
c       doi: https://doi.org/10.1016/0022-5096(90)90006-P
c
      subroutine ummdp_hill1990 ( s,se,dseds,d2seds2,nreq,pryld,ndyld )                    
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
      call ummdp_utility_clear2( a4,3,3 )
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
      end subroutine ummdp_hill1990
c
c
c
************************************************************************
c     HU2005 YIELD FUNCTION AND DERIVATIVES
c
c       doi: https://doi.org/10.1016/j.ijplas.2004.11.004
c
      subroutine ummdp_hu2005 ( s,se,dseds,d2seds2,nreq,pryld,ndyld )
c-----------------------------------------------------------------------
      implicit none
c
      integer,parameter :: maxa = 100
c
      integer,intent(in) :: nreq,ndyld
      real*8 ,intent(in) :: s(6),pryld(ndyld)
c
      real*8,intent(out) :: se
      real*8,intent(out) :: dseds(6)
			real*8,intent(out) :: d2seds2(6,6)
c 
      integer nd0,n,nterms,it,jy,jx,ndmax
      real*8 a(maxa)
      real*8 ipow(maxa,3)
c-----------------------------------------------------------------------
c
      nd0 = 2
c
      n = 0
      do it = 0,nd0
        n = n + (nd0-it)*2 + 1
      end do
      nterms = n
      if ( maxa < nterms ) then
        write (6,*) 'increase maxa :',maxa,nterms
        call ummdp_exit ( 9000 )
      end if
c
      n = 0
      ipow = 0.0d0
      do it = 0,nd0
        ndmax = nd0*2 - it*2
        do jy = 0,ndmax
          jx = ndmax - jy
          n = n + 1
          ipow(n,1) = jx
          ipow(n,2) = jy
          ipow(n,3) = it
        end do
      end do
c
      a    =  0.0d0
      a(1) =  pryld(1+1)    ! X1 
      a(2) =  pryld(1+2)    ! X2
      a(3) =  pryld(1+3)    ! X3
      a(4) =  pryld(1+4)    ! X4
      a(5) =  pryld(1+5)    ! X5
      a(6) =  pryld(1+7)    ! C1 <-
      a(7) = -pryld(1+9)    ! C3 <-
      a(8) =  pryld(1+8)    ! C2 <-
      a(9) =  pryld(1+6)    ! X7
c
      call ummdp_hy_polytype ( s,se,dseds,d2seds2,nreq,nd0,
     1                          a,ipow,maxa,nterms )
c
      return
      end subroutine ummdp_hu2005
c
c
c************************************************************************
c     KARAFILLIS-BOYCE 1993 YIELD FUNCTION AND DERIVATIVES
c
c       doi:
c
      subroutine ummdp_karafillis_boyce ( s,se,dseds,d2seds2,nreq,
     1                                    pryld,ndyld )
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: nreq,ndyld
      real*8 ,intent(in) :: s(6),pryld(ndyld)
c
      real*8,intent(out) :: se
      real*8,intent(out) :: dseds(6)
			real*8,intent(out) :: d2seds2(6,6)
c
      integer k_p,i,j,k,m,n,eqFlag
      real*8 c_p,phi,DseDphi,X12,X13,X23,dum,tol
      real*8 smallS(6),largeS(3),Jinvar(3),phiN(2),coef(2),beta(3),
     1       alpha(2)
      real*8 L(6,6),DjDss(3,6),DphiDs(2,6),DphiDj(2,3),DphiDls(2,3),
     1       DlsDj(3,3),workmat(3,6)
      real*8 DDphiDDs(2,6,6),DDjDDss(3,6,6),DDphiDDj(2,3,3),
     1       DDphiDDls(2,3,3),DDlsDDj(3,3,3),workmat1(2,6,6)
c-----------------------------------------------------------------------
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
c-----------------------------------------------------------------------
c
      tol = 1.0e-7
c                                                        ---- parameters
      do i = 1,6
        do j = 1,6
          L(i,j) = 0.0d0
        end do
      end do
c
      alpha(1) = pryld(1+2)
      alpha(2) = pryld(1+3)
      beta(1)  = (alpha(2)-1.0d0   -alpha(1)) * 0.5d0
      beta(2)  = (alpha(1)-alpha(2)-1.0d0   ) * 0.5d0
      beta(3)  = (1.0d0   -alpha(1)-alpha(2)) * 0.5d0
      L(1,1)   = pryld(1+1) * 1.0d0
      L(1,2)   = pryld(1+1) * beta(1)
      L(1,3)   = pryld(1+1) * beta(2)
      L(2,2)   = pryld(1+1) * alpha(1)
      L(2,3)   = pryld(1+1) * beta(3)
      L(3,3)   = pryld(1+1) * alpha(2)
      L(4,4)   = pryld(1+1) * pryld(1+4)
      L(5,5)   = pryld(1+1) * pryld(1+5)
      L(6,6)   = pryld(1+1) * pryld(1+6)
      L(2,1)   = L(1,2)
      L(3,1)   = L(1,3)
      L(3,2)   = L(2,3)
      k_p      = pryld(1+7)
      c_p      = pryld(1+8)
c                                                 ---- equivalent stress
      call ummdp_utility_mv ( smallS,L,s,6,6 )
      call ummdp_karafillis_boyce_principal_stress( smallS,Jinvar,
     1                                              largeS )
c
      phiN(1) = (largeS(1)-largeS(2))**(2*k_p)
     1        + (largeS(2)-largeS(3))**(2*k_p)
     2        + (largeS(3)-largeS(1))**(2*k_p)
      phiN(2) = largeS(1)**(2*k_p) 
     1        + largeS(2)**(2*k_p)
     2        + largeS(3)**(2*k_p)
      coef(1) = (1d0-c_p)
      coef(2) = c_p * 3d0**(2*k_p)/(2d0**(2*k_p-1)+1d0)
      phi = coef(1)*phiN(1) + coef(2)*phiN(2)
      se = (0.5d0*phi) ** (0.5d0/k_p)
c
c                                            ---- 1st order differential
      if ( nreq >= 1 ) then
c                 ---- check if there are two components of largeS equal 
c                                                          to each other
c               ---- if so, rearrange largeS so that largeS(1)=largeS(2)
        eqFlag = 0
        if ( abs(largeS(1)-largeS(2)) <= tol ) then
          eqFlag = 1
        else if ( abs(largeS(2)-largeS(3)) <= tol ) then
          eqFlag = 2
          dum = largeS(3)
          largeS(3) = largeS(1)
          largeS(1) = dum
        else if ( abs(largeS(1)-largeS(3)) <= tol ) then
          eqFlag = 3
          dum = largeS(3)
          largeS(3) = largeS(2)
          largeS(2) = dum
        end if

        if ( eqFlag == 0 ) then
          DphiDls(1,1) = 2d0 * k_p * ((largeS(1)-largeS(2))**(2*k_p-1)
     1       + (largeS(1)-largeS(3))**(2*k_p-1))
          DphiDls(1,2) = 2d0 * k_p * ((largeS(2)-largeS(1))**(2*k_p-1)
     1       + (largeS(2)-largeS(3))**(2*k_p-1))
          DphiDls(1,3) = 2d0 * k_p * ((largeS(3)-largeS(1))**(2*k_p-1)
     1       + (largeS(3)-largeS(2))**(2*k_p-1))
          do i = 1,3
            DphiDls(2,i) = 2d0 * k_p * largeS(i)**(2*k_p-1)
          end do
c 
          do i = 1,3
            dum = 1d0 /(3d0*largeS(i)**2 - 2d0*Jinvar(1)*largeS(i)
     1         + Jinvar(2))
            DlsDj(i,1) = largeS(i)**2 * dum
            DlsDj(i,2) = -largeS(i) * dum
            DlsDj(i,3) = dum
          end do
c
          do i = 1,2
            do j = 1,3
              DphiDj(i,j) = 0d0
              do k = 1,3
                DphiDj(i,j) = DphiDj(i,j) + DphiDls(i,k) * DlsDj(k,j)
              end do
            end do
          end do
        else
          if ( k_p == 1 ) then
            DphiDj(1,1) = 4d0*(2d0*largeS(1)+largeS(3))
            DphiDj(1,2) = -6d0
            DphiDj(1,3) = 0d0
          else
            dum = (largeS(1)-largeS(3))**(2*k_p-3)
            DphiDj(1,1) = 4d0*k_p*dum
     1         *(k_p*largeS(1)**2-largeS(1)*largeS(3)-largeS(3)**2)
            DphiDj(1,2) = -2d0*k_p*dum
     1         *((2*k_p-1)*largeS(1)-3d0*largeS(3))
            DphiDj(1,3) = 4d0*k_p*(k_p-2)*dum
          end if
c
          if ( abs(largeS(1)-largeS(3)) <= TOL ) then
            DphiDj(2,1) = 2d0*k_p**2*(2*k_p+1)*largeS(1)**(2*k_p-1)
            DphiDj(2,2) = -2d0*k_p**2*(2*k_p-1)*largeS(1)**(2*k_p-2)
            DphiDj(2,3) = 2d0*k_p*(2*k_p-1)*(k_p-1)
     1         *largeS(1)**(2*k_p-3)
          else
            dum = 2d0*k_p/(largeS(1)-largeS(3))**2
            DphiDj(2,1) = dum * ((2*k_p+1)*largeS(1)**(2*k_p)*
     1         (largeS(1)-largeS(3))
     2         - largeS(1)**(2*k_p+1) + largeS(3)**(2*k_p+1))
            DphiDj(2,2) = -dum * (2*k_p*largeS(1)**(2*k_p-1)*
     1         (largeS(1)-largeS(3))
     2         - largeS(1)**(2*k_p) + largeS(3)**(2*k_p))
            DphiDj(2,3) = dum * ((2*k_p-1)*largeS(1)**(2*k_p-2)*
     1         (largeS(1)-largeS(3))
     2         - largeS(1)**(2*k_p-1) + largeS(3)**(2*k_p-1))
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
        do i = 1,3
          do j = 1,6
            workmat(i,j) = 0d0
            do k = 1,6
              workmat(i,j) = workmat(i,j) + DjDss(i,k) * L(k,j)
            end do
          end do
        end do
        do i = 1,2
          do j = 1,6
            DphiDs(i,j) = 0d0
            do k = 1,3
              DphiDs(i,j) = DphiDs(i,j) + DphiDj(i,k) * workmat(k,j)
            end do
          end do
        end do
c
        DseDphi = se / (2d0*k_p*phi)
c
        do i = 1,6
          DphiDs(1,i) = coef(1)*DphiDs(1,i)+coef(2)*DphiDs(2,i)
          dseds(i) = DseDphi * DphiDs(1,i)
        end do
      end if
c                                            ---- 2nd order differential
      if ( nreq >= 2 ) then
        if ( eqFlag == 0) then
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
          do i = 1,3
            do j = i+1,3
              DDphiDDls(1,j,i) = DDphiDDls(1,i,j)
            end do
          end do
c
          do i = 1,3
            do j = 1,3
              DDphiDDls(2,i,j) = 0d0
            end do
            DDphiDDls(2,i,i) = 2*k_p*(2*k_p-1)*largeS(i)**(2*k_p-2)
          end do
c
          do i = 1,3
            dum = 1d0 /(3d0*largeS(i)**2 - 2d0*Jinvar(1)*largeS(i)
     1         + Jinvar(2))**3
            DDlsDDj(i,1,1) = dum*largeS(i)**3*
     1         (6d0*largeS(i)**2-6d0*Jinvar(1)*largeS(i)+4d0*Jinvar(2))
            DDlsDDj(i,1,2) = dum*largeS(i)**2*
     1         (-3d0*largeS(i)**2+4d0*Jinvar(1)*largeS(i)-3d0*Jinvar(2))
            DDlsDDj(i,1,3) = 2d0*dum*(-Jinvar(1)*largeS(i)**2
     1         +Jinvar(2)*largeS(i))
            DDlsDDj(i,2,2) = DDlsDDj(i,1,3)
            DDlsDDj(i,2,3) = dum*(3d0*largeS(i)**2-Jinvar(2))
            DDlsDDj(i,3,3) = -dum*(6d0*largeS(i)-2d0*Jinvar(1))
            DDlsDDj(i,2,1) = DDlsDDj(i,1,2)
            DDlsDDj(i,3,1) = DDlsDDj(i,1,3)
            DDlsDDj(i,3,2) = DDlsDDj(i,2,3)
          end do
c
          do i = 1,2
            do j = 1,3
              do k = j,3
                DDphiDDj(i,j,k) = 0d0
                do m = 1,3
                  do n = 1,3
                    DDphiDDj(i,j,k) = DDphiDDj(i,j,k)
     1                 + DDphiDDls(i,m,n)*DlsDj(m,j)*DlsDj(n,k)
                  end do
                  DDphiDDj(i,j,k) = DDphiDDj(i,j,k)
     1               + DphiDls(i,m)*DDlsDDj(m,j,k)
                end do
              end do
            end do
            do j = 1,3
              do k = j+1,3
                DDphiDDj(i,k,j) = DDphiDDj(i,j,k)
              end do
            end do
          end do
        else
          if ( k_p == 1 ) then
            do i = 1,3
              do j = 1,3
                DDphiDDj(1,i,j) = 0d0
              end do
            end do
            DDphiDDj(1,1,1) = 4d0
          else if ( k_p == 2 ) then
            do i = 1,3
              do j = 1,3
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
     1         -4*k_p*(2*k_p-1)*largeS(1)**3*largeS(3)
     2         -4*(2*k_p-1)*(k_p-2)*largeS(1)**2*largeS(3)**2
     3         +8*(k_p-2)*largeS(1)*largeS(3)**3
     4         +2*(2*k_p-1)*largeS(3)**4) * dum
            DDphiDDj(1,1,2) = (-2*k_p*(2*k_p-1)*(k_p-1)*largeS(1)**3/3d0
     1         +(2*k_p-1)*(5*k_p-4)*largeS(1)**2*largeS(3)
     2         +4*(k_p-2)**2*largeS(1)**1*largeS(3)**2
     3         -6*(k_p-1)*largeS(3)**3) * dum
            DDphiDDj(1,1,3) = dum *
     1         ((2*k_p-1)*(2*k_p**2-11*k_p+6)*largeS(1)**2/3d0
     2         -2*(k_p-2)*(2*k_p-3)*largeS(3)*(largeS(1)+largeS(3)))
            DDphiDDj(1,2,2) = dum *
     1         ((k_p-1)*(2*k_p-1)*(2*k_p-3)*largeS(1)**2/3d0
     2         -(12*k_p**2-22*k_p+14)*largeS(1)*largeS(3)
     3         +(10*k_p-11)*largeS(3)**2)
            DDphiDDj(1,2,3) = dum *
     1         ((-4*k_p**3+30*k_p**2-44*k_p+24)*largeS(1)/3d0
     2         +3*(k_p-2)*(2*k_p-3)*largeS(3))
            DDphiDDj(1,3,3) = dum * (4*k_p**3-48*k_p**2+107*k_p-78)/3d0
            do i = 1,3
              do j = i+1,3
                DDphiDDj(1,j,i) = DDphiDDj(1,i,j)
              end do
            end do
          end if
c
          do i = 1,3
            do j = i,3
              DDphiDDj(2,i,j) = 0d0
            end do
          end do
          do i = 0,2*k_p-2
            DDphiDDj(2,1,1) = DDphiDDj(2,1,1) +
     1         k_p*(i+1)*(i+2)*(i+3)*(2*k_p-i-1)/3*
     2         largeS(1)**i * largeS(3)**(2*k_p-2-i)
          end do
          do i = 0,2*k_p-3
            DDphiDDj(2,1,2) = DDphiDDj(2,1,2) -
     1         k_p*(i+1)*(i+2)*(i+3)*(2*k_p-i-2)/3*
     2         largeS(1)**i * largeS(3)**(2*k_p-3-i)
          end do
          do i = 0,2*k_p-4
            DDphiDDj(2,1,3) = DDphiDDj(2,1,3) +
     1         k_p*(i+1)*(i+2)*(i+3)*(2*k_p-i-3)/3*
     2         largeS(1)**i * largeS(3)**(2*k_p-4-i)
          end do
          DDphiDDj(2,2,2) = DDphiDDj(2,1,3)
          do i = 0,2*k_p-5
            DDphiDDj(2,2,3) = DDphiDDj(2,2,3) -
     1         k_p*(i+1)*(i+2)*(i+3)*(2*k_p-i-4)/3*
     2         largeS(1)**i * largeS(3)**(2*k_p-5-i)
          end do
          do i = 0,2*k_p-6
            DDphiDDj(2,3,3) = DDphiDDj(2,3,3) +
     1         k_p*(i+1)*(i+2)*(i+3)*(2*k_p-i-5)/3*
     2         largeS(1)**i * largeS(3)**(2*k_p-6-i)
          end do
          do i = 1,3
            do j = i+1,3
              DDphiDDj(2,j,i) = DDphiDDj(2,i,j)
            end do
          end do
        end if
c
        do i = 1,3
          do j = 1,6
            do k = 1,6
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
        do i = 1,6
          do j = i+1,6
            DDjDDss(2,j,i) = DDjDDss(2,i,j)
            DDjDDss(3,j,i) = DDjDDss(3,i,j)
          end do
        end do
c
        do i = 1,2
          do j = 1,6
            do k = 1,6
              workmat1(i,j,k) = 0d0
              do m = 1,3
                do n = 1,3
                  workmat1(i,j,k) = workmat1(i,j,k)
     1               + DDphiDDj(i,m,n)*DjDss(m,j)*DjDss(n,k)
                end do
                workmat1(i,j,k) = workmat1(i,j,k)
     1             + DphiDj(i,m)*DDjDDss(m,j,k)
              end do
            end do
          end do
        end do
c
        do i = 1,2
          do j = 1,6
            do k = 1,6
              DDphiDDs(i,j,k) = 0d0
              do m = 1,6
                do n = 1,6
                  DDphiDDs(i,j,k) = DDphiDDs(i,j,k)
     1               + workmat1(i,m,n)*L(m,j)*L(n,k)
                end do
              end do
            end do
          end do
        end do
c
        do j = 1,6
          do k = 1,6
            DDphiDDs(1,j,k) = coef(1)*DDphiDDs(1,j,k)
     1         + coef(2)*DDphiDDs(2,j,k)
          end do
        end do
c
        do i = 1,6
          do j = i,6
            d2seds2(i,j) = (1-2*k_p)*se/(4*k_p**2*phi**2)
     1         *DphiDs(1,i)*DphiDs(1,j) + se/(2*k_p*phi)*DDphiDDs(1,i,j)
          end do
        end do
        do i = 1,6
          do j = i+1,6
            d2seds2(j,i) = d2seds2(i,j)
          end do
        end do
      end if
c
      return
      end subroutine ummdp_karafillis_boyce
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     PRINCIPAL STRESSES AND INVARIANTS BY FRANÇOIS VIETE METHOD
c
      subroutine ummdp_karafillis_boyce_principal_stress ( stress,invar,
     1                                                     pStress )
c-----------------------------------------------------------------------
      implicit none
c
      real*8,intent(in) :: stress(6)
c
      real*8,intent(out) :: invar(3),pStress(3)
c
      real*8 p,q,alpha,c,dum,pi,tol
c-----------------------------------------------------------------------
c
      pi = acos(-1.0d0)
      tol = 1.0e-5
c
      invar(1) = stress(1) + stress(2) + stress(3)
      invar(2) = stress(1)*stress(2) + stress(2)*stress(3)
     1   + stress(1)*stress(3) - stress(4)**2 - stress(5)**2
     2   - stress(6)**2
      invar(3) = stress(1)*stress(2)*stress(3)
     1   + 2d0*stress(4)*stress(5)*stress(6) - stress(1)*stress(5)**2
     2   - stress(2)*stress(6)**2 - stress(3)*stress(4)**2
      p = invar(1)**2/9d0 - invar(2)/3d0
      q = invar(1)**3/27d0 + 0.5d0*invar(3) - invar(1)*invar(2)/6d0
      if ( p <= tol*abs(q) ) then
        pStress(1) = (2d0*q)**(1d0/3d0) + invar(1)/3d0
        pStress(2) = pStress(1)
        pStress(3) = pStress(1)
      else
        dum = q  /sqrt(p)**3
        if ( abs(dum) > 1.0d0 ) then
          if ( abs(abs(dum)-1.0d0) <= tol ) then
            dum = dum / abs(dum)
          else
            call ummdp_exit ( 1000 )
          end if
        end if
        alpha = acos(dum) / 3.0d0
        c = 2.0d0 * sqrt(p)
        pStress(1) = c*cos(alpha) + invar(1)/3d0
        pStress(2) = c*cos(alpha+2d0/3d0*PI) + invar(1)/3d0
        pStress(3) = c*cos(alpha+4d0/3d0*PI) + invar(1)/3d0
      end if
c
      return
      end subroutine ummdp_karafillis_boyce_principal_stress
c
c
c************************************************************************
c     VON MISES YIELD FUNCTION AND DERIVATIVES
c
      subroutine ummdp_mises ( s,se,dseds,d2seds2,nreq )
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: nreq
      real*8 ,intent(in) :: s(6)
c
      real*8,intent(out) :: se
      real*8,intent(out) :: dseds(6)
			real*8,intent(out) :: d2seds2(6,6)
c
      integer i,j
			real*8 phi
      real*8 v(6)
			real*8 c(6,6)
c-----------------------------------------------------------------------
c	
c                                               ---- coefficients matrix
      call ummdp_utility_clear2 ( c,6,6 )
      do i = 1,3
        do j = 1,3
          c(i,j) = -0.5d0
        end do
        c(i,i) = 1.0d0
      end do
      do i = 4,6
        c(i,i) = 3.0d0
      end do
c
      call ummdp_utility_mv  ( v,c,s,6,6 )
      call ummdp_utility_vvs ( phi,s,v,6 )
c                                                 ---- equivalent stress
      se = sqrt(phi)
c                                              ---- 1st order derivative
      if ( nreq >= 1 ) then
        do i = 1,6
          dseds(i) = v(i) / se
        end do
      end if
c                                              ---- 2nd order derivative
      if ( nreq >= 2 ) then
        do i = 1,6
          do j = 1,6
            d2seds2(i,j) = (-v(i)*v(j)/phi+c(i,j)) / se
          end do
        end do
      end if
c
      return
      end subroutine ummdp_mises
c
c
c
************************************************************************
c     YEGTER YIELD FUNCTION AND DERIVATIVES
c
c       doi: https://doi.org/10.1016/j.ijplas.2005.04.009
c
      subroutine ummdp_vegter ( s,se,dseds,d2seds2,nreq,pryld,ndyld )
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
      integer nf
c-----------------------------------------------------------------------
c
      nf = nint(pryld(2)) - 1
      call ummdp_vegter_core ( s,se,dseds,d2seds2,nreq,pryld,ndyld,nf )           
c
      return
      end subroutine ummdp_vegter
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     VEGTER CORE SUBROUTINE
c
      subroutine ummdp_vegter_core ( s,se,dseds,d2seds2,nreq,
     1                               pryld,ndyld,nf )
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: nreq,ndyld,nf
      real*8 ,intent(in) :: s(3),pryld(ndyld)
c
      real*8,intent(out) :: se
      real*8,intent(out) :: dseds(3)
			real*8,intent(out) :: d2seds2(3,3)
c
      integer i,j,k,m,n,isflag,iareaflag,ithetaflag
      real*8 pi,tol0,tol2,tol,vsqrt,vcos2t,vsin2t,theta,theta_rv,f_bi0,
     1       r_bi0,fun1,fun2,run,fsh1,fsh2,rsh,fps1,fps2,rps,fbi1,fbi2,
     2       rbi,fun1r,fun2r,runr,fps1r,fps2r,rpsr,alfa,beta,aa,bb,cc,
     3       dd,dmdctmp,dndctmp,nnmm,d2mdc2tmp,d2ndc2tmp
      real*8 x(4),a(2),b(2),c(2),mm(2),nn(2),mu,f(2),dphidx(3),dfdmu(2),
     1       dadc(2),dcdc(2),dbdc(2),dndc(2),dmdc(2),P(2),dfdc(2),
     2       d2adc2(2),d2bdc2(2),d2cdc2(2),d2ndc2(2),d2mdc2(2)
      real*8 dxds(3,3),d2phidx2(3,3)
      real*8 phi_un(0:nf),phi_sh(0:nf),phi_ps(0:nf),omg(0:nf),
     1       vvtmp(0:6),vvtmp_rv(0:6)
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
c-----------------------------------------------------------------------
c
      pi = acos(-1.0d0)
      tol = 1.0d-4       ! exception treatment tolerance of vsin2t
      tol0 = 1.0d-8      ! exception treatment tolerance of stress state 
      tol2 = 1.0d-2      ! se tolerance change f(1) to f(2)
c
      f_bi0=     pryld(3)
      r_bi0=     pryld(4)
      do i=0,nf
        phi_un(i)=pryld(4+i*4+1)
        phi_sh(i)=pryld(4+i*4+2)
        phi_ps(i)=pryld(4+i*4+3)
        omg(   i)=pryld(4+i*4+4)
      end do
c
      se=0.0
      do i=1,3
        se=se+s(i)**2
      end do
      if ( se<=0.0 ) then
        se=0.0
        return
      end if

c                                       ---- calc x(i) from eq.(14)~(17)
c     x(i)    : Principal stress(i=1^2)
c             : cos2theta(i=3), sin2theta(i=4)
c     isflag  : 0 s(i,i=1,3)=0
c             : 1 not s(i)=0
c                                 ---- exception treatment if all s(i)=0
c
      if(abs(s(1))<=TOL0.and.abs(s(2))<=TOL0.and.
     1                                         abs(s(3))<=TOL0)then
      isflag=0
      call ummdp_utility_clear1 ( x,4 )
      goto 100
c
      else
      isflag=1
      call ummdp_utility_clear1 ( x,4 )
c
c                           ---- exception treatment if s(1)=s(2),s(3)=0
c
      if(abs(s(1)-s(2))<=TOL0.and.abs(s(3))<=TOL0) then
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
        end if
      end if
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
      end do
        fun2=0.0d0
        fun2r=0.0d0
c                                                        ! sh=pure shear
        fsh1=0.0d0
        fsh2=0.0d0
      do m=0,nf
        fsh1=fsh1+phi_sh(m)*cos(2.0d0*dble(m)*theta)
        fsh2=fsh2-phi_sh(m)*cos(2.0d0*dble(m)*theta_rv)
      end do
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
      end do
        rps=-0.0d0
        rpsr=0.0d0
c                                                      ! bi=equi-biaxial
        fbi1=f_bi0
        fbi2=fbi1
        rbi=((r_bi0+1.0d0)+(r_bi0-1.0d0)*vcos2t)/
     &           ((r_bi0+1.0d0)-(r_bi0-1.0d0)*vcos2t)
c
c                                 ---- case distribution by stress state
      if(x(1)/=0.0d0)then
        alfa=x(2)/x(1)
      end if
      if(x(2)/=0.0d0)then
        beta=x(1)/x(2)
      end if
c
c     iareaflag    :stress state flag(i=0~6)
c
      if(x(1)>0.0d0.and.alfa<0.0d0.and.alfa>=fsh2/fsh1) then
        iareaflag=1
      else if(x(1)>0.0d0.and.alfa>=0.0d0
     1                           .and.alfa<fps2/fps1) then
        iareaflag=2
      else if(x(1)>0.0d0.and.alfa>=fps2/fps1
     1                           .and.alfa<=1.0d0) then
        iareaflag=3
c
      else if(x(1)<0.0d0.and.alfa>=1.0d0
     1                           .and.alfa<fps1r/fps2r) then
        iareaflag=4
      else if(x(1)<0.0d0.and.beta<=fps2r/fps1r
     1                                      .and.beta>0.0d0) then
        iareaflag=5
      else if(x(1)>=0.0d0.and.beta<=0.0d0
     1                           .and.beta>fsh1/fsh2) then
        iareaflag=6
c
      else
        go to 100
      end if
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
         call ummdp_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
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
         call ummdp_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
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
         call ummdp_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
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
         call ummdp_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
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
         call ummdp_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
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
         call ummdp_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
c
      case default
         write (6,*) 'iareaflag error :',iareaflag
         call ummdp_exit (9000)
      end select

c                            ---- calc. fourier coefficient mu(0<=mu<=1)
      call ummdp_vegter_calc_mu ( x,a,b,c,mu,iareaflag,s,theta
     1                                            ,aa,bb,cc,dd )
c
c                            ---- calc. normalized yield locus f(i)i=1~2
      call ummdp_vegter_calc_fi ( x,a,b,c,mu,f )
      go to 200
c
c                                                 ---- equivalent stress
  100 continue
      se=0.0d0
      go to 300
  200 continue
      if(f(1)<=TOL2) then
         se=x(2)/f(2)
      else
         se=x(1)/f(1)
      end if
c
      go to 300
c
  300 continue
c
c                                            ---- 1st order differential
c
      if ( nreq>=1 ) then
c                        ---- set dadc,dcdc,dndc,dmdc for eq.(A.7)^(A.9)
c
      call ummdp_utility_clear1 ( dadc,2 )
      call ummdp_utility_clear1 ( dbdc,2 )
      call ummdp_utility_clear1 ( dcdc,2 )
      call ummdp_utility_clear1 ( dndc,2 )
      call ummdp_utility_clear1 ( dmdc,2 )
c
      select case ( iareaflag )
c                                            
      case ( 1 )                                           ! iareaflag=1
       if(abs(vsin2t)>=TOL) then
        do m=0,nf
          dadc(1)=dadc(1)+phi_sh(m)*dble(m)
     1                             *sin(2.0d0*dble(m)*theta)
     2                             /sin(2.0d0*theta)
          dadc(2)=dadc(2)+phi_sh(m)*dble(m)
     1                             *sin(2.0d0*dble(m)*theta_rv)
     2                             /sin(2.0d0*theta)
        end do
       else
          do m=0,nf
              dadc(1)=dadc(1)+phi_sh(m)*dble(m)**2
              dadc(2)=dadc(2)+phi_sh(m)*dble(m)**2
          end do
       end if
c
       if(abs(vsin2t)>=TOL) then
         do m=0,nf
          dcdc(1)=dcdc(1)+phi_un(m)*dble(m)
     1                             *sin(2.0d0*dble(m)*theta)
     2                             /sin(2.0d0*theta)
         end do
        else
          do m=0,nf
              dcdc(1)=dcdc(1)+phi_un(m)*dble(m)**2
          end do
        end if
          dcdc(2)=0.0d0
c
          dndc(1)=0.0d0
          dndc(2)=0.0d0
c
          dmdc(1)=0.0d0
       if(abs(vsin2t)>=TOL) then
         do m=0,nf
          dmdc(2)=dmdc(2)+omg(m)*dble(m)
     1                          *sin(2.0d0*dble(m)*theta)
     2                          /sin(2.0d0*theta)
         end do
        else
          do m=0,nf
              dmdc(2)=dmdc(2)+omg(m)*dble(m)**2
          end do
        end if
      call ummdp_vegter_calc_dbdc ( a,b,c,dadc,dbdc,dcdc,mm,nn,
     1                               dndc,dmdc,iareaflag,P,nnmm)
c
c                                            
      case ( 2 )                                           ! iareaflag=2
       if(abs(vsin2t)>=TOL) then
         do m=0,nf
          dadc(1)=dadc(1)+phi_un(m)*dble(m)
     1                             *sin(2.0d0*dble(m)*theta)
     2                             /sin(2.0d0*theta)
         end do
        else
         do m=0,nf
          dadc(1)=dadc(1)+phi_un(m)*dble(m)**2
         end do
       end if
          dadc(2)=0.0d0
c
       if(abs(vsin2t)>=TOL) then
         do m=0,nf
          dcdc(1)=dcdc(1)+phi_ps(m)*dble(m)
     1                             *sin(2.0d0*dble(m)*theta)
     2                             /sin(2.0d0*theta)
         end do
        else
         do m=0,nf
          dcdc(1)=dcdc(1)+phi_ps(m)*dble(m)**2
         end do
       end if
          dcdc(2)=0.5d0*dcdc(1)
c
          dndc(1)=0.0d0
       if(abs(vsin2t)>=TOL) then
        do m=0,nf
          dndc(2)=dndc(2)+omg(m)*dble(m)
     1                           *sin(2.0d0*dble(m)*theta)
     2                           /sin(2.0d0*theta)
         end do
        else
         do m=0,nf
          dndc(2)=dndc(2)+omg(m)*dble(m)**2
         end do
       end if
c
          dmdc(1)=0.0d0
          dmdc(2)=0.0d0
      call ummdp_vegter_calc_dbdc ( a,b,c,dadc,dbdc,dcdc,mm,nn,
     1                               dndc,dmdc,iareaflag,P,nnmm)
c
c                                           
      case ( 3 )                                           ! iareaflag=3
       if(abs(vsin2t)>=TOL) then
        do m=0,nf
          dadc(1)=dadc(1)+phi_ps(m)*dble(m)
     1                             *sin(2.0d0*dble(m)*theta)
     2                             /sin(2.0d0*theta)
         end do
        else
         do m=0,nf
          dadc(1)=dadc(1)+phi_ps(m)*dble(m)**2
         end do
       end if
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
      call ummdp_vegter_calc_dbdc ( a,b,c,dadc,dbdc,dcdc,mm,nn,
     1                               dndc,dmdc,iareaflag,P,nnmm)
c
c                                            
      case ( 4 )                                           ! iareaflag=4
          dadc(1)=0.0d0
          dadc(2)=0.0d0
c
       if(abs(vsin2t)>=TOL) then
        do m=0,nf
          dcdc(2)=dcdc(2)+phi_ps(m)*dble(m)
     1                             *sin(2.0d0*dble(m)*theta_rv)
     2                             /sin(2.0d0*theta)
         end do
        else
         do m=0,nf
          dcdc(2)=dcdc(2)+phi_ps(m)*dble(m)**2
         end do
       end if
          dcdc(1)=0.5d0*dcdc(2)
c
          dndc(1)=0.0d0
          dndctmp=r_bi0+1.0d0-(r_bi0-1.0d0)*vcos2t
          dndc(2)=-2.0d0*(r_bi0*r_bi0-1.0d0)/(dndctmp*dndctmp)
c
          dmdc(1)=0.0d0
          dmdc(2)=0.0d0
      call ummdp_vegter_calc_dbdc ( a,b,c,dadc,dbdc,dcdc,mm,nn,
     1                               dndc,dmdc,iareaflag,P,nnmm)
c
c                                            
      case ( 5 )                                           ! iareaflag=5
       if(abs(vsin2t)>=TOL) then
        do m=0,nf
          dadc(2)=dadc(2)+phi_ps(m)*dble(m)
     1                             *sin(2.0d0*dble(m)*theta_rv)
     2                             /sin(2.0d0*theta)
         end do
        else
         do m=0,nf
          dadc(2)=dadc(2)+phi_ps(m)*dble(m)**2
         end do
       end if
          dadc(1)=0.5d0*dadc(2)
c
          dcdc(1)=0.0d0
       if(abs(vsin2t)>=TOL) then
        do m=0,nf
          dcdc(2)=dcdc(2)+phi_un(m)*dble(m)
     1                             *sin(2.0d0*dble(m)*theta_rv)
     2                             /sin(2.0d0*theta)
         end do
        else
         do m=0,nf
          dcdc(2)=dcdc(2)+phi_un(m)*dble(m)**2
         end do
       end if
c
          dndc(1)=0.0d0
          dndc(2)=0.0d0
c
       if(abs(vsin2t)>=TOL) then
        do m=0,nf
          dmdc(1)=dmdc(1)+omg(m)*dble(m)
     1                          *sin(2.0d0*dble(m)*theta_rv)
     2                          /sin(2.0d0*theta)
         end do
        else
         do m=0,nf
          dmdc(1)=dmdc(1)+omg(m)*dble(m)**2
         end do
       end if
          dmdc(2)=0.0d0
      call ummdp_vegter_calc_dbdc ( a,b,c,dadc,dbdc,dcdc,mm,nn,
     1                               dndc,dmdc,iareaflag,P,nnmm)
c
c                                            
      case ( 6 )                                           ! iareaflag=6
          dadc(1)=0.0d0
       if(abs(vsin2t)>=TOL) then
        do m=0,nf
          dadc(2)=dadc(2)+phi_un(m)*dble(m)
     1                             *sin(2.0d0*dble(m)*theta_rv)
     2                             /sin(2.0d0*theta)
         end do
        else
         do m=0,nf
          dadc(2)=dadc(2)+phi_un(m)*dble(m)**2
         end do
       end if
c
       if(abs(vsin2t)>=TOL) then
        do m=0,nf
          dcdc(1)=dcdc(1)+phi_sh(m)*dble(m)
     1                             *sin(2.0d0*dble(m)*theta)
     2                             /sin(2.0d0*theta)
          dcdc(2)=dcdc(2)+phi_sh(m)*dble(m)
     1                             *sin(2.0d0*dble(m)*theta_rv)
     2                             /sin(2.0d0*theta)
         end do
        else
         do m=0,nf
          dcdc(1)=dcdc(1)+phi_sh(m)*dble(m)**2
          dcdc(2)=dcdc(2)+phi_sh(m)*dble(m)**2
         end do
       end if
c
       if(abs(vsin2t)>=TOL) then
        do m=0,nf
          dndc(1)=dndc(1)+omg(m)*dble(m)
     1                          *sin(2.0d0*dble(m)*theta_rv)
     2                          /sin(2.0d0*theta)
         end do
        else
         do m=0,nf
          dndc(1)=dndc(1)+omg(m)*dble(m)**2
         end do
       end if
          dndc(2)=0.0d0
c
          dmdc(1)=0.0d0
          dmdc(2)=0.0d0
      call ummdp_vegter_calc_dbdc ( a,b,c,dadc,dbdc,dcdc,mm,nn,
     1                               dndc,dmdc,iareaflag,P,nnmm)
c
      case default
        write (6,*) 'iareaflag error(dseds) :',iareaflag
        call ummdp_exit (9000)
      end select
c
c
c                             ---- calc. dphidx(i) (i=1~3)  eq.(21)~(23)
      call ummdp_vegter_calc_dphidx ( dphidx,se,a,b,c,dadc,dbdc,dcdc,
     1                      mm,nn,dndc,dmdc,f,iareaflag,mu,x,dfdmu,dfdc)
c
c
c                                   ---- calc. dseds(i) (i=1~3)  eq.(20)
      call ummdp_vegter_calc_dseds ( dseds,x,dphidx,vcos2t,vsin2t,
     1                                iareaflag,isflag,dxds)
c
      end if
c
c
c                                            ---- 2nd order differential
      if ( nreq>=2 ) then
c                     ---- set d2adc2,d2cdc2,d2ndc2,d2mdc2 for d2bdc2(2)
c
      call ummdp_utility_clear1 ( d2adc2,2 )
      call ummdp_utility_clear1 ( d2bdc2,2 )
      call ummdp_utility_clear1 ( d2cdc2,2 )
      call ummdp_utility_clear1 ( d2ndc2,2 )
      call ummdp_utility_clear1 ( d2mdc2,2 )
c
c                     ---- define exception treatment condition of theta
c                        ---- if theta<=0.002865deg then apply exception
c
      if(abs(vsin2t)<=TOL) then
         ithetaflag=1
      else 
         ithetaflag=0
      end if
c
      if(ithetaflag==1) then
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
     1          (sin(2.0d0*theta)**3)-dble(m)*cos(2.0d0*dble(m)*theta)/
     2                                           (sin(2.0d0*theta)**2)
        vvtmp_rv(m)=cos(2.0d0*theta)*sin(2.0d0*dble(m)*theta_rv)/
     1       (sin(2.0d0*theta)**3)+dble(m)*cos(2.0d0*dble(m)*theta_rv)/
     2                                           (sin(2.0d0*theta)**2)
       end do
c
      end if
c
      select case ( iareaflag )
c                                            
      case ( 1 )                                           ! iareaflag=1
         do m=0,nf
          d2adc2(1)=d2adc2(1)+phi_sh(m)*dble(m)*vvtmp(m)
          d2adc2(2)=d2adc2(2)+phi_sh(m)*dble(m)*vvtmp_rv(m)
         end do
c
         do m=0,nf
          d2cdc2(1)=d2cdc2(1)+phi_un(m)*dble(m)*vvtmp(m)
         end do
          d2cdc2(2)=0.0d0
c
          d2ndc2(1)=0.0d0
          d2ndc2(2)=0.0d0
c
          d2mdc2(1)=0.0d0
         do m=0,nf
          d2mdc2(2)=d2mdc2(2)+omg(m)*dble(m)*vvtmp(m)
         end do
      call ummdp_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,nn,dndc,
     1     dmdc,iareaflag,d2adc2,d2bdc2,d2cdc2,d2ndc2,d2mdc2,P,nnmm,s)
c
c                                            
      case ( 2 )                                           ! iareaflag=2
         do m=0,nf
          d2adc2(1)=d2adc2(1)+phi_un(m)*dble(m)*vvtmp(m)
         end do
          d2adc2(2)=0.0d0
c
         do m=0,nf
          d2cdc2(1)=d2cdc2(1)+phi_ps(m)*dble(m)*vvtmp(m)
         end do
          d2cdc2(2)=0.5d0*d2cdc2(1)
c
          d2ndc2(1)=0.0d0
         do m=0,nf
          d2ndc2(2)=d2ndc2(2)+omg(m)*dble(m)*vvtmp(m)
         end do
c
          d2mdc2(1)=0.0d0
          d2mdc2(2)=0.0d0
      call ummdp_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,nn,dndc,
     1     dmdc,iareaflag,d2adc2,d2bdc2,d2cdc2,d2ndc2,d2mdc2,P,nnmm,s)
c
c                                            
      case ( 3 )                                           ! iareaflag=3
         do m=0,nf  
          d2adc2(1)=d2adc2(1)+phi_ps(m)*dble(m)*vvtmp(m)
         end do
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
     1                                              (d2mdc2tmp**3)
      call ummdp_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,nn,dndc,
     1     dmdc,iareaflag,d2adc2,d2bdc2,d2cdc2,d2ndc2,d2mdc2,P,nnmm,s)
c
c                                            
      case ( 4 )                                          !  iareaflag=4
          d2adc2(1)=0.0d0
          d2adc2(2)=0.0d0
c
         do m=0,nf
          d2cdc2(2)=d2cdc2(2)+phi_ps(m)*dble(m)*vvtmp_rv(m)
         end do
          d2cdc2(1)=0.5d0*d2cdc2(2)
c
          d2ndc2(1)=0.0d0
             d2ndc2tmp=r_bi0+1.0d0-(r_bi0-1.0d0)*vcos2t
          d2ndc2(2)=-4.0d0*(r_bi0**2-1.0d0)*(r_bi0-1.0d0)/
     1                                              (d2ndc2tmp**3)
c
          d2mdc2(1)=0.0d0
          d2mdc2(2)=0.0d0
      call ummdp_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,nn,dndc,
     1     dmdc,iareaflag,d2adc2,d2bdc2,d2cdc2,d2ndc2,d2mdc2,P,nnmm,s)
c
c                                           
      case ( 5 )                                           ! iareaflag=5
         do m=0,nf
          d2adc2(2)=d2adc2(2)+phi_ps(m)*dble(m)*vvtmp_rv(m)
         end do
          d2adc2(1)=0.5d0*d2adc2(2)
c
          d2cdc2(1)=0.0d0
         do m=0,nf
          d2cdc2(2)=d2cdc2(2)+phi_un(m)*dble(m)*vvtmp_rv(m)
         end do
c
          d2ndc2(1)=0.0d0
          d2ndc2(2)=0.0d0
c
         do m=0,nf
          d2mdc2(1)=d2mdc2(1)+omg(m)*dble(m)*vvtmp_rv(m)
         end do
          d2mdc2(2)=0.0d0
      call ummdp_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,nn,dndc,
     1     dmdc,iareaflag,d2adc2,d2bdc2,d2cdc2,d2ndc2,d2mdc2,P,nnmm,s)
c
c                                            
      case ( 6 )                                           ! iareaflag=6
          d2adc2(1)=0.0d0
         do m=0,nf
          d2adc2(2)=d2adc2(2)+phi_un(m)*dble(m)*vvtmp_rv(m)
         end do
c
         do m=0,nf
          d2cdc2(1)=d2cdc2(1)+phi_sh(m)*dble(m)*vvtmp(m)
          d2cdc2(2)=d2cdc2(2)+phi_sh(m)*dble(m)*vvtmp_rv(m)
         end do
cn
         do m=0,nf
          d2ndc2(1)=d2ndc2(1)+omg(m)*dble(m)*vvtmp_rv(m)
         end do
          d2ndc2(2)=0.0d0
c
          d2mdc2(1)=0.0d0
          d2mdc2(2)=0.0d0
      call ummdp_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,nn,dndc,
     1     dmdc,iareaflag,d2adc2,d2bdc2,d2cdc2,d2ndc2,d2mdc2,P,nnmm,s)
c
      case default
        write (6,*) 'iareaflag error(d2seds2) :',iareaflag
        call ummdp_exit (9000)
      end select
c
c
c                                     ---- calc. d2phidx2(k,l) (k,l=1~3)
      call ummdp_vegter_calc_d2phidx2 (d2phidx2,se,a,b,c,dadc,dbdc,
     1             dcdc,mm,nn,dndc,dmdc,f,iareaflag,mu,x,d2adc2,d2bdc2,
     2            d2cdc2,d2ndc2,d2mdc2,dfdmu,dfdc,s,aa,bb,cc,dd,dphidx)
c
c
c                                      ---- calc. d2seds2(i,j) (i,j=1~3)
      call ummdp_vegter_calc_d2seds2 (d2seds2,d2phidx2,se,a,b,c,mu,x,
     1         vcos2t,vsin2t,iareaflag,dxds,dphidx,isflag,s,dseds,
     2                                 pryld,ndyld)
c
      end if
c
      return
      end subroutine ummdp_vegter_core
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CALCULATE HINGEPOINT b(i,i=1~2)
c
      subroutine ummdp_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: iareaflag
      real*8 ,intent(in) :: a(2),c(2),mm(2),nn(2),s(3)
c
      real*8,intent(out) :: b(2)
c
      real*8 tol,bb,b1u,b2u
c-----------------------------------------------------------------------
c
      tol = 1.0d-8
c
      b1u = mm(2)*(nn(1)*a(1)+nn(2)*a(2))-nn(2)*(mm(1)*c(1)+mm(2)*c(2))
      b2u = nn(1)*(mm(1)*c(1)+mm(2)*c(2))-mm(1)*(nn(1)*a(1)+nn(2)*a(2))
c
      bb = nn(1)*mm(2)-mm(1)*nn(2)
      if ( abs(bb) <= TOL ) then
         write (6,*) 'hingepoint singular error! '
         call ummdp_exit (9000)
      end if
c
      b(1) = b1u / bb
      b(2) = b2u / bb
c
      return
      end subroutine ummdp_vegter_hingepoint
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CALCULATE FOURIER COEFFICIENT mu(0<=mu<=1)
c
      subroutine ummdp_vegter_calc_mu ( x,a,b,c,mu,iareaflag,s,theta,
     1                                  aa,bb,cc,dd)
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: iareaflag
      real*8 ,intent(in) :: theta 
      real*8 ,intent(in) :: x(4),a(2),b(2),c(2),s(3)
c
      real*8,intent(out) :: mu
c
      integer imuflag
      real*8 tol1,tol2,aa,bb,cc,dd
      real*8 xx(2)
c-----------------------------------------------------------------------
c
      tol1 = 1.0d-8
      tol2 = 1.0d-8
c
      aa = x(2)*(a(1)+c(1)-2.0d0*b(1))-x(1)*(a(2)+c(2)-2.0d0*b(2))
      bb = 2.0d0*x(2)*(b(1)-a(1))-2.0d0*x(1)*(b(2)-a(2))
      cc = x(2)*a(1)-x(1)*a(2)
c
      if ( abs(aa) <= tol1 ) then
         write (6,*) 'calc. mu singular error! ',abs(aa),iareaflag
         call ummdp_exit (9000)
      end if
c
      dd = bb*bb - 4.0d0*aa*cc
      if ( dd >= 0.0d0 ) then
        xx(1) = 0.5d0 * (-bb+sign(sqrt(dd),-bb))/aa
        xx(2) = cc / (aa*xx(1))
c
      else
         write (6,*) 'negative dd ! ',dd,iareaflag
         call ummdp_exit (9000)
      end if
c
      if ( xx(1) >= 0.0d0 .and. xx(1) <= 1.0000005d0 ) then
        mu = xx(1)
        imuflag = 1
      else if ( xx(2) >= 0.0d0 .and. xx(2) <= 1.0000005d0 ) then
        mu = xx(2)
        imuflag = 2
      else if ( abs(xx(1) ) <= tol2 .or. abs(xx(2)) <= tol2 ) then
        mu = 0.0d0
      else
        write (6,*) 'can not find mu ! solve error ',iareaflag,xx(1)
     1              ,xx(2)
         call ummdp_exit (9000)
      end if
c
      return
      end subroutine ummdp_vegter_calc_mu
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CALCULATE NORMALIZED YIELD LOCUS
c
      subroutine ummdp_vegter_calc_fi ( x,a,b,c,mu,f )
c-----------------------------------------------------------------------
      implicit none
c
      real*8,intent(in) :: mu
      real*8,intent(in) :: x(4),a(2),b(2),c(2)
c
      real*8,intent(out) :: f(2)
c
      integer i
c-----------------------------------------------------------------------
c
      do i = 1,2
        f(i) = a(i)+ 2.0d0*mu*(b(i)-a(i)) + mu*mu*(a(i)+c(i)-2.0d0*b(i))
      end do
c
      return
      end subroutine ummdp_vegter_calc_fi
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CAKCULATE dbdc(i) (i=1~2) eq.(A.7)
c
      subroutine ummdp_vegter_calc_dbdc ( a,b,c,dadc,dbdc,dcdc,mm,nn,
     1                                    dndc,dmdc,iareaflag,P,nnmm )
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: iareaflag
      real*8 ,intent(in) :: a(2),b(2),c(2),dadc(2),dcdc(2),mm(2),
     1                      nn(2),dndc(2),dmdc(2)
c
      real*8,intent(out) :: nnmm
      real*8,intent(out) :: dbdc(2),P(2)     
c
      integer i
      real*8 tol,nminv
c-----------------------------------------------------------------------
c
      tol = 1.0d-8
c
      P(1) = nn(1)*dadc(1) + dndc(1)*(a(1)-b(1)) + nn(2)*dadc(2)
     1       + dndc(2)*(a(2)-b(2))
c
      P(2) = mm(1)*dcdc(1) + dmdc(1)*(c(1)-b(1)) + mm(2)*dcdc(2)
     1       + dmdc(2)*(c(2)-b(2))
c
      nnmm = nn(1)*mm(2) - mm(1)*nn(2)
         if ( abs(nnmm) < tol ) then
            write (6,*) 'nnmm too small! ',nnmm
            call ummdp_exit (9000)
         end if
      nminv = 1.0d0 / nnmm
c
      dbdc(1) = nminv * (P(1)*mm(2)-P(2)*nn(2))
      dbdc(2) = nminv * (P(2)*nn(1)-P(1)*mm(1))
c
      return
      end subroutine ummdp_vegter_calc_dbdc
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     calc. dphidx(i) (i=1~3)  eq.(21)~(23)
c
      subroutine ummdp_vegter_calc_dphidx ( dphidx,se,a,b,c,dadc,dbdc,
     1                                      dcdc,mm,nn,dndc,dmdc,f,  
     2                                      iareaflag,mu,x,dfdmu,dfdc )
c-----------------------------------------------------------------------   
      implicit none
c
      integer,intent(in) :: iareaflag
      real*8 ,intent(in) :: se
      real*8 ,intent(in) :: a(2),b(2),c(2),dadc(2),dbdc(2),dcdc(2),
     1                      mm(2),nn(2),dndc(2),dmdc(2),f(2),mu,x(4)
c
      real*8,intent(out) :: dphidx(3),dfdmu(2),dfdc(2)
c
      integer i
      real*8 tol1,tol2,tol3,dphidxcoe,dphidxcinv,tmp_u,tmp_b,vtan
      real*8 dphidxtmp(3)
c-----------------------------------------------------------------------       
c
      tol1 = 0.996515d0         ! =44.9deg
      tol2 = 1.003497d0         ! =45.1deg
      tol3 = 1.0d-8
c                                    ---- calc. dfdc(i) (i=1~2)  eq.(23)
      do i = 1,2
        dfdc(i) = dadc(i) + 2.0d0*mu*(dbdc(i)-dadc(i)) + mu*mu*(dadc(i)
     1            + dcdc(i)-2.0d0*dbdc(i))
      end do
c
c                                   ---- calc. dfdmu(i) (i=1~2)  eq.(22)
      do i = 1,2
        dfdmu(i) = 2.0d0*(b(i)-a(i)) + 2.0d0*mu*(a(i)+c(i)-2.0d0*b(i))
      end do
c
c                            ---- calc. dphidx(i) (i=1~3)  eq.(21),(C.1)
      dphidxcinv = f(1)*dfdmu(2) - f(2)*dfdmu(1)
      if ( abs(dphidxcinv) < tol3 ) then
        write (6,*) 'eq.(21) too small! ',dphidxcinv
        call ummdp_exit (9000)
      end if
      dphidxcoe = 1.0d0 / dphidxcinv
c
c                            ---- if condition to avoid singular eq.(20)
c                                            ---- apply 44.9 to 45.1 deg
c
      if ( iareaflag == 3 .or. iareaflag == 4) then
        vtan = x(2) / x(1)
      end if
c
      if ( iareaflag == 4 .and. vtan >= tol1 .and. vtan <= tol2 ) then
c
        tmp_u = 1.0d0*(2.0d0*(1.0d0-mu)*dbdc(2)+mu*dcdc(2))*dfdmu(1)
     1          - 1.0d0*(2.0d0*(1.0d0-mu)*dbdc(1)+mu*dcdc(1))*dfdmu(2)
c
        tmp_b = 2.0d0*(1.0d0-mu)*(b(1)-b(2)) + mu*(c(1)-c(2))

        dphidxtmp(1) =  dfdmu(2)
        dphidxtmp(2) = -dfdmu(1)
        dphidxtmp(3) = tmp_u / tmp_b
c
      else if ( iareaflag == 3 .and. vtan >= tol1 
     1          .and. vtan <= tol2 ) then
c
        tmp_u = 1.0d0*(2.0d0*mu*dbdc(2)+(1.0d0-mu)*dadc(2))*dfdmu(1)
     1          - 1.0d0*(2.0d0*mu*dbdc(1)+(1.0d0-mu)*dadc(1))*dfdmu(2)
c
        tmp_b = 2.0d0*mu*(b(1)-b(2))+(1.0d0-mu)*(a(1)-a(2))

        dphidxtmp(1) =  dfdmu(2)
        dphidxtmp(2) = -dfdmu(1)
        dphidxtmp(3) = tmp_u / tmp_b
c
      else
c
        dphidxtmp(1) =  dfdmu(2)
        dphidxtmp(2) = -dfdmu(1)
        dphidxtmp(3) = se * (dfdc(2)*dfdmu(1)-dfdc(1)*dfdmu(2))
c
      end if
      do i = 1,3
        dphidx(i) = dphidxcoe * dphidxtmp(i)
      end do
c
      return
      end subroutine ummdp_vegter_calc_dphidx
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CALCULATE FIRST ORDER DERIVATIVE
c
      subroutine ummdp_vegter_calc_dseds ( dseds,x,dphidx,vcos2t,
     1                                     vsin2t,iareaflag,isflag,
     2                                     dxds )
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: iareaflag,isflag
      real*8 ,intent(in) :: vcos2t,vsin2t
      real*8 ,intent(in) :: x(4),dphidx(3)
c
      real*8,intent(out) :: dseds(3)
      real*8,intent(out) :: dxds(3,3) 
c
      integer i,j
      real*8 tol1,tol2,vtan
      real*8 dxds_t(3,3)  
c-----------------------------------------------------------------------
c
      tol1 = 0.996515d0       ! =44.9deg
      tol2 = 1.003497d0       ! =45.1deg
c                  ---- set linear transformation matrix dxds(3,3) eq.18
c
c                            ---- if condition to avoid singular eq.(20)
c                                            ---- apply 44.9 to 45.1 deg
c
      if ( iareaflag == 3 .or. iareaflag == 4 ) then
        vtan = x(2) / x(1)
      end if
c
      if ( iareaflag == 3 .and. vtan >= tol1 .and. vtan <= tol2 ) then
        dxds(1,1) = 0.5d0 * (1.0d0+vcos2t)
        dxds(2,1) = 0.5d0 * (1.0d0-vcos2t)
        dxds(3,1) = vsin2t * vsin2t
        dxds(1,2) = 0.5d0 * (1.0d0-vcos2t)
        dxds(2,2) = 0.5d0 * (1.0d0+vcos2t)
        dxds(3,2) = -vsin2t * vsin2t
        dxds(1,3) =  vsin2t
        dxds(2,3) = -vsin2t
        dxds(3,3) = -2.0d0 * vsin2t * vcos2t
        dxds_t = transpose(dxds)
c
      else if ( iareaflag == 4 .and. vtan >= tol1 
     1          .and. vtan <= tol2) then
        dxds(1,1) = 0.5d0 * (1.0d0+vcos2t)
        dxds(2,1) = 0.5d0 * (1.0d0-vcos2t)
        dxds(3,1) = vsin2t * vsin2t
        dxds(1,2) = 0.5d0 * (1.0d0-vcos2t)
        dxds(2,2) = 0.5d0 * (1.0d0+vcos2t)
        dxds(3,2) = -vsin2t * vsin2t
        dxds(1,3) =  vsin2t
        dxds(2,3) = -vsin2t
        dxds(3,3) = -2.0d0 * vsin2t * vcos2t
        dxds_t = transpose(dxds)
c
      else
        dxds(1,1) = 0.5d0 * (1.0d0+vcos2t)
        dxds(2,1) = 0.5d0 * (1.0d0-vcos2t)
        dxds(3,1) = vsin2t * vsin2t / (x(1)-x(2))
        dxds(1,2) = 0.5d0 * (1.0d0-vcos2t)
        dxds(2,2) = 0.5d0 * (1.0d0+vcos2t)
        dxds(3,2) = -vsin2t * vsin2t / (x(1)-x(2))
        dxds(1,3) =  vsin2t
        dxds(2,3) = -vsin2t
        dxds(3,3) = -2.0d0 * vsin2t * vcos2t / (x(1)-x(2))
        dxds_t = transpose(dxds)
      end if
c
c                                   ---- calc. dseds(i) (1=1~3)  eq.(20)
      call ummdp_utility_clear1 ( dseds,3 )
c
      do i = 1,3
        do j = 1,3
          dseds(i) = dseds(i) + dxds_t(i,j)*dphidx(j)
        end do
      end do
c
      return
      end subroutine ummdp_vegter_calc_dseds
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CALCULATE d2bdc2(i) (i=1~2)
c
      subroutine ummdp_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,nn,
     1                                      dndc,dmdc,iareaflag,d2adc2,
     2                                      d2bdc2,d2cdc2,d2ndc2,
     3                                      d2mdc2,P,nnmm,s )      
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: iareaflag
      real*8 ,intent(in) :: nnmm
      real*8 ,intent(in) :: a(2),b(2),c(2),dadc(2),dcdc(2),mm(2),
     1                      nn(2),dndc(2),dmdc(2),dbdc(2),
     2                      d2adc2(2),d2cdc2(2),d2ndc2(2),d2mdc2(2),
     3                      P(2),s(3)
c
      real*8,intent(out) :: d2bdc2(2)
c
      integer i
      real*8 dnnmmdc,dp1m2p2n2,dp2n1p1m1
      real*8 dPdc(2)
c-----------------------------------------------------------------------
c
      dnnmmdc = dndc(1)*mm(2)+nn(1)*dmdc(2)-dmdc(1)*nn(2)-mm(1)*dndc(2)
c
      dPdc(1) = dndc(1)*dadc(1) + nn(1)*d2adc2(1) 
     1          + d2ndc2(1)*(a(1)-b(1)) + dndc(1)*(dadc(1)-dbdc(1))
     2          + dndc(2)*dadc(2) + nn(2)*d2adc2(2)                      
     2          + d2ndc2(2)*(a(2)-b(2)) + dndc(2)*(dadc(2)-dbdc(2))              
c
      dPdc(2) = dmdc(1)*dcdc(1) + mm(1)*d2cdc2(1) 
     1          + d2mdc2(1)*(c(1)-b(1)) + dmdc(1)*(dcdc(1)-dbdc(1))
     2          + dmdc(2)*dcdc(2) + mm(2)*d2cdc2(2)
     3          + d2mdc2(2)*(c(2)-b(2)) + dmdc(2)*(dcdc(2)-dbdc(2))
c
      dp1m2p2n2 = dPdc(1)*mm(2) + P(1)*dmdc(2) - dPdc(2)*nn(2)
     1            - P(2)*dndc(2)
      dp2n1p1m1 = dPdc(2)*nn(1) + P(2)*dndc(1) - dPdc(1)*mm(1)
     1            - P(1)*dmdc(1)
c
      d2bdc2(1) = -1.0d0*dnnmmdc*(P(1)*mm(2)-P(2)*nn(2))/(nnmm*nnmm)
     1            + dp1m2p2n2/nnmm
      d2bdc2(2) = -1.0d0*dnnmmdc*(P(2)*nn(1)-P(1)*mm(1))/(nnmm*nnmm)
     2            + dp2n1p1m1/nnmm
c
      return
      end subroutine ummdp_vegter_calc_d2bdc2
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CALCULATE d2phidx2(k,l) (k,l=1~3)
c
      subroutine ummdp_vegter_calc_d2phidx2 ( d2phidx2,se,a,b,c,dadc,
     1                                        dbdc,dcdc,mm,nn,dndc,
     2                                        dmdc,f,iareaflag,mu,x,
     3                                        d2adc2,d2bdc2,d2cdc2,
     4                                        d2ndc2,d2mdc2,dfdmu,dfdc,
     5                                        s,aa,bb,cc,dd,dphidx )
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: iareaflag
      real*8 ,intent(in) :: se,mu,aa,bb,cc,dd
      real*8 ,intent(in) :: a(2),b(2),c(2),dadc(2),dbdc(2),dcdc(2),
     1                      mm(2),nn(2),dndc(2),dmdc(2),f(2),x(4),
     2                      d2adc2(2),d2cdc2(2),d2ndc2(2),d2mdc2(2),
     3                      d2bdc2(2),dfdmu(2),dfdc(2),s(3),dphidx(3)
c
      real*8,intent(out) :: d2phidx2(3,3)
c
      integer i,j
      real*8 vcommon,vtmp1,vtmp2,vtmp3,vtmp4,va,vc
      real*8 daadx(3),dbbdx(3),dccdx(3),ddddx(3),dmudx(3),d2fdmu2(2),
     1       d2fdmudc(2),d2fdcdmu(2),d2fdc2(2),dvadx(3),dsedx(3),
     2       dvcdx(3)
c-----------------------------------------------------------------------
c
c                                           ---- calc.  dmudx(i) (i=1~3)
      daadx(1) = -a(2) - c(2) + 2.0d0*b(2)
      dbbdx(1) = -2.0d0 * (b(2)-a(2))
      dccdx(1) = -a(2)
c
      daadx(2) = a(1) + c(1) - 2.0d0*b(1)
      dbbdx(2) = 2.0d0 * (b(1)-a(1))
      dccdx(2) = a(1)
c
      daadx(3) = x(2)*(dadc(1)+dcdc(1)-2.0d0*dbdc(1))
     1           - x(1)*(dadc(2)+dcdc(2)-2.0d0*dbdc(2))
      dbbdx(3) = 2.0d0*x(2)*(dbdc(1)-dadc(1)) 
     1           - 2.0d0*x(1)*(dbdc(2)-dadc(2))
      dccdx(3) = x(2)*dadc(1) - x(1)*dadc(2)
c
      do i = 1,3
        dmudx(i) = 0.5d0*daadx(i)*(bb+sqrt(dd))/(aa*aa)
     1          + 0.5d0*(-dbbdx(i)-0.5d0/(sqrt(dd))*(2.0d0*bb*dbbdx(i)
     2             - 4.0d0*daadx(i)*cc-4.0d0*aa*dccdx(i)))/aa
      end do
c
c                                         ---- calc.  d2fdmu2(i) (i=1~2)
      do i = 1,2
        d2fdmu2(i) = 2.0d0 * (a(i)+c(i)-2.0d0*b(i))
      end do
c
c                                        ---- calc.  d2fdmudc(i) (i=1~2)
      do i = 1,2
        d2fdmudc(i) = 2.0d0*(dbdc(i)-dadc(i))
     1                 + 2.0d0*mu*(dadc(i)+dcdc(i)- 2.0d0*dbdc(i))
      end do
c
c                                        ---- calc.  d2fdcdmu(i) (i=1~2)
      do i = 1,2
        d2fdcdmu(i)=2.0d0*(dbdc(i)-dadc(i))+2.0d0*mu*(dadc(i)+dcdc(i)
     1                                             -2.0d0*dbdc(i))
      end do
c
c                                          ---- calc.  d2fdc2(i) (i=1~2)
      do i = 1,2
        d2fdc2(i)=d2adc2(i)+2.0d0*mu*(d2bdc2(i)-d2adc2(i))
     &            +mu*mu*(d2adc2(i)+d2cdc2(i)-2.0d0*d2bdc2(i))
      end do
c
c                                                 ---- for d2phidx2(k,l)
c
      vcommon = 1.0d0/(f(1)*dfdmu(2)-f(2)*dfdmu(1))
      vtmp1 = dfdc(1)*dfdmu(2)+f(1)*d2fdmudc(2)
     &                           -dfdc(2)*dfdmu(1)-f(2)*d2fdmudc(1)
      vtmp2 = dfdmu(1)*dfdmu(2)+f(1)*d2fdmu2(2)
     &                           -dfdmu(2)*dfdmu(1)-f(2)*d2fdmu2(1)
      vtmp3 = d2fdcdmu(2)*dfdmu(1)+dfdc(2)*d2fdmu2(1)
     &                     -d2fdcdmu(1)*dfdmu(2)-dfdc(1)*d2fdmu2(2)
      vtmp4 = d2fdc2(2)*dfdmu(1)+dfdc(2)*d2fdmudc(1)
     &                      -d2fdc2(1)*dfdmu(2)-dfdc(1)*d2fdmudc(2)
c
      va=vcommon
      vc=dfdc(2)*dfdmu(1)-dfdc(1)*dfdmu(2)
c
      do i=1,2
        dvadx(i)=-vtmp2*vcommon*vcommon*dmudx(i)
      end do
        dvadx(3)=-vtmp1*vcommon*vcommon
     &           -vtmp2*vcommon*vcommon*dmudx(3)
c
      do i=1,3
        dsedx(i)=dphidx(i)
      end do
c
      do i=1,2
        dvcdx(i)=vtmp3*dmudx(i)
      end do
        dvcdx(3)=vtmp4+vtmp3*dmudx(3)
c
c                                    ---- calc.  d2phidx2(i,j) (i,j=1~3)
      do j=1,2
        d2phidx2(1,j)=(-dfdmu(2)*vtmp2*vcommon*vcommon
     &                             +d2fdmu2(2)*vcommon)*dmudx(j)
      end do
        d2phidx2(1,3)=-dfdmu(2)*vtmp1*vcommon*vcommon
     &                +d2fdmudc(2)*vcommon+(-dfdmu(2)*vtmp2
     &                  *vcommon*vcommon+d2fdmu2(2)*vcommon)*dmudx(3)
c
      do j=1,2
        d2phidx2(2,j)=(dfdmu(1)*vtmp2*vcommon*vcommon
     &                             -d2fdmu2(1)*vcommon)*dmudx(j)
      end do
        d2phidx2(2,3)=dfdmu(1)*vtmp1*vcommon*vcommon
     &                -d2fdmudc(1)*vcommon+(dfdmu(1)*vtmp2
     &                  *vcommon*vcommon-d2fdmu2(1)*vcommon)*dmudx(3)
c
      do j=1,3
        d2phidx2(3,j)=dvadx(j)*se*vc+va*dsedx(j)*vc+va*se*dvcdx(j)
      end do
c
      return
      end subroutine ummdp_vegter_calc_d2phidx2
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CALCULATE d2seds2(i,j) (i,j=1~3)
c
      subroutine ummdp_vegter_calc_d2seds2 ( d2seds2,d2phidx2,se,a,b,c,
     1                                       mu,x,vcos2t,vsin2t,
     2                                       iareaflag,dxds,dphidx,
     3                                       isflag,s,dseds,pryld,
     4                                       ndyld )      
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: iareaflag,isflag,ndyld
      real*8 ,intent(in) :: se,mu,vcos2t,vsin2t
      real*8 ,intent(in) :: a(2),b(2),c(2),x(4),dphidx(3),s(3),
     1                       dseds(3),pryld(ndyld)
      real*8 ,intent(in) :: d2phidx2(3,3),dxds(3,3)
c            
      real*8,intent(out) :: d2seds2(3,3)
c
      integer i,j,k,l,iflag 
      real*8 tol1,tol2,tol3a,tol3b,tol4,tol4a,tol4b,vtan,vx1x2
      real*8 d2xds2(3,3,3)
c-----------------------------------------------------------------------
c
      tol1 = 0.996515d0          ! =44.9deg
      tol2 = 1.003497d0          ! =45.1deg
      tol3a = 1.0d-6
      tol3b = 1.0d0 - tol3a
      tol4 = 1.0d-7
      tol4a = 1.0d0 - tol4       !thera<0.012812deg
      tol4b = - 1.0d0 + tol4     !thera>89.98719deg
c                             ---- if condition to apply numerical diff.
c
      if(iareaflag==3.or.iareaflag==4) then
      vtan=x(2)/x(1)
      end if
c
      if(iareaflag==4.and.vtan>=TOL1.and.vtan<=TOL2) then
      iflag=1
c
      else if(iareaflag==3.and.vtan>=TOL1.and.vtan<=TOL2) then
      iflag=2
c
      else if(abs(mu)<=TOL3a.or.abs(mu)>=TOL3b) then
      iflag=3
c
      else if(vcos2t>=TOL4a.or.vcos2t<=TOL4b) then
      iflag=4
c
      else
      vx1x2=x(1)-x(2)
      iflag=0
      end if
c
c                                    ---- set d2xds2(k,i,j)  (k,i,j=1~3)
c
      if (iflag==0) then
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
     1                                           (vx1x2*vx1x2)
c
      d2xds2(3,2,1)=3.0d0*vcos2t*vsin2t*vsin2t/(vx1x2*vx1x2)
      d2xds2(3,2,2)=-3.0d0*vcos2t*vsin2t*vsin2t/(vx1x2*vx1x2)
      d2xds2(3,2,3)=-2.0d0*vsin2t*(2.0d0-3.0d0*vsin2t*vsin2t)/
     1                                           (vx1x2*vx1x2)
c
      d2xds2(3,3,1)=2.0d0*vsin2t*(3.0d0*vcos2t*vcos2t-1.0d0)/
     1                                           (vx1x2*vx1x2)
      d2xds2(3,3,2)=-2.0d0*vsin2t*(3.0d0*vcos2t*vcos2t-1.0d0)/
     1                                           (vx1x2*vx1x2)
      d2xds2(3,3,3)=4.0d0*vcos2t*(3.0d0*vsin2t*vsin2t-1.0d0)/
     1                                           (vx1x2*vx1x2)
      end if
c
c                                      ---- calc. d2seds2(i,j) (i,j=1~3)
c
      call ummdp_utility_clear2 ( d2seds2,3,3 )
c
      if (iflag/=0) then
      call ummdp_vegter_d2seds2n(d2seds2,s,dseds,pryld,ndyld,se)
c
      else
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              d2seds2(i,j)=d2seds2(i,j)
     1                    +d2phidx2(k,l)*dxds(l,j)*dxds(k,i)
            end do
              d2seds2(i,j)=d2seds2(i,j)
     1                    +dphidx(k)*d2xds2(k,i,j)
          end do
        end do
      end do
      end if
c
      return
      end subroutine ummdp_vegter_calc_d2seds2
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     numerical differential for 2nd order differentials
c
      subroutine ummdp_vegter_d2seds2n ( d2seds2,s,dseds,
     1                                   pryld,ndyld,se )
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: ndyld
      real*8 ,intent(in) :: se
      real*8 ,intent(in) :: dseds(3),pryld(ndyld),s(3)
c
      real*8,intent(out) :: d2seds2(3,3)
c
      integer j,k
      real*8 delta,sea,seb,a,b,seba,seaa,sebb,seab,se0
      real*8 s0(3),ss(3)     
c-----------------------------------------------------------------------
c
      delta = 1.0d-3
c
      s0(:) = s(:)
      ss(:) = s(:)
      do j=1,3
        do k=1,3
          if ( j==k ) then
            se0=se
            ss(j)=s0(j)-delta
            call ummdp_vegter_yieldfunc(3,ss,sea,dseds,d2seds2,0,
     1                                   pryld,ndyld)
            ss(j)=s0(j)+delta
            call ummdp_vegter_yieldfunc(3,ss,seb,dseds,d2seds2,0,
     1                                    pryld,ndyld)
            ss(j)=s0(j)
            a=(se0-sea)/delta
            b=(seb-se0)/delta
            d2seds2(j,k)=(b-a)/delta
          else
            ss(j)=s0(j)-delta
            ss(k)=s0(k)-delta
            call ummdp_vegter_yieldfunc(3,ss,seaa,dseds,d2seds2,0,
     1                                    pryld,ndyld)
            ss(j)=s0(j)+delta
            ss(k)=s0(k)-delta
            call ummdp_vegter_yieldfunc(3,ss,seba,dseds,d2seds2,0,
     1                                    pryld,ndyld)
            ss(j)=s0(j)-delta
            ss(k)=s0(k)+delta
            call ummdp_vegter_yieldfunc(3,ss,seab,dseds,d2seds2,0,
     1                                    pryld,ndyld)
            ss(j)=s0(j)+delta
            ss(k)=s0(k)+delta
            call ummdp_vegter_yieldfunc(3,ss,sebb,dseds,d2seds2,0,
     1                                    pryld,ndyld)
            ss(j)=s0(j)
            ss(k)=s0(k)
            a=(seba-seaa)/(2.0d0*delta)
            b=(sebb-seab)/(2.0d0*delta)
            d2seds2(j,k)=(b-a)/(2.0d0*delta)
          end if
        end do
      end do
c
      return
      end subroutine ummdp_vegter_d2seds2n
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     calc. equivalent stress for d2seds2n
c
      subroutine ummdp_vegter_yieldfunc ( nttl,s,se,dseds,d2seds2,
     1                                    nreq,pryld,ndyld )
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: nttl,nreq,ndyld
      real*8 ,intent(in) :: s(3),pryld(ndyld)
c
      real*8,intent(out) :: se
      real*8,intent(out) :: dseds(3)
      real*8,intent(out) :: d2seds2(3,3)
c-----------------------------------------------------------------------
c
      call ummdp_vegter ( s,se,dseds,d2seds2,nreq,pryld,ndyld )
c
      return
      end subroutine ummdp_vegter_yieldfunc
c
c
c************************************************************************
c
c     YLD2000-2D YIELD FUNCTION AND DERIVATIVES
c
c       doi: https://doi.org/10.1016/S0749-6419(02)00019-0
c
      subroutine ummdp_yld2000 ( s,se,dseds,d2seds2,nreq,pryld,ndyld )
c     
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
      integer i,j,k,l,m,n,nd,nd1,nd2
      real*8 em,q
      real*8 a(8),phi(2),dsedphi(2),dphidx(2,2)
      real*8 x(2,2),y(2,3),d2sedphi2(2,2)
      real*8 am(2,3,3),dxdy(2,2,3),dyds(2,3,3),d2phidx2(2,2,2)
      real*8 d2xdy2(2,2,3,3)
c-----------------------------------------------------------------------         
c
c     >>> Arguments List
c
c     s        | stress tensor                                      (in)
c     pryld    | yield function properties                          (in)
c     ndyld    | number of yield function properties                (in)
c     nreq     | flag for required outputs                          (in)
c
c     se       | equivalent stress                                 (out)
c     dseds    | 1st order derivative of yield function wrt stress (out)
c     d2seds2  | 2nd order derivative of yield function wrt stress (out)
c
c     >>> Local Variables List
c
c     x(1,i)     | X'i   (i=1~2)
c     x(2,i)     | X"i   (i=1~2)
c     y(1,i)     | X'xx,X'yy,X'xy (i=1~3)
c     y(2,i)     | X"xx,X"yy,X"xy (i=1~3)
c     phi(1)     | phi'
c     phi(2)     | phi"
c     
c     am      | linear transformation matrix for stress tensor
c     a       | anisotropic parameters
c
c-----------------------------------------------------------------------  
c
c                                            ---- anisotropic parameters
      do i = 1,8
        a(i) = pryld(i+1)
      end do
      em = pryld(9+1)
c                                  ---- set linear transformation matrix
      call ummdp_yld2000_2d_am ( a,am )
c                                                 ---- equivalent stress
      call ummdp_yld2000_2d_xyphi ( s,em,am,x,y,phi )
      q = phi(1) + phi(2)
      if ( q <= 0.0 ) q = 0.0
      se = (0.5d0*q) ** (1.0d0/em)
c                                              ---- 1st order derivative
      if ( nreq >= 1 ) then
        call ummdp_yld2000_2d_ds1 ( em,am,x,y,phi,dsedphi,dphidx,dxdy,
     1                              dyds,se )
        call ummdp_utility_clear1 ( dseds,3 )
        do nd = 1,2
          do m = 1,2
            do k = 1,3
              do i = 1,3
                dseds(i) = dseds(i) + (dsedphi(nd)*dphidx(nd,m)
     1                                 * dxdy(nd,m,k)*dyds(nd,k,i))
              end do
            end do
          end do
        end do
      end if
c                                              ---- 2nd order derivative
      if ( nreq >= 2 ) then
        call ummdp_yld2000_2d_ds2 ( phi,x,y,em,d2sedphi2,d2phidx2,
     1                              d2xdy2,se )                   
        call ummdp_utility_clear2 ( d2seds2,3,3 )
        do i = 1,3
        do j = 1,3
          do nd1 = 1,2
          do nd2 = 1,2
            do k = 1,2
            do l = 1,2
              do m = 1,3
              do n = 1,3
                d2seds2(i,j) = d2seds2(i,j) 
     1                         + (d2sedphi2(nd1,nd2)*dphidx(nd1,k)                     
     2                            * dxdy(nd1,k,m)*dyds(nd1,m,i)
     3                            * dphidx(nd2,l)*dxdy(nd2,l,n)
     4                            * dyds(nd2,n,j))
              end do
              end do
            end do
            end do
          end do
          end do
          do nd = 1,2
            do k = 1,2
            do l = 1,2
              do m = 1,3
              do n = 1,3
                d2seds2(i,j) = d2seds2(i,j) 
     1                          + (dsedphi(nd)*d2phidx2(nd,k,l)
     2                             * dxdy(nd,k,m)*dyds(nd,m,i)
     3                             * dxdy(nd,l,n)*dyds(nd,n,j))
              end do
              end do
            end do
            end do
          end do
          do nd = 1,2
            do k = 1,2
              do m = 1,3
              do n = 1,3
                d2seds2(i,j) = d2seds2(i,j) 
     1                         + (dsedphi(nd)*dphidx(nd,k)
     2                            * d2xdy2(nd,k,m,n)*dyds(nd,m,i)
     3                            * dyds(nd,n,j))
              end do
              end do
            end do
          end do
        end do
        end do
      end if
c
      return
      end subroutine ummdp_yld2000
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     SET LINEAR TRANSFORMATION MATRIX
c
      subroutine ummdp_yld2000_2d_am ( a,am )
c
c-----------------------------------------------------------------------
      implicit none
c
      real*8,intent(in) :: a(8)
c
      real*8,intent(out) :: am(2,3,3)
c
      integer i,j
c-----------------------------------------------------------------------
c
c                                      ---- linear transformation matrix
      am(1,1,1) =  2.0d0*a(1)
      am(1,1,2) = -1.0d0*a(1)
      am(1,1,3) =  0.0d0
      am(1,2,1) = -1.0d0*a(2)
      am(1,2,2) =  2.0d0*a(2)
      am(1,2,3) =  0.0d0
      am(1,3,1) =  0.0d0
      am(1,3,2) =  0.0d0
      am(1,3,3) =  3.0d0*a(7)
c
      am(2,1,1) = -2.0d0*a(3) + 2.0d0*a(4) + 8.0d0*a(5) - 2.0d0*a(6)
      am(2,1,2) =        a(3) - 4.0d0*a(4) - 4.0d0*a(5) + 4.0d0*a(6)
      am(2,1,3) =  0.0d0
      am(2,2,1) =  4.0d0*a(3) - 4.0d0*a(4) - 4.0d0*a(5) +       a(6)
      am(2,2,2) = -2.0d0*a(3) + 8.0d0*a(4) + 2.0d0*a(5) - 2.0d0*a(6)
      am(2,2,3) =  0.0d0
      am(2,3,1) =  0.0d0
      am(2,3,2) =  0.0d0
      am(2,3,3) =  9.0d0*a(8)
c
      do i = 1,3
        do j = 1,3
          am(1,i,j) = am(1,i,j) / 3.0d0
          am(2,i,j) = am(2,i,j) / 9.0d0
        end do
      end do
c
      return
      end subroutine ummdp_yld2000_2d_am
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     CALCULATE barlat-yld2k function x,y,phi
c
      subroutine ummdp_yld2000_2d_xyphi ( s,em,am,x,y,phi )
c
c-----------------------------------------------------------------------
      implicit none
c
      real*8,intent(in) :: em
      real*8,intent(in) :: s(3)
      real*8,intent(in) :: am(2,3,3)
c
      real*8,intent(out) :: phi(2)
      real*8,intent(out) :: x(2,2),y(2,3)
c
      integer i,j,nd
      real*8 a
      real*8 p(2)
c-----------------------------------------------------------------------
c
      p(1) =  1.0d0
      p(2) = -1.0d0
c                                                       ---- {y}=[am]{s}
      call ummdp_utility_clear2 ( y,2,3 )
      do nd = 1,2
        do i = 1,3
          do j = 1,3
            y(nd,i) = y(nd,i) + am(nd,i,j)*s(j)
          end do
        end do
      end do
c                                        ---- {x}=principle value of {y}
      do nd = 1,2
        a = (y(nd,1)-y(nd,2))**2.d0 + 4.0d0*y(nd,3)**2.d0
        a = sqrt(a)
        do i = 1,2
          x(nd,i) = 0.5d0 * (y(nd,1)+y(nd,2)+p(i)*a)
        end do
      end do
c                                                 ---- phi(1) and phi(2)
      nd = 1
      phi(nd) = abs(x(nd,1)-x(nd,2))**em
      nd = 2
      phi(nd) = abs(2.0d0*x(nd,2)+x(nd,1))**em
     1          + abs(2.0d0*x(nd,1)+x(nd,2))**em
c
      return
      end subroutine ummdp_yld2000_2d_xyphi
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     SET 1ST ORDER DERIVATIVE OF PARAMETERS
c
      subroutine ummdp_yld2000_2d_ds1 ( em,am,x,y,phi,dsedphi,dphidx,
     1                                  dxdy,dyds,se )     
c  
c-----------------------------------------------------------------------
      implicit none
c
      real*8,intent(in) :: em,se
      real*8,intent(in) :: phi(2)
      real*8,intent(in) :: x(2,2),y(2,3)
      real*8,intent(in) :: am(2,3,3)
c
      real*8,intent(out) :: dsedphi(2)
      real*8,intent(out) :: dphidx(2,2)
      real*8,intent(out) :: dxdy(2,2,3),dyds(2,3,3)
c
      integer i,j,nd
      real*8 eps,emi,q,a,a0,a1,a2,b0,b1,b2,sgn0,sgn1,sgn2
      real*8 p(2)
c-----------------------------------------------------------------------
c
      eps = 1.0d-16
c
      p(1) =  1.0d0
      p(2) = -1.0d0
      emi = 1.0d0 / em
c                                                          ---- dse/dphi
      q = phi(1) + phi(2)
      if ( q <= 0.0 ) q = 0.0
      do nd = 1,2
        dsedphi(nd) = (0.5d0**emi) * emi * q**(emi-1.0d0)
      end do
c                                                           ---- dphi/dx
      nd = 1
      a0 = x(nd,1) - x(nd,2)
      b0 = abs(a0)
      sgn0 = 0
      if ( b0 >= eps*se ) sgn0 = a0 / b0
      dphidx(nd,1) =  em * b0**(em-1.0d0) * sgn0
      dphidx(nd,2) = -em * b0**(em-1.0d0) * sgn0
c
      nd = 2
      a1 = 2.0d0*x(nd,1) +       x(nd,2)
      a2 =       x(nd,1) + 2.0d0*x(nd,2)
      b1 = abs(a1)
      b2 = abs(a2)
      sgn1 = 0.0d0
      sgn2 = 0.0d0
      if ( b1 >= eps*se ) sgn1 = a1 / b1
      if ( b2 >= eps*se ) sgn2 = a2 / b2
      dphidx(nd,1) = em*( 2.0d0*b1**(em-1.0d0)*sgn1
     1                    + b2**(em-1.0d0)*sgn2 )
      dphidx(nd,2) = em*( b1**(em-1.0d0)*sgn1
     1                    + 2.0d0*b2**(em-1.0d0)*sgn2 )
c
      do nd = 1,2
        a = (y(nd,1)-y(nd,2))*(y(nd,1)-y(nd,2)) + 4.0d0*y(nd,3)*y(nd,3)
        a = sqrt(a)
        if ( a > eps*se ) then
          do j = 1,2
            dxdy(nd,j,1) = 0.5d0 * (1.0d0+p(j)*(y(nd,1)-y(nd,2))/a)
            dxdy(nd,j,2) = 0.5d0 * (1.0d0-p(j)*(y(nd,1)-y(nd,2))/a)
            dxdy(nd,j,3) = 2.0d0 *        p(j)* y(nd,3)         /a
          end do
        else
          do j = 1,2
            dxdy(nd,j,1) = 0.5d0 * (1.0d0+0.0d0)
            dxdy(nd,j,2) = 0.5d0 * (1.0d0-0.0d0)
            dxdy(nd,j,3) = 2.0d0 *        0.0d0
          end do
        end if
      end do
c
      do nd = 1,2
        do i = 1,3
          do j = 1,3
            dyds(nd,i,j) = am(nd,i,j)
          end do
        end do
      end do
c
      return
      end subroutine ummdp_yld2000_2d_ds1
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     SET 2ND ORDER DERIVATIVE OF PARAMETERS
c
      subroutine ummdp_yld2000_2d_ds2 ( phi,x,y,em,d2sedphi2,d2phidx2,
     1                                  d2xdy2,se )
c
c-----------------------------------------------------------------------
      implicit none
c
      real*8,intent(in) :: em,se
      real*8,intent(in) :: phi(2)
      real*8,intent(in) :: x(2,2),y(2,3)
c
      real*8,intent(out) :: d2sedphi2(2,2)
      real*8,intent(out) :: d2phidx2(2,2,2)
      real*8,intent(out) :: d2xdy2(2,2,3,3)
c
      integer i,j,m,nd,nd1,nd2,ij
      real*8 eps,emi,q,a
      real*8 p(2)
c-----------------------------------------------------------------------
c
      eps = 1.0d-16
c
      p(1) =  1.0d0
      p(2) = -1.0d0
      emi = 1.0d0 / em
c                                                        ---- d2se/dphi2
      q = phi(1) + phi(2)
      if ( q <= 0.0d0 ) q = 0.0d0
      do nd1 = 1,2
        do nd2 = 1,2
          a = 0.5d0**emi * emi * (emi-1.0d0) * q**(emi-2.0d0)
          d2sedphi2(nd1,nd2) = a
        end do
      end do
c                                                         ---- d2phi/dx2
      nd = 1
      do i = 1,2
        do j = 1,2
          a = (em-1.0d0) * em * (abs(x(nd,1)-x(nd,2)))**(em-2.0d0)
          if ( i /= j ) a = -a
          d2phidx2(nd,i,j) = a
        end do
      end do
      nd = 2
      do i = 1,2
        do j = 1,2
          if ( i == j ) then
            if ( i == 1 ) then
              a = (em-1.0d0) * em
     1            * (4.0d0*(abs(2.0d0*x(nd,1)+x(nd,2)))**(em-2.0d0)
     2               + (abs(x(nd,1)+2.0d0*x(nd,2)))**(em-2.0d0))
            else
              a = (em-1.0d0) * em
     1             * ((abs(2.0d0*x(nd,1)+x(nd,2)))**(em-2.0d0)
     2                + 4.0d0*(abs(x(nd,1)+2.0d0*x(nd,2)))**(em-2.0d0))
            end if
          else
            a = (em-1.0d0) * em
     1           * (2.0d0*(abs(2.0d0*x(nd,1)+x(nd,2)))**(em-2.0d0)
     2              + 2.0d0*(abs(x(nd,1)+2.0d0*x(nd,2)))**(em-2.0d0))
          end if
          d2phidx2(nd,i,j) = a
        end do
      end do
c                                                           ---- d2x/dy2
      do nd = 1,2
        a = (y(nd,1)-y(nd,2))*(y(nd,1)-y(nd,2)) + 4.0d0*y(nd,3)*y(nd,3)       
        if ( a > eps*se ) then
          a = 1.0d0 / sqrt(a**3)
          do m = 1,2
            do i = 1,3
              do j = 1,3
                ij = i*10+j
                if ( (ij == 11) .or. (ij == 22) ) then
                  q = y(nd,3) * y(nd,3)
                else if ( ij == 33 ) then
                  q = (y(nd,1)-y(nd,2)) * (y(nd,1)-y(nd,2))
                else if ( (ij == 12) .or. (ij == 21) ) then
                  q = -y(nd,3) * y(nd,3)
                else if ( (ij == 23) .or. (ij == 32) ) then
                  q = y(nd,3) * (y(nd,1)-y(nd,2))
                else
                  q = -y(nd,3) * (y(nd,1)-y(nd,2))
                end if
                d2xdy2(nd,m,i,j) = 2.0d0 * a * p(m) * q
              end do
            end do
          end do
        else
          do m = 1,2
            do i = 1,3
              do j = 1,3
                d2xdy2(nd,m,i,j) = 0.0d0
              end do
            end do
          end do
        end if
      end do
c
      return
      end subroutine ummdp_yld2000_2d_ds2
c
c
c
************************************************************************
c     YLD2004-18p YIELD FUNCTION AND DERIVATIVES
c
c       doi: 
c
      subroutine ummdp_yld2004_18p ( s,se,dseds,d2seds2,nreq,
     1                                pryld,ndyld )
c-----------------------------------------------------------------------
      implicit none
c
      integer nreq,ndyld
      real*8 se
      real*8 s(6),dseds(6),pryld(ndyld)
      real*8 d2seds2(6,6)
c
      integer i,j,k,l,m,n,ip,iq,ir
      real*8 am,ami,dc,pi,eps2,eps3,del,cetpq1,cetpq2,fai,dsedfa,
     1       d2sedfa2,dummy
      real*8 sp1(6),sp2(6),psp1(3),psp2(3),hp1(3),hp2(3),
     1       dfadpsp1(3),dfadpsp2(3),dfadhp1(3),dfadhp2(3),dfads(6)
      real*8 cp1(6,6),cp2(6,6),cl(6,6),ctp1(6,6),ctp2(6,6),
     1       dpsdhp1(3,3),dpsdhp2(3,3),dhdsp1(3,6),dhdsp2(3,6),
     2       dsdsp1(6,6),dsdsp2(6,6),d2fadpsp11(3,3),d2fadpsp22(3,3),
     3       d2fadpsp12(3,3),d2fadpsp21(3,3),d2fadhp11(3,3),
     4       d2fadhp22(3,3),d2fadhp12(3,3),d2fadhp21(3,3),d2fads2(6,6),
     5       delta(3,3),xx1(3,6),xx2(3,6)
      real*8 d2psdhp11(3,3,3),d2psdhp22(3,3,3),d2hdsp11(3,6,6),
     1       d2hdsp22(3,6,6)
c-----------------------------------------------------------------------
c
c     sp1(6),sp2(6)     : linear-transformed stress
c     cp1(6,6),cp2(6,6) : matrix for anisotropic parameters
c     cl(6,6)           : matrix for transforming Cauchy stress to deviator
c     ctp1(6,6)         : matrix for transforming Cauchy stress to sp1
c     ctp2(6,6)         : matrix for transforming Cauchy stress to sp2
c     dc                : coef. of equivalent stress (se=(fai/dc)**(1/am))
c     aap1(3),aap2(3)   : for dc
c     ppp1,ppp2         : for dc
c     qqp1,qqp2         : for dc
c     ttp1,ttp2         : for dc
c     bbp1(3),bbp2(3)   : for dc
c     fai               : yield fuction
c     psp1(3)           : principal values of sp1
c     psp2(3)           : principal values of sp2
c     hp1(3)            : invariants of sp1
c     hp2(3)            : invariants of sp2
c     cep1,cep2         : coef. of characteristic equation p
c     ceq1,ceq2         : coef. of characteristic equation q
c     cet1,cet2         : arccos ( q/p^(3/2) )
c
c     dfadpsp1(3),dfadpsp2(3)   : d(fai)/d(psp)
c     dfadhp1(3),dfadhp2(3)     : d(fai)/d(hp)
c     dpsdhp1(3,3),dpsdhp2(3,3) : d(psp)/d(hp)
c     dhdsp1(3,6),dhdsp2(3,6)   : d(hp)/d(sp)
c     dsdsp1(6,6), dsdsp2(6,6)  : d(sp)/d(s)
c     dfads(6)                  : d(fai)/d(s)
c
c     d2fadpsp11(3,3),d2fadpsp22(3,3)   : d2(fai)/d(psp)2
c     d2fadpsp12(3,3),d2fadpsp21(3,3)   : d2(fai)/d(psp)2
c     d2psdhp11(3,3,3),d2psdhp22(3,3,3) : d2(psp)/d(hp)2
c     d2hdsp11(3,6,6),d2hdsp22(3,6,6)   : d2(hp)/d(sp)2
c     d2fads2(6,6)                      : d2(fai)/d(s)2
c
c     eps2,eps3 : values for calculation of d(psp)/d(hp),d2(psp)/d(hp)2
c
      pi = acos(-1.0d0)
      eps2 = 1.0d-15
      eps3 = 1.0d-8
      del = 1.0d-4
c                                                   ---- Kronecker Delta
      call ummdp_utility_clear2 ( delta,3,3 )
      do i = 1,3
        delta(i,i) = 1.0d0
      end do
c                                        ---- set anisotropic parameters
      call ummdp_utility_clear2 ( cp1,6,6 )
      call ummdp_utility_clear2 ( cp2,6,6 )
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
      dc       =  4
      ami = 1.0d0 / am
c
c             ---- set matrix for transforming Cauchy stress to deviator
      call ummdp_utility_clear2 ( cl,6,6 )
      do i = 1,3
        do j = 1,3
          if ( i == j ) then
            cl(i,j) = 2.0d0
          else
            cl(i,j) = -1.0d0
          end if
        end do
      end do
      do i = 4,6
        cl(i,i) = 3.0d0
      end do
      do i = 1,6
        do j = 1,6
          cl(i,j) = cl(i,j) / 3.0d0
        end do
      end do
c
c                  ---- matrix for transforming Cauchy stress to sp1,sp2
      call ummdp_utility_mm ( ctp1,cp1,cl,6,6,6 )
      call ummdp_utility_mm ( ctp2,cp2,cl,6,6,6 )
c                               ---- coefficient of equivalent stress dc
c      call ummdp_yld2004_18p_coef ( cp1,cp2,pi,am,dc )
c                                  ---- calculation of equivalent stress
      call ummdp_yld2004_18p_yf ( ctp1,ctp2,s,am,ami,dc,pi,
     &                             sp1,sp2,psp1,psp2,hp1,hp2,
     &                             cetpq1,cetpq2,fai,se )
c
c                                            ---- 1st order differential
      if ( nreq >= 1 ) then
c                                                     ---- d(fai)/d(psp)
        do i = 1,3
          dfadpsp1(i) = 0.0d0
          dfadpsp2(i) = 0.0d0
          do j = 1,3
            dfadpsp1(i) = dfadpsp1(i) + (psp1(i)-psp2(j)) *
     &                    abs(psp1(i)-psp2(j))**(am-2.0d0)
            dfadpsp2(i) = dfadpsp2(i) + (psp1(j)-psp2(i)) *
     &                    abs(psp1(j)-psp2(i))**(am-2.0d0)
          end do
          dfadpsp1(i) = dfadpsp1(i) * am
          dfadpsp2(i) = dfadpsp2(i) * (-am)
        end do
c                                         ---- d(psp)/d(hp)&d(fai)/d(hp)
        call ummdp_utility_clear2 ( dpsdhp1,3,3 )
        call ummdp_utility_clear2 ( dpsdhp2,3,3 )
        call ummdp_utility_clear1 ( dfadhp1,3 )
        call ummdp_utility_clear1 ( dfadhp2,3 )
c                                                  ---- theta'<>0 & <>pi
        if ( abs(cetpq1-1.0d0) >= eps2 .and. 
     &       abs(cetpq1+1.0d0) >= eps2 ) then
          do i = 1,3
            call ummdp_yld2004_18p_dpsdhp ( i,psp1,hp1,dpsdhp1 )
          end do
c                                                          ---- theta'=0
        else if ( abs(cetpq1-1.0d0) < eps2 ) then
          i = 1
          call ummdp_yld2004_18p_dpsdhp ( i,psp1,hp1,dpsdhp1 )
          do i = 2,3
            do j = 1,3
              dpsdhp1(i,j) = -0.5d0 * (dpsdhp1(1,j)-3.0d0*delta(1,j))
            end do
           end do
c                                                         ---- theta'=pi
        else
          i = 3
          call ummdp_yld2004_18p_dpsdhp ( i,psp1,hp1,dpsdhp1 )
          do i = 1,2
            do j = 1,3
              dpsdhp1(i,j) = -0.5d0 * (dpsdhp1(3,j)-3.0d0*delta(1,j))
            end do
           end do
        end if
c                                                 ---- theta''<>0 & <>pi
        if ( abs(cetpq2-1.0d0) >= eps2 .and. 
     &       abs(cetpq2+1.0d0) >= eps2 ) then
          do i = 1,3
            call ummdp_yld2004_18p_dpsdhp ( i,psp2,hp2,dpsdhp2 ) 
          end do
c                                                         ---- theta''=0
        else if ( abs(cetpq2-1.0d0) < eps2 ) then
          i = 1
          call ummdp_yld2004_18p_dpsdhp ( i,psp2,hp2,dpsdhp2 )
          do i = 2,3
            do j = 1,3
              dpsdhp2(i,j) = -0.5d0 * (dpsdhp2(1,j)-3.0d0*delta(1,j))
            end do
           end do
c                                                        ---- theta''=pi
        else
          i = 3
          call ummdp_yld2004_18p_dpsdhp ( i,psp2,hp2,dpsdhp2 )
          do i = 1,2
            do j = 1,3
              dpsdhp2(i,j) = -0.5d0 * (dpsdhp2(3,j)-3.0d0*delta(1,j))
            end do
           end do
        end if
c
        do i = 1,3
          do j = 1,3
            dfadhp1(i) = dfadhp1(i) + dfadpsp1(j)*dpsdhp1(j,i)
            dfadhp2(i) = dfadhp2(i) + dfadpsp2(j)*dpsdhp2(j,i)
          end do
        end do
c                                                       ---- d(hp)/d(sp)
        call ummdp_utility_clear2 ( dhdsp1,3,6 )
        call ummdp_utility_clear2 ( dhdsp2,3,6 )
        do i = 1,3
          j = mod(i,  3) + 1
          k = mod(i+1,3) + 1
          l = mod(i,  3) + 4
          if ( i == 1 ) l = 6
          if ( i == 2 ) l = 5
          dhdsp1(1,i) =  1.0d0/3.0d0
          dhdsp2(1,i) =  1.0d0/3.0d0
          dhdsp1(2,i) = -1.0d0/3.0d0 * (sp1(j)+sp1(k))
          dhdsp2(2,i) = -1.0d0/3.0d0 * (sp2(j)+sp2(k))
          dhdsp1(3,i) =  1.0d0/2.0d0 * (sp1(j)*sp1(k)-sp1(l)**2)
          dhdsp2(3,i) =  1.0d0/2.0d0 * (sp2(j)*sp2(k)-sp2(l)**2)
        end do
        do i = 4,6
          k = mod(i+1,3) + 1
          l = mod(i,  3) + 4
          m = mod(i+1,3) + 4
          if ( i == 5 ) k = 2
          if ( i == 6 ) k = 1
          dhdsp1(2,i) = 2.0d0/3.0d0 * sp1(i)
          dhdsp2(2,i) = 2.0d0/3.0d0 * sp2(i)
          dhdsp1(3,i) = sp1(l)*sp1(m) - sp1(k)*sp1(i)
          dhdsp2(3,i) = sp2(l)*sp2(m) - sp2(k)*sp2(i)
        end do
c                                                        ---- d(sp)/d(s)
        do i = 1,6
          do j = 1,6
            dsdsp1(i,j) = ctp1(i,j)
            dsdsp2(i,j) = ctp2(i,j)
          end do
        end do
c                                                       ---- d(fai)/d(s)
        call ummdp_utility_clear1 ( dfads,6 )
        call ummdp_utility_clear2 ( xx1,3,6 )
        call ummdp_utility_clear2 ( xx2,3,6 )
        do l = 1,6
          do j = 1,3
            do k = 1,6
              xx1(j,l) = xx1(j,l) + dhdsp1(j,k)*dsdsp1(k,l)
              xx2(j,l) = xx2(j,l) + dhdsp2(j,k)*dsdsp2(k,l)
            end do
            dfads(l) = dfads(l) +
     &                 dfadhp1(j)*xx1(j,l) + 
     &                 dfadhp2(j)*xx2(j,l)
          end do
        end do
c                                                        ---- d(se)/d(s)
        dsedfa = fai**(ami-1.0d0) / am / dc**ami
        do i = 1,6
          dseds(i) = dsedfa*dfads(i)
        end do
c
      endif
c                                            ---- 2nd order differential
      if ( nreq >= 2 ) then
c                                                   ---- d2(fai)/d(psp)2
        call ummdp_utility_clear2 ( d2fadpsp11,3,3 )
        call ummdp_utility_clear2 ( d2fadpsp22,3,3 )
        call ummdp_utility_clear2 ( d2fadpsp12,3,3 )
        call ummdp_utility_clear2 ( d2fadpsp21,3,3 )
        do i = 1,3
          d2fadpsp11(i,i) = am*(am-1.0d0) *
     1                      ( abs(psp1(i)-psp2(1))**(am-2.0d0) +
     2                        abs(psp1(i)-psp2(2))**(am-2.0d0) +
     3                        abs(psp1(i)-psp2(3))**(am-2.0d0) )
          d2fadpsp22(i,i) = am*(am-1.0d0) * 
     1                      ( abs(psp1(1)-psp2(i))**(am-2.0d0) +
     2                        abs(psp1(2)-psp2(i))**(am-2.0d0) +
     3                        abs(psp1(3)-psp2(i))**(am-2.0d0) )
          do j = 1,3
            d2fadpsp12(i,j) = -am*(am-1.0d0) * 
     1                         abs(psp1(i)-psp2(j))**(am-2.0d0)
            d2fadpsp21(i,j) = -am*(am-1.0d0) *
     1                         abs(psp1(j)-psp2(i))**(am-2.0d0)
          end do
        end do
c                                                    ---- d2(psp)/d(hp)2
        call ummdp_utility_clear3 ( d2psdhp11,3,3,3 )
        call ummdp_utility_clear3 ( d2psdhp22,3,3,3 )
c
        if ( abs(cetpq1-1.0d0) >= eps3 .and. 
     1       abs(cetpq2-1.0d0) >= eps3 .and.
     2       abs(cetpq1+1.0d0) >= eps3 .and. 
     3       abs(cetpq2+1.0d0) >= eps3 ) then
c
          do i = 1,3
            call ummdp_yld2004_18p_d2psdhp ( i,psp1,hp1,d2psdhp11 )
            call ummdp_yld2004_18p_d2psdhp ( i,psp2,hp2,d2psdhp22 )
          end do
c                                                    ---- d2(fai)/d(hp)2
          call ummdp_utility_clear2 ( d2fadhp11,3,3 )
          call ummdp_utility_clear2 ( d2fadhp12,3,3 )
          call ummdp_utility_clear2 ( d2fadhp21,3,3 )
          call ummdp_utility_clear2 ( d2fadhp22,3,3 )
c                                                  -- d2(fai)/d(hd)d(hd)
          do iq = 1,3
            do m = 1,3
              do ip = 1,3
                do l = 1,3
                  d2fadhp11(iq,m) = d2fadhp11(iq,m) + d2fadpsp11(ip,l) *
     1                              dpsdhp1(l,m) * dpsdhp1(ip,iq)
                end do
                d2fadhp11(iq,m) = d2fadhp11(iq,m) + 
     1                            dfadpsp1(ip)*d2psdhp11(ip,iq,m)
              end do
            end do
          end do
c                                              ---- d2(fai)/d(hdd)d(hdd)
          do iq = 1,3
            do m = 1,3
              do ip = 1,3
                do l = 1,3
                  d2fadhp22(iq,m) = d2fadhp22(iq,m) + (
     1                              d2fadpsp22(ip,l) * 
     2                              dpsdhp2(l,m) * dpsdhp2(ip,iq) )
                end do
                d2fadhp22(iq,m) = d2fadhp22(iq,m) + 
     1                            dfadpsp2(ip)*d2psdhp22(ip,iq,m)
              end do
            end do
          end do
c                         ---- d2(fai)/d(hdd)d(hd) & d2(fai)/d(hd)d(hdd)
          do iq = 1,3
            do m = 1,3
              do ip = 1,3
                do l = 1,3
                  d2fadhp12(iq,m) = d2fadhp12(iq,m) + (
     1                              d2fadpsp12(ip,l)*
     2                              dpsdhp2(l,m)*dpsdhp1(ip,iq) )
                  d2fadhp21(iq,m) = d2fadhp21(iq,m) + 
     1                              d2fadpsp21(ip,l)*
     2                              dpsdhp1(l,m)*dpsdhp2(ip,iq)
                end do
               end do
            end do
          end do
c                                                     ---- d2(hp)/d(sp)2
          call ummdp_utility_clear3 ( d2hdsp11,3,6,6 )
          call ummdp_utility_clear3 ( d2hdsp22,3,6,6 )
          do i = 1,3
            j = mod(i,  3) + 1
            k = mod(i+1,3) + 1
            dummy = -1.0d0 / 3.0d0
            d2hdsp11(2,i,j) = dummy
            d2hdsp11(2,j,i) = dummy
            d2hdsp22(2,i,j) = dummy
            d2hdsp22(2,j,i) = dummy
            d2hdsp11(3,i,j) = sp1(k) / 2.0d0
            d2hdsp11(3,j,i) = sp1(k) / 2.0d0
            d2hdsp22(3,i,j) = sp2(k) / 2.0d0
            d2hdsp22(3,j,i) = sp2(k) / 2.0d0
          end do
          do i = 4,6
            j = mod(i,  3) + 4
            k = mod(i+1,3) + 1
            m = mod(i+1,3) + 4
            if ( i == 5 ) k = 2
            if ( i == 6 ) k = 1
            dummy = 2.0d0 / 3.0d0
            d2hdsp11(2,i,i) = dummy
            d2hdsp22(2,i,i) = dummy
            d2hdsp11(3,i,i) = -sp1(k)
            d2hdsp22(3,i,i) = -sp2(k)
            d2hdsp11(3,i,j) =  sp1(m)
            d2hdsp11(3,j,i) =  sp1(m)
            d2hdsp22(3,i,j) =  sp2(m)
            d2hdsp22(3,j,i) =  sp2(m)
          end do
          do i = 1,3
            j = mod(i,3) + 4
            if ( i == 1 ) j = 6
            if ( i == 2 ) j = 5
            d2hdsp11(3,i,j) = -sp1(j)
            d2hdsp11(3,j,i) = -sp1(j)
            d2hdsp22(3,i,j) = -sp2(j)
            d2hdsp22(3,j,i) = -sp2(j)
          end do
c                                                     ---- d2(fai)/d(s)2
          call ummdp_utility_clear2 ( d2fads2,6,6 )
          call ummdp_utility_clear2 ( xx1,3,6 )
          call ummdp_utility_clear2 ( xx2,3,6 )
          do i = 1,3
            do j = 1,6
              do ip = 1,6
                xx1(i,j) = xx1(i,j) + dhdsp1(i,ip)*dsdsp1(ip,j)
                xx2(i,j) = xx2(i,j) + dhdsp2(i,ip)*dsdsp2(ip,j)
              end do
            end do
          end do
          do i = 1,6
            do j = 1,6
              do iq = 1,3
                do m = 1,3
                  d2fads2(i,j) = d2fads2(i,j) +
     1                           d2fadhp11(iq,m)*xx1(iq,i)*xx1(m,j) +
     2                           d2fadhp12(iq,m)*xx1(iq,i)*xx2(m,j) +
     3                           d2fadhp21(iq,m)*xx2(iq,i)*xx1(m,j) +
     4                           d2fadhp22(iq,m)*xx2(iq,i)*xx2(m,j)
                end do
                do n = 1,6
                  do ir = 1,6
                    d2fads2(i,j) = d2fads2(i,j) +
     1                             d2hdsp11(iq,ir,n)*dfadhp1(iq)*
     2                             dsdsp1(ir,i)*dsdsp1(n,j) +
     3                             d2hdsp22(iq,ir,n)*dfadhp2(iq)*
     4                             dsdsp2(ir,i)*dsdsp2(n,j)
                  end do
                end do
              end do
            end do
          end do
c                                                      ---- d2(se)/d(s)2
          d2sedfa2 = ami * (ami-1.0d0) * fai**(ami-2.0d0) / dc**ami
          do i = 1,6
            do j = 1,6
              d2seds2(i,j) = d2sedfa2*dfads(i)*dfads(j) +
     1                       dsedfa*d2fads2(i,j)
            end do
          end do
        else
c                                            ---- numerical differential
          call ummdp_yld2004_18p_nu2 ( ctp1,ctp2,s,se,am,ami,
     1                                  dc,pi,del,d2seds2 )
        end if
      end if
c
      return
      end
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CALCULATE COEFFICIENT OF EQUIVALENT STRESS 1
c
      subroutine ummdp_yld2004_18p_coef ( cp1,cp2,pi,am,dc )
c-----------------------------------------------------------------------
      implicit none
c
      real*8 pi,am,dc
      real*8 cp1(6,6),cp2(6,6)
c
      integer i,j
      real*8 bbp1(3),bbp2(3)
c-----------------------------------------------------------------------
c
      call ummdp_yld2004_18p_coef_sub ( cp1,pi,bbp1 )
      call ummdp_yld2004_18p_coef_sub ( cp2,pi,bbp2 )
      dc = 0.0d0
      do i = 1,3
        do j = 1,3
          dc = dc + (abs(bbp1(i)-bbp2(j)))**am
        end do
      end do
c
      return
      end subroutine ummdp_yld2004_18p_coef
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CALCULATE COEFFICIENT OF EQUIVALENT STRESS 2
c
      subroutine ummdp_yld2004_18p_coef_sub ( cp,pi,bbp )
c-----------------------------------------------------------------------
      implicit none
c
      real*8 pi
      real*8 cp(6,6)
c
      real*8 ppp,qqp,ttp
      real*8 aap(3),bbp(3)
c-----------------------------------------------------------------------
c                                                  ---- coefficients aap
      aap(1) = (cp(1,2)+cp(1,3)-2.0d0*cp(2,1) +
     1          cp(2,3)-2.0d0*cp(3,1)+cp(3,2))/9.0d0
      aap(2) = ((2.0d0*cp(2,1)-cp(2,3))*(cp(3,2)-2.0d0*cp(3,1)) +
     2          (2.0d0*cp(3,1)-cp(3,2))*(cp(1,2)+cp(1,3)) +
     3          (cp(1,2)+cp(1,3))*(2.0d0*cp(2,1)-cp(2,3)))/2.7d1
      aap(3) = (cp(1,2)+cp(1,3))*(cp(2,3)-2.0d0*cp(2,1))*
     1         (cp(3,2)-2.0d0*cp(3,1))/5.4d1
c                                          ---- coefficients ppp,qqp,ttp
      ppp = aap(1)**2 + aap(2)
      qqp = (2.0d0*aap(1)**3+3.0d0*aap(1)*aap(2)+2.0d0*aap(3)) / 2.0d0
      ttp = acos(qqp / ppp**(3.0d0/2.0d0))
c                                                  ---- coefficients bbp
      bbp(1) = 2.0d0*sqrt(ppp)*cos(ttp/3.0d0) + aap(1)
      bbp(2) = 2.0d0*sqrt(ppp)*cos((ttp+4.0d0*pi)/3.0d0) + aap(1)
      bbp(3) = 2.0d0*sqrt(ppp)*cos((ttp+2.0d0*pi)/3.0d0) + aap(1)
c
      return
      end subroutine ummdp_yld2004_18p_coef_sub
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CALCULATE YIELD FUNCTION 1
c
      subroutine ummdp_yld2004_18p_yf ( ctp1,ctp2,s,am,ami,dc,pi,
     1                                   sp1,sp2,psp1,psp2,hp1,hp2,
     2                                   cetpq1,cetpq2,fai,se )
c-----------------------------------------------------------------------
      implicit none
c
      real*8 am,ami,dc,pi,cetpq1,cetpq2,fai,se
      real*8 s(6),sp1(6),sp2(6),psp1(3),psp2(3),hp1(3),hp2(3)
      real*8 ctp1(6,6),ctp2(6,6)
c
      integer i,j
c-----------------------------------------------------------------------
      call ummdp_yld2004_18p_yfsub ( ctp1,s,pi,sp1,psp1,hp1,cetpq1 )
      call ummdp_yld2004_18p_yfsub ( ctp2,s,pi,sp2,psp2,hp2,cetpq2 )
c                                                    ---- yield function
      fai = 0.0d0
      do i = 1,3
        do j = 1,3
          fai = fai + (abs(psp1(i)-psp2(j)))**am
        end do
      end do
c                                                 ---- equivalent stress
      se = (fai/dc) ** ami
c
      return
      end subroutine ummdp_yld2004_18p_yf
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CALCULATE YIELD FUNCTION 2
c
      subroutine ummdp_yld2004_18p_yfsub ( ctp,s,pi,sp,psp,hp,cetpq )
c-----------------------------------------------------------------------
      implicit none
c
      real*8 pi,cetpq
      real*8 s(6),sp(6),psp(3),hp(3)
      real*8 ctp(6,6)
c
      integer i
      real*8 hpq,cep,ceq,cet
c-----------------------------------------------------------------------
c
c                                         ---- linear-transformed stress
      call ummdp_utility_mv ( sp,ctp,s,6,6 )
c
c                                                  ---- invariants of sp
      hp(1) = (sp(1)+sp(2)+sp(3)) / 3.0d0
      hp(2) = (sp(5)**2+sp(6)**2+sp(4)**2 - 
     1         sp(2)*sp(3)-sp(3)*sp(1)-sp(1)*sp(2)) / 3.0d0
      hp(3) = (2.0d0*sp(5)*sp(6)*sp(4)+sp(1)*sp(2)*sp(3) -
     1         sp(1)*sp(6)**2-sp(2)*sp(5)**2-sp(3)*sp(4)**2) / 2.0d0
c
c                           ---- coefficients of characteristic equation
      hpq = sqrt(hp(1)**2 + hp(2)**2 + hp(3)**2)
      if ( hpq > 1.0e-16 ) then
        cep = hp(1)**2 + hp(2)
        ceq = (2.0d0*hp(1)**3+3.0d0*hp(1)*hp(2)+2.0d0*hp(3)) / 2.0d0
        cetpq = ceq / cep**(3.0d0/2.0d0)
        if ( cetpq >  1.0d0 ) cetpq =  1.0d0
        if ( cetpq < -1.0d0 ) cetpq = -1.0d0
        cet = acos(cetpq)
c
c                                           ---- principal values of sp1
        psp(1) = 2.0d0*sqrt(cep)*cos(cet/3.0d0) + hp(1)
        psp(2) = 2.0d0*sqrt(cep)*cos((cet+4.0d0*pi)/3.0d0) + hp(1)
        psp(3) = 2.0d0*sqrt(cep)*cos((cet+2.0d0*pi)/3.0d0) + hp(1)
      else
        cetpq = 0.0
        do i = 1,3
          psp(i) = 0.0
        end do
      end if
c
      return
      end subroutine ummdp_yld2004_18p_yfsub
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     NUMERICAL DIFFERENTIATION FOR 2ND ORDER DERIVATIVES
c
      subroutine ummdp_yld2004_18p_nu2 ( ctp1,ctp2,s,se,am,ami,
     1                                    dc,pi,del,d2seds2 )
c-----------------------------------------------------------------------
      implicit none
c
      real*8 se,am,ami,dc,pi,del
      real*8 s(6)
      real*8 ctp1(6,6),ctp2(6,6),d2seds2(6,6)
c
      integer i,j
      real*8 cetpq1,cetpq2,fai,sea,seb,seaa,seba,seab,sebb,abc1,abc2
      real*8 sp1(6),sp2(6),psp1(3),psp2(3),hp1(3),hp2(3),s0(6)
c-----------------------------------------------------------------------
c
      s0(:) = s(:)
      do i = 1, 6
        do j = 1, 6
          if ( i == j ) then
            s0(i) = s(i) - del
            call ummdp_yld2004_18p_yf ( ctp1,ctp2,s0,am,ami,dc,pi,
     1                                   sp1,sp2,psp1,psp2,hp1,hp2,
     2                                   cetpq1,cetpq2,fai,sea )
            s0(i) = s(i) + del
            call ummdp_yld2004_18p_yf ( ctp1,ctp2,s0,am,ami,dc,pi,
     1                                   sp1,sp2,psp1,psp2,hp1,hp2,
     2                                   cetpq1,cetpq2,fai,seb )
            s0(i) = s(i)
            abc1 = (se-sea) / del
            abc2 = (seb-se) / del
            d2seds2(i,j) = (abc2-abc1) / del
          else
            s0(i) = s(i) - del
            s0(j) = s(j) - del
            call ummdp_yld2004_18p_yf ( ctp1,ctp2,s0,am,ami,dc,pi,
     1                                   sp1,sp2,psp1,psp2,hp1,hp2,
     2                                   cetpq1,cetpq2,fai,seaa )
            s0(i) = s(i) + del
            s0(j) = s(j) - del
            call ummdp_yld2004_18p_yf ( ctp1,ctp2,s0,am,ami,dc,pi,
     1                                   sp1,sp2,psp1,psp2,hp1,hp2,
     2                                   cetpq1,cetpq2,fai,seba )
            s0(i) = s(i) - del
            s0(j) = s(j) + del
            call ummdp_yld2004_18p_yf ( ctp1,ctp2,s0,am,ami,dc,pi,
     1                                   sp1,sp2,psp1,psp2,hp1,hp2,
     2                                   cetpq1,cetpq2,fai,seab )
            s0(i) = s(i) + del
            s0(j) = s(j) + del
            call ummdp_yld2004_18p_yf ( ctp1,ctp2,s0,am,ami,dc,pi,
     1                                   sp1,sp2,psp1,psp2,hp1,hp2,
     2                                   cetpq1,cetpq2,fai,sebb )
            s0(i) = s(i)
            s0(j) = s(j)
            abc1 = (seba-seaa) / (2.0d0*del)
            abc2 = (sebb-seab) / (2.0d0*del)
            d2seds2(i,j) = (abc2-abc1) / (2.0d0*del)
          end if
        end do
      end do
c
      return
      end subroutine ummdp_yld2004_18p_nu2
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CALCULATE d(psp)/d(hp)
c
      subroutine ummdp_yld2004_18p_dpsdhp ( i,psp,hp,dpsdhp )
c-----------------------------------------------------------------------
      implicit none
c
      integer i
      real*8 psp(3),hp(3)
      real*8 dpsdhp(3,3)
c
      real*8 dummy
c-----------------------------------------------------------------------
c
      dummy = psp(i)**2-2.0d0*hp(1)*psp(i) - hp(2)
      dpsdhp(i,1) = psp(i)**2/dummy
      dpsdhp(i,2) = psp(i) / dummy
      dpsdhp(i,3) = 2.0d0/3.0d0/dummy
c
      return
      end subroutine ummdp_yld2004_18p_dpsdhp
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CALCULATE d2(psp)/d(hp)2
c
      subroutine ummdp_yld2004_18p_d2psdhp ( i,psp,hp,d2psdhp )
c-----------------------------------------------------------------------
      implicit none
c
      integer i
      real*8 psp(3),hp(3)
      real*8 d2psdhp(3,3,3)
c
      real*8 dummy
c-----------------------------------------------------------------------
      dummy = (psp(i)**2-2.0d0*hp(1)*psp(i)-hp(2)) ** 3
      d2psdhp(i,1,1) = 2.0d0*psp(i)**3*
     1                 (psp(i)**2-3.0d0*hp(1)*psp(i)-2.0d0*hp(2))/dummy
      d2psdhp(i,2,2) = -2.0d0*psp(i) * (hp(1)*psp(i)+hp(2)) / dummy
      d2psdhp(i,3,3) = -8.0d0/9.0d0 * (psp(i)-hp(1)) / dummy
      d2psdhp(i,1,2) = psp(i)**2 * 
     1                (psp(i)**2-4.0d0*hp(1)*psp(i)-3.0d0*hp(2)) / dummy
      d2psdhp(i,2,1) = d2psdhp(i,1,2)
      d2psdhp(i,2,3) = -2.0d0 * (psp(i)**2+hp(2)) / 3.0d0 / dummy
      d2psdhp(i,3,2) = d2psdhp(i,2,3)
      d2psdhp(i,3,1) = -4.0d0 * psp(i) * (hp(1)*psp(i)+hp(2)) /
     1                  3.0d0 / dummy
      d2psdhp(i,1,3) = d2psdhp(i,3,1)
c
      return
      end subroutine ummdp_yld2004_18p_d2psdhp
c
c
c
************************************************************************
c     YLD89 YIELD FUNCTION AND DERIVATIVES
c
c       doi: 
c
      subroutine ummdp_yld89 ( s,se,dseds,d2seds2,nreq,pryld,ndyld )
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
      call ummdp_yld89_branch ( s,se,dseds,d2seds2,nreq,
     1                           pryld,ndyld )
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
        call ummdp_yld89_branch ( s0,se0,dsedsz,d2seds2z,1,
     1                             pryld,ndyld )
        ta1 = dsedsz(1)
        ta2 = dsedsz(2)
        ta3 = dsedsz(3)
c
        s0(i) = s(i) - delta
        call ummdp_yld89_branch ( s0,se0,dsedsz,d2seds2z,1,
     2                             pryld,ndyld )
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
      end subroutine ummdp_yld89
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     BRANCH OF YLD89
c
      subroutine ummdp_yld89_branch ( s,se,dseds,d2seds2,nreq,
     1                                pryld,ndyld )
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
c                                                       ---- clear start
      call ummdp_utility_clear1 ( dfdK,2 )
      call ummdp_utility_clear2 ( dKds,2,3 )
c
      call ummdp_utility_clear2 ( d2fdK2,2,2 )
      call ummdp_utility_clear2 ( dfdKdfdK,2,2 )
      call ummdp_utility_clear2 ( d2K1ds2,3,3 )
      call ummdp_utility_clear2 ( d2K2ds2,3,3 )
c
      call ummdp_utility_clear1 ( df1dK,2 )
      call ummdp_utility_clear1 ( df2dK,2 )
      call ummdp_utility_clear1 ( df3dK,2 )
c
      call ummdp_utility_clear2 ( d2f1dK2,2,2 )
      call ummdp_utility_clear2 ( d2f2dK2,2,2 )
      call ummdp_utility_clear2 ( d2f3dK2,2,2 )
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
        call ummdp_utility_clear2 ( d2K1ds2,3,3 )
c                                                      ---- d2K2ds2(3,3)
        if ( pK2 == 0.0d0 ) then
          DpK22 = 1.0d-32
c
          d2K2ds2(1,1) = (DpK22**(-0.5d0)-pK3*DpK22**(-1.5d0)) / 4.0d0
          d2K2ds2(1,2) = -h * d2K2ds2(1,1)
          d2K2ds2(1,3) = -p*p*s(3)*pK3/2.0d0*DpK22**(-1.5d0)

          d2K2ds2(2,1) = d2K2ds2(1,2)
          d2K2ds2(2,2) = h * h * d2K2ds2(1,1)
          d2K2ds2(2,3) = -h * d2K2ds2(1,3)

          d2K2ds2(3,1) = d2K2ds2(1,3)
          d2K2ds2(3,2) = d2K2ds2(2,3)
          d2K2ds2(3,3) = p*p*(DpK22**(-0.5d0) - 
     1                   p*p*s(3)*s(3)*DpK22**(-1.5d0))
c
        else
c
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
c
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
        do i = 1,2
          dfdK(i) = df1dK(i)+df2dK(i)+df3dK(i)
        end do
c                                                             ---- dsedf
        dsedf = 0.0d0
        dsedf = (0.5d0*f) ** ((1.0d0-pM)/pM)
        dsedf = dsedf / (2.0d0*pM)
c                                                             ---- dseds
        do i = 1,3
          dseds(i) = 0.0d0
          do j = 1,2
            dseds(i) = dseds(i)+dfdK(j)*dKds(j,i)
          end do
          dseds(i) = dseds(i)*dsedf
        end do
      end if
c
c                                                           ---- d2seds2
      if ( nreq >= 2 ) then
c                                                            ---- d2fdK2
c
c                                                           ---- d2f1dK2
c
        do i = 1,2
          do j = 1,2
            d2f1dK2(i,j) = a*pM*(pM-1.0d0)*abs(pK1+pK2)**(pM-2.0d0)
          end do
        end do
c
c                                                           ---- d2f2dK2
c
        d2f2dK2(1,1) = a*pM*(pM-1.0d0)*abs(pK1-pK2)**(pM-2.0d0)
        d2f2dK2(1,2) = -a*pM*(pM-1.0d0)*abs(pK1-pK2)**(pM-2.0d0)
        d2f2dK2(2,1) = d2f2dK2(1,2)
        d2f2dK2(2,2) = d2f2dK2(1,1)
c
c                                                           ---- d2f3dK2
c
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
        d2sedf2 = 0.0d0
        d2sedf2 = (0.5d0*f)**((1.0d0-2.0d0*pM)/pM)
        d2sedf2 = d2sedf2*(1.0d0-pM) / (4.0d0*pM*pM)
c
        do i = 1,3
          do j = 1,3
c
            call ummdp_utility_clear1 (w,2)
            call ummdp_utility_clear1 (x,2)
c
            do k = 1,2
              do l = 1,2
                w(k) = w(k) + dfdKdfdK(k,l)*dKds(l,j)
                x(k) = x(k) + d2fdK2(k,l)*dKds(l,j)
              end do
            end do
c
            sc1 = 0.0d0
            sc2 = 0.0d0
            sc3 = 0.0d0
c
            do k = 1,2
              sc1 = sc1 + dKds(k,i)*w(k)
              sc2 = sc2 + dKds(k,i)*x(k)
            end do
            sc3 = dfdK(1)*d2K1ds2(j,i) + dfdK(2)*d2K2ds2(j,i)
c
            d2seds2(i,j) = d2sedf2*sc1 + dsedf*sc2 + d2sedf2*sc3
c
          end do
        end do
      end if
c
      return
      end subroutine ummdp_yld89_branch
c
c
c
************************************************************************
c     YOSHIDA2011 YIELD FUNCTION AND DERIVATIVES
c
c       doi: https://doi.org/10.1016/j.ijplas.2013.01.010
c
      subroutine ummdp_yoshida2011 ( s,se,dseds,d2seds2,nreq,
     1                               pryld,ndyld )
c-----------------------------------------------------------------------
      implicit none
c
      integer,parameter :: maxa = 100
c
      integer,intent(in) :: nreq,ndyld
      real*8 ,intent(in) :: s(6),pryld(ndyld)
c
      real*8,intent(out) :: se
      real*8,intent(out) :: dseds(6)
      real*8,intent(out) :: d2seds2(6,6)
c
      integer nd0,n,it,nterms,ndmax,jy,jx
      integer ipow(maxa,3)
      real*8 a(maxa)
c-----------------------------------------------------------------------
c
      nd0 = 3
c
      n = 0
      do it = 0,nd0
        n = n + (nd0-it)*2 + 1
      end do
      nterms = n
      if ( maxa < nterms ) then
        write (6,*) 'increase maxa :',maxa,nterms
        call ummdp_exit ( 9000 )
      end if
c
      n = 0
      ipow = 0
      do it = 0,nd0
        ndmax = nd0*2 - it*2
        do jy = 0,ndmax
          jx = ndmax - jy
          n = n + 1
          ipow(n,1) = jx
          ipow(n,2) = jy
          ipow(n,3) = it
        end do
      end do
c
      a = 0.0d0
      a(1) =  1.0d0 *          pryld(1+1)           !       c1
      a(2) = -3.0d0 *          pryld(1+2)           !    -3*c2
      a(3) =  6.0d0 *          pryld(1+3)           !     6*c3
      a(4) = -7.0d0 *          pryld(1+4)           !    -7*c4
      a(5) =  6.0d0 *          pryld(1+5)           !     6*c5
      a(6) =- 3.0d0 *          pryld(1+6)           !    -3*c6
      a(7) =  1.0d0 *          pryld(1+7)           !       c7
      a(8) =  1.0d0 * 9.0d0  * pryld(1+8)           !  9   *c8
      a(9) = -2.0d0 * 9.0d0  * pryld(1+9)           !  9*-2*c9
      a(10) = 3.0d0 * 9.0d0  * pryld(1+10)          !  9* 3*c10
      a(11) =-2.0d0 * 9.0d0  * pryld(1+11)          !  9*-2*c11
      a(12) = 1.0d0 * 9.0d0  * pryld(1+12)          !  9   *c12
      a(13) = 1.0d0 * 27.0d0 * pryld(1+13)          ! 27   *c13
      a(14) =-1.0d0 * 27.0d0 * pryld(1+14)          ! 27*-1*c14
      a(15) = 1.0d0 * 27.0d0 * pryld(1+15)          ! 27   *c15
      a(16) = 1.0d0 * 27.0d0 * pryld(1+16)          ! 27   *c16
c
      call ummdp_hy_polytype ( s,se,dseds,d2seds2,nreq,nd0,
     1                          a,ipow,maxa,nterms )
c
      return
      end subroutine ummdp_yoshida2011
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     HU2005 & YOSHIDA2011 STYLE POLYNOMIAL TYPE YIELD FUNCTION
c
      subroutine ummdp_hy_polytype ( s,se,dseds,d2seds2,nreq,nd0,
     1                                a,ipow,maxa,nterms )
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: nreq,nd0,maxa,nterms
      integer,intent(in) :: ipow(maxa,3)
      real*8 ,intent(in) :: s(6),a(maxa)
c
      real*8,intent(out) :: se
      real*8,intent(out) :: dseds(6)
      real*8,intent(out) :: d2seds2(6,6)
c    
      integer i,j,k,n,id,idmax,jdmax,jd,nd
      integer ii(3)
      real*8 dinv,fai,q,dd,ff,fff,ddi,ddj
      real*8 sterm(3),v(6)
c-----------------------------------------------------------------------
c     nd        : order of polynominal nd=2*nd0
c     a(n)      : constants of function
c     ipow(n,i) : power of terms
c-----------------------------------------------------------------------
c
      nd = nd0 * 2
      dinv = 1.0d0 / float(nd)
c
      sterm(1) = s(1) - s(3)                  ! sx-sz
      sterm(2) = s(2) - s(3)                  ! sy-sz
      sterm(3) = s(4)**2 + s(5)**2 + s(6)**2  ! txy^2+tyz^2+tzx^2
c
      fai = 0.0d0
      do n = 1,nterms
        q = a(n)
        do k = 1,3
          if ( ipow(n,k) > 0 ) then
            q = q * sterm(k)**ipow(n,k)
          end if
        end do
        fai = fai + q
      end do
      se = fai ** dinv
      if ( nreq == 0 ) return
c
      v = 0.0d0
      do i = 1,6
        idmax = 1
        if ( i == 3 ) idmax = 2
        do id = 1,idmax
          do n = 1,nterms
            do k = 1,3
              ii(k) = ipow(n,k)
            end do
            select case ( i )
            case ( 1,2 )
              dd = float(ii(i))
              ii(i) = ii(i) - 1
            case ( 3 )
              dd = -1.0d0 * float(ii(id))
              ii(id) = ii(id) - 1
            case default
              dd = 2.0d0 * s(i) * float(ii(3))
              ii(3) = ii(3) - 1
            end select
            q = dd * a(n)
            do k = 1,3
              if ( ii(k) > 0 ) then
                q = q * sterm(k)**ii(k)
              else if ( ii(k) < 0 ) then
                q = 0.0
              end if
            end do
            v(i) = v(i) + q
          end do
        end do
      end do
      ff = dinv * fai**(dinv-1.0d0)
      do i = 1,6
        dseds(i) = ff * v(i)
      end do
      if ( nreq == 1 ) return
c
      fff = dinv * (dinv-1.0d0) * fai**(dinv-2.0d0)
      do i = 1,6
        do j = 1,6
          d2seds2(i,j) = fff * v(i) * v(j)
        end do
      end do
      do i = 1,6
        idmax = 1
        if ( i == 3 ) idmax = 2
        do id = 1,idmax
          do j = 1,6
            jdmax = 1
            if (  j == 3 ) jdmax = 2
            if ( ( j > 3 ) .and. ( i == j ) ) jdmax = 2
            do jd = 1,jdmax
              do n = 1,nterms
                do k = 1,3
                  ii(k) = ipow(n,k)
                end do
                select case ( i )
                case ( 1,2 )
                  ddi = float(ii(i))
                  ii(i) = ii(i) - 1
                case ( 3 )
                  ddi = -1.0d0 * float(ii(id))
                  ii(id) = ii(id) - 1
                case default
                  ddi = 2.0d0 * s(i) * float(ii(3))
                  ii(3) = ii(3) - 1
                end select 
                select case ( j )
                case ( 1,2 )
                  ddj = float(ii(j))
                  ii(j) = ii(j) - 1
                case ( 3 )
                  ddj = -1.0d0 * float(ii(jd))
                  ii(jd) = ii(jd) - 1
                case default
                  if ( jd == 1 ) then
                    ddj = 2.0d0 * s(j) * float(ii(3))
                    ii(3) = ii(3) - 1
                  else
                    ddi = 2.0d0 * float(ipow(n,3))
                    ddj = 1.0d0
                  end if
                end select
                q = a(n) * ddi * ddj
                do k = 1,3
                  if ( ii(k) > 0 ) then
                    q = q * sterm(k)**ii(k)
                  else if ( ii(k) < 0 ) then
                    q = 0.0d0
                  end if
                end do
                d2seds2(i,j) = d2seds2(i,j) + ff*q
              end do
            end do
          end do
        end do
      end do
c
      return
      end subroutine ummdp_hy_polytype
c
c
c
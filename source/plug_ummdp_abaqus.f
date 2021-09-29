c***********************************************************************
c
c     UMMDp - Unified Material Model Driver for Plasticity
c
c***********************************************************************
c
c     > Copyright (c) 2018 JANCAE
c       . This software includes code originally developed by the  
c       Material Modeling Working group of JANCAE.
c
c     > Extended by M.G. Oliveira from University of Aveiro, Portugal
c       . Added additional isotropic hardening laws
c       . Corrected order of Voigt notation for Yld2004-18p with Abaqus
c       . Linked kinematic hardening laws to the core of UMMDp
c       . Added Chaboche kinematic hardening law as used by Abaqus
c     	. Implemented uncouple rupture criteria
c
c***********************************************************************
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
      call jancae_debugmode ( nvbs,nvbs0 )
c                                       ---- output detailed information
      if ( nvbs >= 4 ) then
        call jancae_printinfo  ( kinc,ndi,nshr )
        call jancae_printinout ( 0,stress,dstran,ddsdde,ntens,
     &                           statev,nstatv )
      end if
c
c                                           ---- set material properties
      do i = 2,nprops
        prop(i-1) = props(i)
      end do
c
      call jancae_prop_dim ( prop,nprop,propdim,
     &                       ndela,ndyld,ndihd,ndkin,
     &                       npbs,ndrup )
      if ( npbs > mxpbs ) then
        write (6,*) 'npbs > mxpbs error in umat'
        write (6,*) 'npbs =',npbs
        write (6,*) 'mxpbs=',mxpbs
        call jancae_exit ( 9000 )
      end if
c                                                      ---- check nstatv
      call jancae_check_nisv ( nstatv,ntens,npbs )
c                             ---- copy current internal state variables
      call jancae_isvprof ( isvrsvd,isvsclr )
      call jancae_isv2pex ( isvrsvd,isvsclr,
     &                      statev,nstatv,
     &                      p,pe,x1,ntens,mxpbs,npbs )
c
c                             ---- update stress and set tangent modulus
      mjac = 1
      call jancae_plasticity ( stress,s2,dstran,
     &                         p,dp,dpe,de33,
     &                         x1,x2,mxpbs,
     &                         ddsdde,
     &                         ndi,nshr,ntens,
     &                         nvbs,mjac,
     &                         prop,nprop,propdim )
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
        call jancae_printinout ( 1,stress,dstran,ddsdde,ntens,
     &                           statev,nstatv )
      end if
c
      return
      end
c
c
c
c***********************************************************************
c
      SUBROUTINE SDVINI ( STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,
     &                    LAYER,KSPT )
c
      INCLUDE 'ABA_PARAM.INC'
c
      DIMENSION STATEV(NSTATV),COORDS(NCRDS)
c
c***********************************************************************
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
      end
c
c
c
c***********************************************************************
c
      SUBROUTINE UVARM ( UVAR,DIRECT,T,TIME,DTIME,CMNAME,ORNAME,
     &    NUVARM,NOEL,NPT,LAYER,KSPT,KSTEP,KINC,NDI,NSHR,COORD,
     &    jmac,jmatyp,matlayo,laccfla )
c
      INCLUDE 'ABA_PARAM.INC'
c
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  flgray(15)
      DIMENSION UVAR(NUVARM),DIRECT(3,3),T(3,3),TIME(2)
      DIMENSION array(15),jarray(15),jmac(*),jmatyp(*),COORD(*)
c
c***********************************************************************
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
      call getvrm ( 'S',array,jarray,flgray,jrcd,jmac,jmatyp,
     &                  matlayo,laccfla )
      if ( jrcd /= 0 ) then
        write (6,*) 'request error in uvarm for s'
        write (6,*) 'stop in uvrm.'
        call jancae_exit ( 9000 )
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
        call jancae_exit ( 9000 )
      end if
      call getvrm ( 'SDV',ARRAY2,JARRAY2,FLGRAY2,jrcd,jmac,jmatyp,
     &                    matlayo,laccfla)
      if ( jrcd /= 0 ) then
        write (6,*) 'request error in uvarm for sdv'
        write (6,*) 'stop in uvrm.'
        call jancae_exit ( 9000 )
      end if
      do i = 1,nsdv
        sdv(i) = array2(i)
      end do
c                                           ---- set material properties
      call jancae_prop_dim ( prop,nprop,propdim,
     &                       ndela,ndyld,ndihd,ndkin,
     &                       npbs,ndrup )
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
      call jancae_isvprof ( isvrsvd,isvsclr )
      call jancae_isv2pex ( isvrsvd,isvsclr,
     &                      sdv,maxsdv,
     &                      p,pe,x,ntens,mxpbs,npbs )
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
        call jancae_yfunc ( se,dseds,d2seds2,0,
     &                      eta,ntens,ndi,nshr,
     &                      pryld,ndyld )
        uvar(1) = se
      end if
c                                                       ---- flow stress
      if ( nuvarm >= 2 ) then
        call jancae_hardencurve ( sy,dsydp,d2sydp2,0,p,prihd,ndihd )
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
          call jancae_rupture ( ntens,sdv,nsdv,uvar,uvar1,nuvarm,
     &                          jrcd,jmac,jmatyp,matlayo,laccfla,
     &                          nt,ndrup,prrup)
        end if
      end if
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
c
c-----------------------------------------------------------------------
      INCLUDE 'ABA_PARAM.INC'
c
c
      isvrsvd = 0             ! no reserved variables
c
      isvsclr = 1             ! statev(1) is for eq.plast.strain
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
c
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
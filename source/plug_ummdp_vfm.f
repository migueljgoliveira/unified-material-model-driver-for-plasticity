c***********************************************************************
c
c     UMMDp - Unified Material Model Driver for Plasticity
c
c     VFM Edition
c
c***********************************************************************
c
c     > Copyright (c) 2018 JANCAE
c       . This software includes code originally developed by the  
c       Material Modeling Working group of JANCAE.
c
c     > Developed by M.G. Oliveira from University of Aveiro, Portugal
c
c***********************************************************************
c
      SUBROUTINE UMAT ( STRESS,STATEV,STRAN,DSTRAN,DROT,
     &                  NDI,NSHR,NTENS,NSTATV,
     &                  PROPS,NPROPS,
     &                  NOEL,NPT,KSPT,KINC )
c
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     &          STRAN(NTENS),DSTRAN(NTENS),
     &          PROPS(NPROPS),DROT(3,3)
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

      dimension ddsdde(ntens,ntens)
c
      dimension ustatev(6)
c
      parameter (mxprop=100,nrot=3)
      dimension prop(mxprop)
c-----------------------------------------------------------------------
c
cf2py intent(in) stress,statev,stran,dstran,drot
cf2py intent(in) ndi,nshr,ntens,nstatv
cf2py intent(in) props,nprops
cf2py intent(in) noel,npt,kspt,kinc
cf2py intent(out) stress,statev
cf2py depend(ntens) stress,stran,dstran
cf2py depend(nstatv) statev
cf2py depend(nrot,nrot) drot
cf2py depend(nprops) props
c
c-----------------------------------------------------------------------
c
c                                                   ---- open debug file
      ! open(unit=6,file=trim('ummdp.log'),status='NEW')
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
      ! call jancae_debugmode ( nvbs,nvbs0 )
! c                                       ---- output detailed information
!       if ( nvbs >= 4 ) then
!         call jancae_printinfo  ( kinc,ndi,nshr )
!         call jancae_printinout ( 0,stress,dstran,ddsdde,ntens,
!      &                           statev,nstatv )
!       end if
! c
! c                                           ---- set material properties
!       do i = 2,nprops
!         prop(i-1) = props(i)
!       end do
! c		
!       call jancae_prop_dim ( prop,nprop,propdim,
!      &                       ndela,ndyld,ndihd,ndkin,
!      &                       npbs,ndrup )
!       if ( npbs > mxpbs ) then
!         write (6,*) 'npbs > mxpbs error in umat'
!         write (6,*) 'npbs =',npbs
!         write (6,*) 'mxpbs=',mxpbs
!         call jancae_exit ( 9000 )
!       end if
! c                                                      ---- check nstatv
!       call jancae_check_nisv ( nstatv,ntens,npbs )
! c                             ---- copy current internal state variables
      ! call jancae_isvprof ( isvrsvd,isvsclr )
!       call jancae_isv2pex ( isvrsvd,isvsclr,
!      &                      statev,nstatv,
!      &                      p,pe,x1,ntens,mxpbs,npbs )
! c
! c                             ---- update stress and set tangent modulus
!       mjac = 0
!       call jancae_plasticity ( stress,s2,dstran,
!      &                         p,dp,dpe,de33,
!      &                         x1,x2,mxpbs,
!      &                         ddsdde,
!      &                         ndi,nshr,ntens,
!      &                         nvbs,mjac,
!      &                         prop,nprop,propdim )
! c                                                     ---- update stress
!       do i = 1,ntens
!         stress(i) = s2(i)
!       end do
! c                                  ---- update equivalent plastic strain
!       statev(isvrsvd+1) = p + dp
! c                                  ---- update plastic strain components
!       call rotsig ( statev(isvrsvd+2),drot,ustatev,2,ndi,nshr )
! c
!       do i = 1,ntens
!         is = isvrsvd + isvsclr + i
!         statev(is) = ustatev(i) + dpe(i)
!       end do
! c                                  ---- update of back stress components
!       if ( npbs /= 0 ) then
!         do n = 1,npbs
!           do i = 1,ntens
!             is = isvrsvd + isvsclr + ntens*n + i
!             statev(is) = x2(n,i)
!           end do
!         end do
!       end if
! c                           ----  if debug mode, output return arguments
!       if ( nvbs >= 4 ) then
!         call jancae_printinout ( 1,stress,dstran,ddsdde,ntens,
!      &                           statev,nstatv )
!       end if
c                                                  ---- close debug file
      ! close(6)
c
      return
      end
c
c
c
! c-----------------------------------------------------------------------
! c     set internal state variables profile
! c
!       subroutine jancae_isvprof ( isvrsvd,isvsclr )
! c
! c-----------------------------------------------------------------------
! c
!       isvrsvd = 0           ! no reserved variables
! c
!       isvsclr = 1           ! statev(1) is for equivalent plastic strain
! c
!       return
!       end
! c
! c
! c
! c-----------------------------------------------------------------------
! c     rotate a tensor
! c
!       subroutine rotsig ( statev,drot,ustatev,lstr,ndi,nshr )
! c
! c-----------------------------------------------------------------------
!       dimension statev(6),ustatev(6),drot(3,3)
!       dimension aux1(ndi,ndi),aux2(ndi,ndi),aux3(ndi,ndi)
!       dimension auxrot(ndi,ndi)
! c                                              ---- set statev to tensor
!       call jancae_clear2( aux1,ndi,ndi )
!       do i = 1,ndi
!         do j = 1,ndi
!           if ( i == j ) then
!             aux1(i,j) = statev(i)
!           else
!             if ( lstr == 1 ) then
!               aux1(i,j) = statev(i+j+1)
!             else if ( lstr == 2 ) then
!               aux1(i,j) = statev(i+j+1)/2.0d0
!             end if
!           end if
!         end do
!       end do
! c                                               ---- copy drot to auxrot
!       call jancae_clear2( aux1,ndi,ndi )
!       do i = 1,ndi
!         do j = 1,ndi
!           auxrot(i,j) = drot(i,j)
!         end do
!       end do
! c                                                     ---- rotate statev
!       call jancae_mm ( aux2,auxrot,aux1,ndi,ndi,ndi )
!       call jancae_mm ( aux3,aux2,transpose(auxrot),ndi,ndi,ndi )
! c                                     ---- set rotated tensor to ustatev
!       do i = 1,ndi
!         do j = 1,ndi
!           if ( i == j ) then
!             ustatev(i) = aux3(i,j)
!           else
!             if ( lstr == 1 ) then
!               ustatev(i+j) = aux3(i,j)
!             else if ( lstr == 2 ) then
!               ustatev(i+j) = aux3(i,j)*2.0d0
!             end if
!           end if
!         end do
!       end do
! c
!       return
!       end
! c
! c
! c
! c-----------------------------------------------------------------------
! c     exit program by error
! c
!       subroutine jancae_exit (nexit)
! c
! c-----------------------------------------------------------------------
!       common /jancae1/ne,ip,lay
! c                                  nexit : exit code
!       write (6,*) 'error code :',nexit
!       write (6,*) 'element no.           :',ne
!       write (6,*) 'integration point no. :',ip
!       write (6,*) 'layer no.             :',lay
! c
!       ! stop
! c
!       return
!       end
c
c
c
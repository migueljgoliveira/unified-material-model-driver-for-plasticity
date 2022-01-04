************************************************************************
*                                                                      *
*                                  UMMDp                               *
*                                                                      *
*                             <><><><><><><>                           *
*                                                                      *
*              UNIFIED MATERIAL MODEL DRIVER FOR PLASTICITY            *
*                                                                      *
*                  < PLUG-IN FOR VIRTUAL FIELDS METHOD >               *
*                                                                      *
************************************************************************
*                                                                      *
*     > Copyright (c) 2018 JANCAE                                      *
*       . This software includes code originally developed by the      *
*       Material Modeling Working group of JANCAE.                     *
*                                                                      *
*     > Extended by M.G. Oliveira from University of Aveiro, Portugal  *
*       . Added additional isotropic hardening laws                    *
*       . Corrected Voigt notation for Yld2004-18p with Abaqus         *
*       . Linked kinematic hardening laws to the core of UMMDp         *
*       . Added Chaboche kinematic hardening law as used by Abaqus     *
*     	. Implemented uncoupled rupture criteria                       *
*       . Modified code to explicit declaration                        *
*                                                                      *
************************************************************************
c
c     UMMDP-VFM MAIN SUBROUTINE
c
      subroutine ummdp_vfm ( stress1,statev1,strain,dstrain,ndi,nshr,
     1                       ntens,nstatev,props,nprops,noel,npt,kinc,
     2                       stress2,statev2,nexit )
c
c-----------------------------------------------------------------------
      implicit none
c
      common /ummdp1/ne,ip,lay
      common /ummdp2/prop
      common /ummdp3/nsdv
      common /ummdp4/propdim
      common /ummdp5/nexito
c
      integer,intent(in) :: ndi,nshr,ntens,nstatev,nprops,noel,npt,kinc                      
      real*8 ,intent(in) :: props(nprops)
      real*8 ,intent(in) :: stress1(ntens),statev1(nstatev),
     1                      strain(ntens),dstrain(ntens)
c
      integer,intent(out) :: nexit
      real*8 ,intent(out) :: stress2(ntens),statev2(nstatev)
c 
      integer mxpbs,mxprop,nrot
      parameter (mxpbs=10,mxprop=100)
      integer ne,ip,lay,nsdv,propdim,nexito,i,k,n,is,nprop,nvbs0,nvbs,
     1        ndela,ndyld,ndihd,ndkin,npbs,ndrup,mjac,isvrsvd,isvsclr,
     2        maxsdv
      real*8 de33,p,dp
      real*8 ustatev(ntens),s2(ntens),dpe(ntens),pe(ntens),prop(mxprop)
      real*8 x1(mxpbs,ntens),x2(mxpbs,ntens),ddsdde(ntens,ntens)
      character*100 text
c-----------------------------------------------------------------------
c
cf2py intent(in) stress1,statev1,strain,dstrain
cf2py intent(in) ndi,nshr,ntens,nstatev
cf2py intent(in) props,nprops
cf2py intent(in) noel,npt,kspt,kinc
cf2py intent(out) stress2,statev2,nexit
cf2py depend(ntens) stress1,stress2,strain,dstrain
cf2py depend(nstatev) statev1,statev2
cf2py depend(nprops) props
c
c-----------------------------------------------------------------------
c
c                                                   ---- open debug file
      if ( kinc == 1 ) then
        open(6,file='ummdp_vfm.log')
      else
        open(6,file='ummdp_vfm.log',access='APPEND',status='OLD')
      end if
c 
      ne = noel
      ip = npt
      lay = 1
      nsdv = nstatev
      nprop = mxprop
      propdim = nprops - 1
c                                        ---- set debug and verbose mode
      nvbs0 = props(1)
      call ummdp_debugmode ( nvbs,nvbs0 )
c                                       ---- print detailed information
      if ( nvbs >= 1 ) then
        call ummdp_print_info  ( kinc,ndi,nshr )
      end if
c                                             ---- print input arguments
      if ( nvbs >= 4 ) then
        ddsdde = 0.0d0
        call ummdp_print_inout ( 0,stress1,dstrain,ddsdde,ntens,statev1,
     1                           nstatev )
      end if
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
        call ummdp_exit ( 301 )
      end if
c                                                     ---- check nstatev
      call ummdp_check_nisv ( nstatev,ntens,npbs )
c                             ---- copy current internal state variables
      call ummdp_isvprof ( isvrsvd,isvsclr )
      call ummdp_isv2pex ( isvrsvd,isvsclr,statev1,nstatev,p,pe,x1,
     1                     ntens,mxpbs,npbs )
c
c                             ---- update stress and set tangent modulus
      mjac = 0
      call ummdp_plasticity ( stress1,s2,dstrain,p,dp,dpe,de33,x1,x2,
     1                        mxpbs,ddsdde,ndi,nshr,ntens,nvbs,mjac,
     2                        prop,nprop,propdim )
c                                                     ---- update stress
      do i = 1,ntens
        stress2(i) = s2(i)
      end do
c                                  ---- update equivalent plastic strain
      statev2(isvrsvd+1) = p + dp
c                                  ---- update plastic strain components
      do i = 1,ntens
        is = isvrsvd + isvsclr + i
        statev2(is) = statev1(is) + dpe(i)
      end do
c                                     ---- update back stress components
      if ( npbs /= 0 ) then
        do n = 1,npbs
          do i = 1,ntens
            is = isvrsvd + isvsclr + ntens*n + i
            statev2(is) = x2(n,i)
          end do
        end do
      end if
c                                            ---- print output arguments
      if ( nvbs >= 4 ) then
        call ummdp_print_inout ( 1,stress2,dstrain,ddsdde,ntens,
     1                           statev2,nstatev )
      end if
c                                                  ---- close debug file
      close(6)
c                                                 ---- return error code
      nexit = nexito
c
      return
      end
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     SET INTERNAL STATE VARIABLES PROFILE
c
      subroutine ummdp_isvprof ( isvrsvd,isvsclr )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(out) :: isvrsvd,isvsclr
c-----------------------------------------------------------------------

      isvrsvd = 0           ! no reserved variables
c
      isvsclr = 1           ! statev(1) is for equivalent plastic strain
c
      return
      end subroutine ummdp_isvprof
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     EXIT PROGRAM BY ERROR
c
      subroutine ummdp_exit ( nexit )
c
c-----------------------------------------------------------------------
      implicit none
c
      common /ummdp1/ne,ip,lay
      common /ummdp5/nexito
c
      integer,intent(in) :: nexit
c
      integer ne,ip,lay,nexito
      character*50 fmt1,fmt2,tmp
c-----------------------------------------------------------------------
c
      write(6,'(4/8xA)') '!!!!!!!!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!!!!!!'
c
      write (tmp,'(I)') nexit
      write (6,'(/12xA,A)') '             Code : ',adjustl(tmp)
c
      call ummdp_print_element ( )
c
      write(6, '(/8xA)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
c
c
      nexito = -1
c
      return
      end subroutine ummdp_exit
c
c
c
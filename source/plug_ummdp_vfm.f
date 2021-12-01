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
*       . Modified code to use only explicit variables                 *
*                                                                      *
************************************************************************
c
c     UMMDP-VFM MAIN SUBROUTINE
c
      subroutine ummdp_vfm ( stress,statev,strain,dstrain,ndi,nshr,
     1                       ntens,nstatev,props,nprops,noel,npt,kinc,
     2                       nexit )
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
      real*8 ,intent(in) :: strain(ntens),dstrain(ntens)
c
      real*8 ,intent(inout) :: stress(ntens),statev(nstatev)
c
      integer,intent(out) :: nexit
c 
      integer mxpbs,mxprop,nrot
      parameter (mxpbs=10,mxprop=100)
      integer ne,ip,lay,nsdv,propdim,nexito,i,k,n,is,nprop,nvbs0,nvbs,
     1        ndela,ndyld,ndihd,ndkin,npbs,ndrup,mjac,isvrsvd,isvsclr,
     2        maxsdv
      real*8 de33,p,dp
      real*8 ustatev(ntens),s2(ntens),dpe(ntens),pe(ntens),prop(mxprop)
      real*8 x1(mxpbs,ntens),x2(mxpbs,ntens),ddsdde(ntens,ntens)
c-----------------------------------------------------------------------
c
cf2py intent(in,out) stress,statev
cf2py intent(in) strain,dstrain
cf2py intent(in) ndi,nshr,ntens,nstatev
cf2py intent(in) props,nprops
cf2py intent(in) noel,npt,kspt,kinc
cf2py intent(out) nexit
cf2py depend(ntens) stress,strain,dstrain
cf2py depend(nstatev) statev
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
        call ummdp_print_inout ( 0,stress,dstrain,ddsdde,ntens,statev,
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
      call ummdp_isv2pex ( isvrsvd,isvsclr,statev,nstatev,p,pe,x1,ntens,
     1                     mxpbs,npbs )
c
c                             ---- update stress and set tangent modulus
      mjac = 1
      call ummdp_plasticity ( stress,s2,dstrain,p,dp,dpe,de33,x1,x2,
     1                        mxpbs,ddsdde,ndi,nshr,ntens,nvbs,mjac,
     2                        prop,nprop,propdim )
c                                                     ---- update stress
      do i = 1,ntens
        stress(i) = s2(i)
      end do
c                                  ---- update equivalent plastic strain
      statev(isvrsvd+1) = p + dp
c                                  ---- update plastic strain components
      do i = 1,ntens
        is = isvrsvd + isvsclr + i
        statev(is) = statev(is) + dpe(i)
      end do
c                                     ---- update back stress components
      if ( npbs /= 0 ) then
        do n = 1,npbs
          do i = 1,ntens
            is = isvrsvd + isvsclr + ntens*n + i
            statev(is) = x2(n,i)
          end do
        end do
      end if
c                                            ---- print output arguments
      if ( nvbs >= 4 ) then
        call ummdp_print_inout ( 1,stress,dstrain,ddsdde,ntens,statev,
     1                           nstatev )
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
c-----------------------------------------------------------------------
c
      write (6,*) 'error code :',nexit
      write (6,*) 'element no.           :',ne
      write (6,*) 'integration point no. :',ip
      write (6,*) 'layer no.             :',lay           
c
      nexito = -1
c
      return
      end subroutine ummdp_exit
c
c
c
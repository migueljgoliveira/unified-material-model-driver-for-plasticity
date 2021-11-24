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
c
      subroutine ummdp_vfm ( stress,statev,strain,dstrain,drot,ndi,nshr,
     1                       ntens,nstatev,props,nprops,noel,npt,kspt,
     2                       kinc )
c        
c-----------------------------------------------------------------------
      implicit none
c
      common /jancae1/ne,ip,lay
      common /jancae3/prop
      common /jancaea/nsdv
      common /jancaeb/propdim
      integer ne,ip,lay,nsdv,propdim
c
      integer,intent(in) :: ndi,nshr,ntens,nstatev,nprops,noel,npt,kspt,
     1                      kinc
      real*8 ,intent(in) :: props(nprops)
      real*8 ,intent(in) :: strain(ntens),dstrain(ntens)      
      real*8 ,intent(in) :: drot(3,3)
c
      real*8 ,intent(inout) :: stress(ntens),statev(nstatev)
c
      integer mxpbs,mxprop,nrot
      parameter (mxpbs=10,mxprop=100,nrot=3)
      integer i,k,n,is,nprop,nvbs0,nvbs,ndela,ndyld,ndihd,ndkin,npbs,
     1        ndrup,mjac,isvrsvd,isvsclr,maxsdv
      real*8 de33,p,dp
      real*8 ustatev(6),s2(ntens),dpe(ntens),pe(ntens),prop(mxprop)
      real*8 x1(mxpbs,ntens),x2(mxpbs,ntens),ddsdde(ntens,ntens)     
c-----------------------------------------------------------------------
c
cf2py intent(in) stress,statev,strain,dstrain,drot
cf2py intent(in) ndi,nshr,ntens,nstatev
cf2py intent(in) props,nprops
cf2py intent(in) noel,npt,kspt,kinc
cf2py intent(out) stress,statev
cf2py depend(ntens) stress,strain,dstrain
cf2py depend(nstatev) statev
cf2py depend(nrot,nrot) drot
cf2py depend(nprops) props
c
c-----------------------------------------------------------------------
c
c                                                   ---- open debug file
      if ( kinc == 1 ) then
        open(6,file=trim('ummdp.log'),status='NEW')
      else
        open(6,file=trim('ummdp.log'),access='APPEND',status='OLD')
      end if
c
      ne = noel
      ip = npt
      lay = kspt
      if ( lay == 0 ) lay = 1
      nsdv = nstatev
      nprop = mxprop
      propdim = nprops - 1
c                                        ---- set debug and verbose mode
      nvbs0 = props(1)
      call ummdp_debugmode ( nvbs,nvbs0 )
c                                       ---- output detailed information
      if ( nvbs >= 4 ) then
        call ummdp_print_info  ( kinc,ndi,nshr )
        call ummdp_print_inout ( 0,stress,dstrain,ddsdde,ntens,statev,
     1                           nstatev )
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
c                                                     ---- check nstatev
      call ummdp_check_nisv ( nstatev,ntens,npbs )
c                             ---- copy current internal state variables
      call ummdp_isvprof ( isvrsvd,isvsclr )
      call ummdp_isv2pex ( isvrsvd,isvsclr,statev,nstatev,p,pe,x1,ntens,
     1                     mxpbs,npbs )                      
c
c                             ---- update stress and set tangent modulus
      mjac = 0
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
      call ummdp_rotsig ( statev(isvrsvd+2),drot,ustatev,2,ndi,nshr )
c
      do i = 1,ntens
        is = isvrsvd + isvsclr + i
        statev(is) = ustatev(i) + dpe(i)
      end do
c                                  ---- update of back stress components
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
        call ummdp_print_inout ( 1,stress,dstrain,ddsdde,ntens,statev,
     1                           nstatev )
      end if
c                                                  ---- close debug file
      close(6)
c
      return
      end
c
c
c
************************************************************************
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
************************************************************************
c     ROTATE A TENSOR
c
      subroutine ummdp_rotsig ( statev,drot,ustatev,lstr,ndi,nshr )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: lstr,ndi,nshr
      real*8 ,intent(in) :: statev(6)
      real*8 ,intent(in) :: drot(3,3)
c
      real*8,intent(out) :: ustatev(6)
c
      integer i,j
      real*8 aux1(ndi,ndi),aux2(ndi,ndi),aux3(ndi,ndi),auxrot(ndi,ndi)
c-----------------------------------------------------------------------
c
c                                              ---- set statev to tensor
      call ummdp_utility_clear2( aux1,ndi,ndi )
      do i = 1,ndi
        do j = 1,ndi
          if ( i == j ) then
            aux1(i,j) = statev(i)
          else
            if ( lstr == 1 ) then
              aux1(i,j) = statev(i+j+1)
            else if ( lstr == 2 ) then
              aux1(i,j) = statev(i+j+1)/2.0d0
            end if
          end if
        end do
      end do
c                                               ---- copy drot to auxrot
      call ummdp_utility_clear2( aux1,ndi,ndi )
      do i = 1,ndi
        do j = 1,ndi
          auxrot(i,j) = drot(i,j)
        end do
      end do
c                                                     ---- rotate statev
      call ummdp_utility_mm ( aux2,auxrot,aux1,ndi,ndi,ndi )
      call ummdp_utility_mm ( aux3,aux2,transpose(auxrot),ndi,ndi,ndi )
c                                     ---- set rotated tensor to ustatev
      do i = 1,ndi
        do j = 1,ndi
          if ( i == j ) then
            ustatev(i) = aux3(i,j)
          else
            if ( lstr == 1 ) then
              ustatev(i+j) = aux3(i,j)
            else if ( lstr == 2 ) then
              ustatev(i+j) = aux3(i,j)*2.0d0
            end if
          end if
        end do
      end do
c
      return
      end subroutine ummdp_rotsig
c
c
c
************************************************************************
c     EXIT PROGRAM BY ERROR
c
      subroutine ummdp_exit ( nexit )
c
c-----------------------------------------------------------------------
      implicit none
c
      common /jancae1/ne,ip,lay
      integer ne,ip,lay
c
      integer,intent(in) :: nexit
c-----------------------------------------------------------------------
c
      write (6,*) 'error code :',nexit
      write (6,*) 'element no.           :',ne
      write (6,*) 'integration point no. :',ip
      write (6,*) 'layer no.             :',lay
c
      ! stop
c
      return
      end subroutine ummdp_exit
c
c
c
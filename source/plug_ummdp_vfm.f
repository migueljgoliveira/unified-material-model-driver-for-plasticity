
      subroutine ummdp_vfm ( strain,ne,ndi,nshr,ntens,nstatev,props,
     1                       nprops,nf,fout,stress,statev,de33 )
c
c-----------------------------------------------------------------------
      implicit none
c
      common /ummdpa/n1234
c
      integer      ,intent(in) :: ne,ndi,nshr,ntens,nstatev,nprops,nf
      real*8       ,intent(in) :: props(nprops)
      real*8       ,intent(in) :: strain(nf,ne,ntens)
      character*100,intent(in) :: fout
c
      real*8 ,intent(out) :: de33(nf,ne)
      real*8 ,intent(out) :: stress(nf,ne,ntens),statev(nf,ne,nstatev)
c
      integer n1234,i,j,kinc,noel
      real*8 dstrain(nf,ne,ntens)
c-----------------------------------------------------------------------
c
cf2py intent(in) strain,ne,ndi,nshr,ntens,nstatev,props,nprops,nf,fout
cf2py intent(out) stress,statev,de33
cf2py depend(nprops) props
cf2py depend(nf,ne) de33
cf2py depend(nf,ne,ntens) stress,strain
cf2py depend(nf,ne,nstatev) statev
c
c-----------------------------------------------------------------------
c
c                                                   ---- open debug file
      if ( nint(props(1)) /= -1 ) then
        open(6,file='output/'//trim(fout)//'/'//trim(fout)//'.ummdp')
        n1234 = 0
      else
        n1234 = 1234
      end if
c
c                                      ---- stress integration variables
      stress = 0.0d0
      statev = 0.0d0
      de33 = 0.0d0
c                                         ---- loop over time increments
      do i = 2,nf
        kinc = i - 1
c                                                ---- loop over elements
c$omp parallel do
        do j = 1,ne
          noel = j
c                                            ---- total strain increment
          dstrain(i,j,:) = strain(i,j,:) - strain(i-1,j,:)
c
c                  ---- stress integration in corotational material csys
          call plug_ummdp_vfm ( stress(i-1,j,:),statev(i-1,j,:),
     1                          dstrain(i,j,:),ndi,nshr,ntens,nstatev,
     2                          props,nprops,noel,kinc,stress(i,j,:),
     3                          statev(i,j,:),de33(i,j) )
c
c                               ---- total strain in thickness direction
          de33(i,j) = de33(i,j) + de33(i-1,j)
c
        end do
c$omp end parallel do
      end do
c
c                                                  ---- close debug file
      if ( nint(props(1)) /= -1 ) then
        close(6)
      end if
c
      return
      end subroutine ummdp_vfm

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
c     UMMDp-VFM Plug-In
c
      subroutine plug_ummdp_vfm ( stress1,statev1,dstrain,ndi,nshr,
     1                            ntens,nstatev,props,nprops,noel,kinc,
     2                            stress2,statev2,de33 )
c
c-----------------------------------------------------------------------
      implicit none
c
      common /ummdp1/ne,ip,lay
      common /ummdp2/prop
      common /ummdp3/nsdv
      common /ummdp4/propdim
      common /ummdpa/n1234
c
      integer,intent(in) :: ndi,nshr,ntens,nstatev,nprops,noel,kinc
      real*8 ,intent(in) :: props(nprops)
      real*8 ,intent(in) :: stress1(ntens),statev1(nstatev),
     1                      dstrain(ntens)
c
      real*8 ,intent(out) :: de33
      real*8 ,intent(out) :: stress2(ntens),statev2(nstatev)
c
      integer n1234,mxpbs,mxprop,nrot
      parameter (mxpbs=10,mxprop=100)
      integer ne,ip,lay,nsdv,propdim,i,k,n,is,nprop,nvbs0,nvbs,ndela,
     1        ndyld,ndihd,ndkin,npbs,ndrup,mjac,isvrsvd,isvsclr,maxsdv
      real*8 p,dp
      real*8 ustatev(ntens),s2(ntens),dpe(ntens),pe(ntens),prop(mxprop)
      real*8 x1(mxpbs,ntens),x2(mxpbs,ntens),ddsdde(ntens,ntens)
      character*100 text
c-----------------------------------------------------------------------
c
cf2py intent(in) stress1,statev1,dstrain
cf2py intent(in) ndi,nshr,ntens,nstatev
cf2py intent(in) props,nprops
cf2py intent(in) noel,kinc
cf2py intent(out) stress2,statev2,de33
cf2py depend(ntens) stress1,stress2,dstrain
cf2py depend(nstatev) statev1,statev2
cf2py depend(nprops) props
c
c-----------------------------------------------------------------------
c
      ne = noel
      ip = 1
      lay = 1
      nsdv = nstatev
      nprop = mxprop
      propdim = nprops - 1
c
c                                        ---- set debug and verbose mode
      nvbs0 = 0
      call ummdp_debugmode ( nvbs,nvbs0 )
c
c                                           ---- set material properties
      do i = 1,nprops
        prop(i) = props(i)
      end do
c
      call ummdp_prop_dim ( prop,nprop,propdim,ndela,ndyld,ndihd,ndkin,
     1                      npbs,ndrup )
c
      if ( npbs > mxpbs ) then
        write (6,*) 'npbs > mxpbs error in umat'
        write (6,*) 'npbs =',npbs
        write (6,*) 'mxpbs=',mxpbs
        call ummdp_exit ( 301 )
      end if
c
c                                                     ---- check nstatev
      call ummdp_check_nisv ( nstatev,ntens,npbs )
c
c                             ---- copy current internal state variables
      call ummdp_isvprof ( isvrsvd,isvsclr )
      call ummdp_isv2pex ( isvrsvd,isvsclr,statev1,nstatev,p,pe,x1,
     1                     ntens,mxpbs,npbs )
c
c                             ---- update stress and set tangent modulus
      mjac = 0
      call ummdp_plasticity ( stress1,stress2,dstrain,p,dp,dpe,de33,x1,
     1                        x2,mxpbs,ddsdde,ndi,nshr,ntens,nvbs,mjac,
     2                        prop,nprop,propdim )
c
c                                  ---- update equivalent plastic strain
      statev2(isvrsvd+1) = p + dp
c
c                                  ---- update plastic strain components
      do i = 1,ntens
        is = isvrsvd + isvsclr + i
        statev2(is) = statev1(is) + dpe(i)
      end do
c
c                                     ---- update back stress components
      if ( npbs /= 0 ) then
        do n = 1,npbs
          do i = 1,ntens
            is = isvrsvd + isvsclr + ntens*n + i
            statev2(is) = x2(n,i)
          end do
        end do
      end if
c
      return
      end subroutine plug_ummdp_vfm
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
c
      integer,intent(in) :: nexit
c
      integer ne,ip,lay
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
      call f2py_stop ( )
c
      return
      end subroutine ummdp_exit
c
c
c
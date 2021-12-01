************************************************************************
*
*     PRINT SUBROUTINES
*
************************************************************************
c
c     ummdp_print_ummdp ( )
c       print ummmdp separator
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
c     PRINT UMMDP SEPARATOR
c
      subroutine ummdp_print_ummdp ( )
c-----------------------------------------------------------------------
      implicit none
c
      character*20 fmt
c-----------------------------------------------------------------------
c
      fmt = '(1/,4xA,1/)'
      write (6,fmt) '~~~~~~~~~~~~~~~~~~ UMMDp ~~~~~~~~~~~~~~~~~~'
c
      return
      end subroutine ummdp_print_ummdp
c
c
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
      integer ntela,i
      character*50 fmtid,fmtpr
c-----------------------------------------------------------------------
c
      fmtid = '(16xA,I1)'
      fmtpr = '(16xA10,I1,A3,E20.12)'
c
      write (6,'(/12xA)') '> Elasticity'
c
      ntela = nint(prela(1))
      select case ( ntela )
      case ( 0 )
        write (6,fmtid) '. Young Modulus & Poisson Ratio | ',ntela
      case ( 1 )
        write (6,fmtid) '. Bulk Modulus & Modulus of Rigidity | ',ntela
      end select
c
      do i = 1,ndela-1
        write (6,fmtpr) '. prela(1+',i,') =',prela(i+1)
      end do
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
      character*50 fmtid1,fmtid2,fmtpr
c-----------------------------------------------------------------------
c
      fmtid1 = '(16xA,I1)'
      fmtid2 = '(16xA,I2)'
      fmtpr = '(16xA10,I1,A3,E20.12)'
c
      write (6,'(/12XA)') '> Yield Function'
c
      ntyld = pryld(1)
      select case ( ntyld )
      case ( 0 )
        write (6,fmtid1) '. von Mises | ',ntyld
      case ( 1 )
        write (6,fmtid1) '. Hill 1948 | ',ntyld
      case ( 2 )
        write (6,fmtid1) '. Yld2004-18p | ',ntyld
      case ( 3 )
        write (6,fmtid1) '. CPB 2006 | ',ntyld
      case ( 4 )
        write (6,fmtid1) '. Karafillis-Boyce 1993 | ',ntyld
      case ( 5 )
        write (6,fmtid1) '. Hu 2005 | ',ntyld
      case ( 6 )
        write (6,fmtid1) '. Yoshida 2011 | ',ntyld
c
      case ( -1 )
        write (6,fmtid2) '. Gotoh | ',ntyld
      case ( -2 )
        write (6,fmtid2) '. Yld2000-2d | ',ntyld
      case ( -3 )
        write (6,fmtid2) '. Vegter | ',ntyld
      case ( -4 )
        write (6,fmtid2) '. BBC 2005 | ',ntyld
      case ( -5 )
        write (6,fmtid2) '. Yld89 | ',ntyld
      case ( -6 )
        write (6,fmtid2) '. BBC 2008 | ',ntyld
      case ( -7 )
        write (6,fmtid2) '. Hill 1990 | ',ntyld
      end select
c
      select case ( ntyld )
      case ( -3 )                                               ! Vegter
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
      case ( -6 )                                             ! BBC 2008
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
      case default
        do i = 1,ndyld-1
          write (6,fmtpr) '. pryld(1+',i,') =',pryld(i+1)
        end do
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
      integer ntihd,i
      character*50 fmtid,fmtpr
c-----------------------------------------------------------------------
c
      fmtid = '(16xA,I1)'
      fmtpr = '(16xA10,I1,A3,E20.12)'
c
      write (6,'(/12xA)') '> Isotropic Hardening Law'
c
      ntihd = nint(prihd(1))
      select case ( ntihd )
      case ( 0 )
        write (6,fmtid) '. Perfect Plasticity | ',ntihd
      case ( 1 )
        write (6,fmtid) '. Linear | ',ntihd
      case ( 2 )
        write (6,fmtid) '. Swift | ',ntihd
      case ( 3 )
        write (6,fmtid) '. Ludwick | ',ntihd
      case ( 4 )
        write (6,fmtid) '. Voce | ',ntihd
      case ( 5 )
        write (6,fmtid) '. Voce & Linear | ',ntihd
      case ( 6 )
        write (6,fmtid) '. Voce & Swift | ',ntihd
      end select
c
      do i = 1,ndihd-1
        write (6,fmtpr) '. prihd(1+',i,') =',prihd(i+1)
      end do
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
      character*50 fmtid,fmtpr
c-----------------------------------------------------------------------
c
      fmtid = '(16xA,I1)'
      fmtpr = '(16xA10,I1,A3,E20.12)'
c
      write (6,'(/12xA)') '> Kinematic Hardening Law'
c
      ntkin = nint(prkin(1))
      select case ( ntkin )
      case ( 0 )
        write (6,fmtid) '. None | ',ntkin
      case ( 1 )
        write (6,fmtid) '. Prager | ',ntkin
      case ( 2 )
        write (6,fmtid) '. Ziegler | ',ntkin
      case ( 3 )
        write (6,fmtid) '. Armstrong-Frederick | ',ntkin
      case ( 4 )
        write (6,fmtid) '. Chaboche I | ',ntkin
      case ( 5 )
        write (6,fmtid) '. Chaboche II | ',ntkin
      case ( 6 )
        write (6,fmtid) '. Yoshida-Uemori | ',ntkin
      end select
c
      do i = 1,ndkin-1
        write (6,fmtpr) '. prkin(1+',i,') =',prkin(i+1)
      end do
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
      integer ntrup,i
      character*50 fmtid,fmtpr
c-----------------------------------------------------------------------
c
      fmtid = '(16xA,I1)'
      fmtpr = '(16xA10,I1,A3,E20.12)'
c
      write (6,'(/12xA)') '>> Uncoupled Rupture Criterion'
c
      ntrup = nint(prrup(1))
      select case ( ntrup )
      case ( 0 )
        write (6,fmtid) '. None | ',ntrup
      case ( 1 )
        write (6,fmtid) '. Equivalent Plastic Strain | ',ntrup
      case ( 2 )
        write (6,fmtid) '. Cockroft and Latham | ',ntrup
      case ( 3 )
        write (6,fmtid) '. Rice and Tracey | ',ntrup
      case ( 4 )
        write (6,fmtid) '. Ayada | ',ntrup
      case ( 5 )
        write (6,fmtid) '. Brozzo | ',ntrup
      case ( 6 )
        write (6,fmtid) '. Forming Limit Diagram | ',ntrup
      end select
c
      do i = 1,ndrup-1
        write (6,fmtpr) '. prrup(1+',i,') =',prrup(i+1)
      end do
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
      common /ummdp1/ne,ip,lay
c
			integer inc,nnrm,nshr
c
			integer ne,ip,lay,nttl,nerr
      character*50 fmt1,fmt2,fmt3,fmt4,ptype,tmp
c-----------------------------------------------------------------------
      fmt1 = '(/12xA,A)'
      fmt2 =  '(12xA,A)'
      fmt3 = '(/12xA,I1)'
      fmt4 =  '(12xA,I1)'
c
      nttl = nnrm + nshr
c
      write(6,'(4/8xA)') '>> Info'
c
      write (tmp,'(I)') inc
      write (6,fmt1) '        Increment : ',adjustl(tmp)
c
      write (tmp,'(I)') ne
      write (6,fmt1) '          Element : ',adjustl(tmp)
      write (tmp,'(I)') ip
      write (6,fmt2) 'Integration Point : ',adjustl(tmp)
      write (tmp,'(I)') lay
      write (6,fmt2) '            Layer : ',adjustl(tmp)
! c
      write (6,fmt3) ' Total Components : ',nttl
      write (6,fmt4) 'Normal Components : ',nnrm
      write (6,fmt4) ' Shear Components : ',nshr
! c
      nerr = 0
      if ( nnrm == 3 ) then
        if ( nshr == 3 ) then
          ptype = '3D Solid'
        else if ( nshr == 1 ) then
          ptype = 'Plane Strain or Axisymmetric Solid'
        else
          nerr = 1
        end if
      else if ( nnrm == 2 ) then
        if ( nshr == 1 ) then
          ptype = 'Plane Stress or Thin Shell'
        else if ( nshr == 3 ) then
          ptype = 'Thick Shell'
        else
          nerr = 1
        end if
      else
        nerr = 1
      end if
c
      if ( nerr == 0 ) then
        write (6,'(/12xA,A)') '     Element Type : ',ptype
      else
        ptype = 'Not Supported'
        write (6,'(/12xA,A)') '     Element Type : ',ptype
        call ummdp_exit ( 100 )
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
      character*100 text
c-----------------------------------------------------------------------
c
      if ( io == 0 ) then
        write(6,'(//8xA)') '>> Input'
      else
        write(6,'(//8xA)') '>> Output'
      end if
c
      text = 'Stress'
      call ummdp_utility_print1 ( text,s,nttl,0 )
c
      text = 'Internal State Variables'
      call ummdp_utility_print1 ( text,stv,nstv,0 )
c
      text = 'Strain Increment'
      call ummdp_utility_print1 ( text,de,nttl,0 )
c
      text = 'Tangent Modulus Matrix'
      call ummdp_utility_print2 ( text,d,nttl,nttl,0 )
c
      return
      end subroutine ummdp_print_inout
c
c
c
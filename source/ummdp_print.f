c***********************************************************************
c
c     UMMDp: Print Subroutines
c
c***********************************************************************
c
c     jancae_elast_print ( prela,ndela )
c       print elasticity parameters
c
c     jancae_yfunc_print ( pryld,ndyld )
c       print yield criteria parameters
c
c     jancae_harden_print ( prihd,ndihd )
c       print isotropic hardening law parameters
c
c     jancae_kinematic_print ( prkin,ndkin,npbs )
c       print kinematic hardening law parameters
c
c     jancae_rupture_print ( prrup,ndrup )
c       print uncoupled rupture criterion parameters
c     
************************************************************************
c     PRINT ELASTICITY PARAMETERS
c
      subroutine jancae_elast_print ( prela,ndela )
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
      write (6,*) '>> Elasticity',ntela
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
      end subroutine jancae_elast_print
c
c
c
************************************************************************
c     PRINT YIELD CRITERIA PARAMETERS
c
      subroutine jancae_yfunc_print ( pryld,ndyld )
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
      write (6,*) '>> Yield Criterion',ntyld
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
      end subroutine jancae_yfunc_print
c
c
c
************************************************************************
c     PRINT ISOTROPIC HARDENING LAW PARAMETERS
c
      subroutine jancae_harden_print ( prihd,ndihd )
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
      end subroutine jancae_harden_print
c
c
c
************************************************************************
c     PRINT KINEMATIC HARDENING LAW PARAMETERS
c
      subroutine jancae_kinematic_print ( prkin,ndkin,npbs )
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
      end subroutine jancae_kinematic_print
c
c
c
************************************************************************
c     PRINT UNCOUPLED RUPTURE CRITERION PARAMETERS
c
      subroutine jancae_rupture_print ( prrup,ndrup )
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
      end subroutine jancae_rupture_print
c
c
c
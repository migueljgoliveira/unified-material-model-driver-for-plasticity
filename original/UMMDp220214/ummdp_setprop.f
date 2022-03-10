c-------------------------------------------------------------
c     set material parameters
c
      subroutine jancae_set_prop ( prop,nprop )
c-------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension prop(nprop)
      logical   maninp
c
      maninp=.false.  ! manual input
c
c
      call jancae_clear1 ( prop,nprop )
      n0=0
c                                  ---------------- elastic properties
c
      nela=0
      prop(n0+1  )=float(nela)         ! Isotropic elasticity
      prop(n0+1+1)=200.0d0*1000.0d0    !   E  young's modulus
      prop(n0+1+2)=0.3d0               !   nu poisson's ratio
      ndela=3
c
c*****160130 TUAT Nagai spherical deep drawing
      prop(n0+1+1)=69.20d0*1000.0d0  !   E  young's modulus
      prop(n0+1+2)= 0.33d0           !   nu poisson's ratio
      ndela=3
c
      n0=n0+ndela
c
c
c                                      ---------------- yield function
      nyf=2
c
   10 if ( maninp ) then
        write (6,1010)
        write (6,*   ) 'type of yield function'
        read  (5,*   ) nyf
      endif
 1010 format( ' 0 : von Mises isotropic (1913)',/,
     &        ' 1 : Hill quadratic (1948)',/,
     &        ' 2 : Barlat yld2004 (2005)',/,
     &        ' 3 : Cazacu (2006)',/,
     &        ' 4 : Karafillis-Boyce (1993)',/,
     &        ' 5 : Hu (2005)',/,
     &        ' 6 : Yoshida (2011)',/,
     &        '-1 : Gotoh biquadratic (1978)',/,
     &        '-2 : Barlat YLD2000-2d (2000)',/,
     &        '-3 : Vegter',/,
     &        '-4 : Banabic BBC2005',/,
     &        '-5 : Barlat YLD89',/,
     &        '-6 : Banabic BBC2008',/,
     &        '-7 : Hill 1990',/)
c
c
c
      select case ( nyf )
c
      case ( 0 )
      prop(n0+1   )=float(0) ! von Mises (1913)
      ndyld=1
c
      case ( 1 )
      prop(n0+1   )=float(1) ! Hill quadratic (1948)
      prop(n0+1+1)=1.0d0       !   F
      prop(n0+1+2)=1.0d0       !   G  
      prop(n0+1+3)=1.0d0       !   H
      prop(n0+1+4)=3.0d0       !   L
      prop(n0+1+5)=3.0d0       !   M
      prop(n0+1+6)=3.0d0       !   N
c                                     * 151220 TUAT HoleExp
      prop(n0+1+1)=0.806105214d0 ! F
      prop(n0+1+2)=0.833482177d0 ! G
      prop(n0+1+3)=1.166583671d0 ! H
      prop(n0+1+6)=2.806543046d0 ! N
c
      ndyld=7
c
      case ( 2 )
      prop(n0+1   )=float(2) ! Barlat yld2004-18p (2005)
c     sample : IJP v.21(2005) p1009-1039. 2090-T3
      prop(n0+1+ 1)=-0.069888d0  ! c'12
      prop(n0+1+ 2)= 0.936408d0  ! c'13
      prop(n0+1+ 3)= 0.079143d0  ! c'21
      prop(n0+1+ 4)= 1.003060d0  ! c'23
      prop(n0+1+ 5)= 0.524741d0  ! c'31
      prop(n0+1+ 6)= 1.363180d0  ! c'32
      prop(n0+1+ 7)= 0.954322d0  ! c'44
      prop(n0+1+ 8)= 1.023770d0  ! c'55
      prop(n0+1+ 9)= 1.069060d0  ! c"66
      prop(n0+1+10)= 0.981171d0  ! c"12
      prop(n0+1+11)= 0.476741d0  ! c"13
      prop(n0+1+12)= 0.575316d0  ! c"21
      prop(n0+1+13)= 0.866827d0  ! c"23
      prop(n0+1+14)= 1.145010d0  ! c"31
      prop(n0+1+15)=-0.079294d0  ! c"32
      prop(n0+1+16)= 1.404620d0  ! c"44
      prop(n0+1+17)= 1.051660d0  ! c"55
      prop(n0+1+18)= 1.147100d0  ! c"66
      prop(n0+1+19)=float(8)     ! order
c
c
c*****160130 TUAT Nagai Case-1(Case-1b)
      do i=1,18
        prop(n0+1+i)=1.0d0
      enddo
      prop(n0+1+ 1)= 0.567726464d0  ! c'12
      prop(n0+1+ 2)= 1.011906010d0  ! c'13
      prop(n0+1+ 3)= 0.918738294d0  ! c'21
      prop(n0+1+ 4)= 0.949253189d0  ! c'23
      prop(n0+1+ 5)= 0.804482945d0  ! c'31
      prop(n0+1+ 6)= 0.541426215d0  ! c'32
c     c'44=c'55=1
      prop(n0+1+ 9)= 0.641800051d0  ! c'66
      prop(n0+1+10)= 0.881538993d0  ! c"12
      prop(n0+1+11)= 1.281298612d0  ! c"13
      prop(n0+1+12)= 1.468510772d0  ! c"21
      prop(n0+1+13)= 1.694334282d0  ! c"23
      prop(n0+1+14)= 0.874665642d0  ! c"31
      prop(n0+1+15)= 0.562106847d0  ! c"32
c     c"44=c"55=1
      prop(n0+1+18)= 1.275309821d0  ! c"66
      prop(n0+1+19)= 6.565182435d0  ! exponent a
c
      ndyld=20
c
c
      case ( 3 )
      prop(n0+1   )=float(3) ! Cazacu (2006)
c     sample :  IJP v.22(2006) p1171-1194. Mg-Th 1%
      prop(n0+1+ 1)= 1.0000d0      ! C11
      prop(n0+1+ 2)= 0.4802d0      ! C12
      prop(n0+1+ 3)= 0.2592d0      ! C13
      prop(n0+1+ 4)= prop(n0+1+ 2) ! C21
      prop(n0+1+ 5)= 0.9517d0      ! C22
      prop(n0+1+ 6)= 0.2071d0      ! C23
      prop(n0+1+ 7)= prop(n0+1+ 3) ! C31
      prop(n0+1+ 8)= prop(n0+1+ 6) ! C32
      prop(n0+1+ 9)= 0.4654d0      ! C33
      prop(n0+1+10)= 1.0000d0      ! C44
      prop(n0+1+11)= 1.0000d0      ! C55
      prop(n0+1+12)= 1.0000d0      ! C66
      prop(n0+1+13)= 2.0000d0      ! a
      prop(n0+1+14)= 0.3539d0      ! k
      ndyld=15
c                                     * 151220 TUAT HoleExp
      prop(n0+1+ 1)=  0.916990888d0      ! C11
      prop(n0+1+ 2)= -0.150679720d0      ! C12
      prop(n0+1+ 3)= -0.438968742d0      ! C13
      prop(n0+1+ 4)=  prop(n0+1+ 2)      ! C21
      prop(n0+1+ 5)=  0.989590667d0      ! C22
      prop(n0+1+ 6)=  0.237646841d0      ! C23
      prop(n0+1+ 7)=  prop(n0+1+ 3)      ! C31
      prop(n0+1+ 8)=  prop(n0+1+ 6)      ! C32
      prop(n0+1+ 9)= -1.009587128d0      ! C33
      prop(n0+1+10)=  1.0000d0           ! C44
      prop(n0+1+11)=  1.0000d0           ! C55
      prop(n0+1+12)=  1.165424349d0      ! C66
      prop(n0+1+13)=  1.900988379d0      ! a
      prop(n0+1+14)= -0.068235547d0      ! k
c
c
      case ( 4 )
      prop(n0+1   )=float(4)       ! Karafillis-Boyce (1993)
      ndyld=1+8
      prop(n0+1+ 1)=0.599 ! C
      prop(n0+1+ 2)=1.103 ! alpha1
      prop(n0+1+ 3)=1.120 ! alpha2
      prop(n0+1+ 4)=1.0   ! gamma1
      prop(n0+1+ 5)=1.0   ! gamma2
      prop(n0+1+ 6)=0.730 ! gamma3
      prop(n0+1+ 7)=15.0  ! k
      prop(n0+1+ 8)=0.810 ! c
c
c     c=0.599
c     alpha1=1.103
c     alpha2=1.120
c     gamma3=0.730
c     beta1=(alpha2-alpha1-1d0)*0.5d0
c     beta2=(alpha1-alpha2-1d0)*0.5d0
c     beta3=(1d0-alpha1-alpha2)*0.5d0
c     prop(n0+1+ 1)=c           ! L(1,1)
c     prop(n0+1+ 2)=c*beta1     ! L(1,2)
c     prop(n0+1+ 3)=c*beta2     ! L(1,3)
c     prop(n0+1+ 7)=c*alpha1    ! L(2,2)
c     prop(n0+1+ 8)=c*beta3     ! L(2,3)
c     prop(n0+1+12)=c*alpha2    ! L(3,3)
c     prop(n0+1+16)=1d0         ! L(4,4)
c     prop(n0+1+19)=1d0         ! L(5,5)
c     prop(n0+1+21)=c*gamma3    ! L(6,6)
c     prop(n0+1+22)=float(15)   ! k
c     prop(n0+1+23)=0.81d0      ! c
c
      case ( 5 )
      prop(n0+1   )=float(5)    ! Hu (2005)
      ndyld=1+9
c                               ! coefs for von Mises 
      prop(n0+1+ 1)= 1.000d0    ! X1
      prop(n0+1+ 2)=-2.000d0    ! X2
      prop(n0+1+ 3)= 3.000d0    ! X3
      prop(n0+1+ 4)=-2.000d0    ! X4
      prop(n0+1+ 5)= 1.000d0    ! X5
      prop(n0+1+ 6)= 9.000d0    ! X7
      prop(n0+1+ 7)= 6.000d0    ! C1
      prop(n0+1+ 8)= 6.000d0    ! C2
      prop(n0+1+ 9)= 6.000d0    ! C3
c
      case ( 6 )
      prop(n0+1   )=float(6)    ! Yoshida (2005)
      ndyld=1+16
c     sample : IJP v.45 (2013) p119-139 590HSS 
      prop(n0+1+ 1)= 1.0000d0      ! C1
      prop(n0+1+ 2)= 0.6014d0      ! C2
      prop(n0+1+ 3)= 0.4841d0      ! C3
      prop(n0+1+ 4)= 0.3982d0      ! C4
      prop(n0+1+ 5)= 0.4375d0      ! C5
      prop(n0+1+ 6)= 0.5752d0      ! C6
      prop(n0+1+ 7)= 0.7591d0      ! C7
      prop(n0+1+ 8)= 1.0350d0      ! C8
      prop(n0+1+ 9)= 0.6789d0      ! C9
      prop(n0+1+10)= 0.6664d0      ! C10
      prop(n0+1+11)= 0.7540d0      ! C11
      prop(n0+1+12)= 0.9423d0      ! C12
      prop(n0+1+13)= 1.1601d0      ! C13
      prop(n0+1+14)= 1.0615d0      ! C14
      prop(n0+1+15)= 1.2470d0      ! C15
      prop(n0+1+16)= 1.7732d0      ! C16
c                                     * 151220 TUAT HoleExp
      prop(n0+1+ 1)= 0.999670012d0      ! C1
      prop(n0+1+ 2)= 1.153640645d0      ! C2
      prop(n0+1+ 3)= 1.242357637d0      ! C3
      prop(n0+1+ 4)= 1.365842839d0      ! C4
      prop(n0+1+ 5)= 1.272444524d0      ! C5
      prop(n0+1+ 6)= 1.180478701d0      ! C6
      prop(n0+1+ 7)= 0.981892260d0      ! C7
      prop(n0+1+ 8)= 1.035989203d0      ! C8
      prop(n0+1+ 9)= 1.153577495d0      ! C9
      prop(n0+1+10)= 0.989009593d0      ! C10
      prop(n0+1+11)= 1.216987726d0      ! C11
      prop(n0+1+12)= 1.099239973d0      ! C12
      prop(n0+1+13)= 1.066042025d0      ! C13
      prop(n0+1+14)= 0.932677556d0      ! C14
      prop(n0+1+15)= 1.195813771d0      ! C15
      prop(n0+1+16)= 0.841010957d0      ! C16
c
      ndyld=1+16
c
c
      case ( -1 )
      prop(n0+1   )=float(-1) ! Gotoh biquadratic (1978)
c     sample : J.JSTP v19 no210 p599 Al-killed steel
      prop(n0+1+1)= 1.00d0     !   A1 (eq.strs=ut in x : 1.0)
      prop(n0+1+2)=-2.60d0    !   A2  
      prop(n0+1+3)= 3.75d0    !   A3
      prop(n0+1+4)=-2.79d0    !   A4
      prop(n0+1+5)= 0.99d0    !   A5
      prop(n0+1+6)= 6.29d0    !   A6
      prop(n0+1+7)=-7.72d0    !   A7
      prop(n0+1+8)= 6.33d0    !   A8
      prop(n0+1+9)= 8.96d0    !   A9
      ndyld=10
c
      case ( -2 )
      prop(n0+1   )=float(-2) ! Barlat yld2000-2d (2003)
c     sample : IJP v.19(2003)p1297-1319 2090-T3 Exp.rb
      prop(n0+1+1)=0.4865d0   !   alpha1
      prop(n0+1+2)=1.3783d0   !   alpha2
      prop(n0+1+3)=0.7536d0   !   alpha3
      prop(n0+1+4)=1.0246d0   !   alpha4
      prop(n0+1+5)=1.0363d0   !   alpha5
      prop(n0+1+6)=0.9036d0   !   alpha6
      prop(n0+1+7)=1.2321d0   !   alpha7
      prop(n0+1+8)=1.4858d0   !   alpha8
      prop(n0+1+9)=float(8)   !   order a 
c                                    * 151220 TUAT HoleExp
      prop(n0+1+1)=0.973737830d0
      prop(n0+1+2)=1.062061767d0
      prop(n0+1+3)=0.843006204d0
      prop(n0+1+4)=0.927158106d0
      prop(n0+1+5)=0.941646566d0
      prop(n0+1+6)=0.776575885d0
      prop(n0+1+7)=0.983219298d0
      prop(n0+1+8)=1.121953229d0
      prop(n0+1+9)=4.893604841d0
c
      ndyld=10
c
      case ( -3 )
      prop(n0+1   )=float(-3)       ! Vegter (2006) almost pass
c     sample : AL 3 tests (0,45,90deg)
      nf=3
      prop(n0+1+1)= float(nf)         ! nf
      prop(n0+1+2)= 1.004000d0        ! f_bi0
      prop(n0+1+3)= 0.889000d0        ! r_bi0
      i=0
      prop(n0+1+3+i*4+1)=  1.001000d0  ! phi_un(0)
      prop(n0+1+3+i*4+2)=  0.600000d0  ! phi_sh(0)
      prop(n0+1+3+i*4+3)=  1.045750d0  ! phi_ps(0)
      prop(n0+1+3+i*4+4)= -0.367678d0  ! omg(0)
      i=1
      prop(n0+1+3+i*4+1)=  0.006000d0  ! phi_un(1)
      prop(n0+1+3+i*4+2)=  0.000000d0  ! phi_sh(1)
      prop(n0+1+3+i*4+3)=  0.006500d0  ! phi_ps(1)
      prop(n0+1+3+i*4+4)=  0.020787d0  ! omg(1)
      i=2
      prop(n0+1+3+i*4+1)=  0.014000d0  ! phi_un(2)
      prop(n0+1+3+i*4+2)= -0.040000d0  ! phi_sh(2)
      prop(n0+1+3+i*4+3)=  0.008750d0  ! phi_ps(2)
      prop(n0+1+3+i*4+4)= -0.043353d0  ! omg(2)
      ndyld=1+3+nf*4
c     sample IF steel 3 tests (0,45,90deg)
c     nf=3
c     prop(n0+1+ 1)= float(nf)        ! nf
c     prop(n0+1+ 2)= 1.157000d0       ! f_bi0
c     prop(n0+1+ 3)= 0.777000d0       ! r_bi0
c     i=0
c     prop(n0+1+3+i*4+1)= 0.999250d0  ! phi_un(0)
c     prop(n0+1+3+i*4+2)= 0.541000d0  ! phi_sh(0)
c     prop(n0+1+3+i*4+3)= 1.250250d0  ! phi_ps(0)
c     prop(n0+1+3+i*4+4)=-0.677657d0  ! omg(0)
c     i=1
c     prop(n0+1+3+i*4+1)= 0.003500d0  ! phi_un(1)
c     prop(n0+1+3+i*4+2)= 0.000000d0  ! phi_sh(1)
c     prop(n0+1+3+i*4+3)=-0.001500d0  ! phi_ps(1)
c     prop(n0+1+3+i*4+4)= 0.032988d0  ! omg(1)
c     i=2
c     prop(n0+1+3+i*4+1)= 0.001250d0  ! phi_un(2)
c     prop(n0+1+3+i*4+2)=-0.004000d0  ! phi_sh(2)
c     prop(n0+1+3+i*4+3)=-0.001750d0  ! phi_ps(2)
c     prop(n0+1+3+i*4+4)= 0.032988d0  ! omg(2)
c     ndyld=1+3+(nf+1)*4
c
      case ( -4 )
      prop(n0+1  )=float(-4)      ! Banabic (BBC2005)  
c                                   sample : von Mises
      k=2                         ! order 2*k
      prop(n0+1+1)=float(k)
      prop(n0+1+2)=0.5d0          ! a
      prop(n0+1+3)=0.5d0          ! b
      prop(n0+1+4)=1.0d0          ! L
      prop(n0+1+5)=1.0d0          ! M
      prop(n0+1+6)=1.0d0          ! N
      prop(n0+1+7)=1.0d0          ! P
      prop(n0+1+8)=1.0d0          ! Q
      prop(n0+1+9)=1.0d0          ! R
      ndyld=1+9
c
      case ( -5 )
      prop(n0+1  )=float(-5)      ! Barlat (yld89) 
      prop(n0+1+1)=float(8)       ! order
      prop(n0+1+2)=1.0d0          ! a
      prop(n0+1+3)=2.0d0/3.0d0    ! h
      prop(n0+1+4)=1.0d0          ! p
c     prop(n0+1+1)=float(2)       ! order
c     prop(n0+1+2)=1.0d0          ! a
c     prop(n0+1+3)=1.0d0          ! h
c     prop(n0+1+4)=1.0d0          ! p
      ndyld=1+4
c
      case ( -6 )
      prop(n0+1   ) = float(-6)      ! Banabic (BBC2008)
c                                      sample : AA2090-T3 (s=2)
      prop(n0+1+ 1) = float( 2)      ! s
      prop(n0+1+ 2) = float( 4)      ! k
c     --------------------------
      prop(n0+1+ 3) = 0.130866d0     ! l_1
      prop(n0+1+ 4) = 0.621742d0     ! l_2
      prop(n0+1+ 5) = 0.783422d0     ! m_1     
      prop(n0+1+ 6) = 0.660402d0     ! m_2
      prop(n0+1+ 7) = 0.000079d0     ! m_3 
      prop(n0+1+ 8) = 0.110991d0     ! n_1 
      prop(n0+1+ 9) = 0.048245d0     ! n_2 
      prop(n0+1+10) = 0.307522d0     ! n_3 
c     -------------------------
      prop(n0+1+11) = 1.033922d0     ! l_1
      prop(n0+1+12) =-0.071963d0     ! l_2
      prop(n0+1+13) = 0.000113d0     ! m_1 
      prop(n0+1+14) = 0.000077d0     ! m_2 
      prop(n0+1+15) = 0.538047d0     ! m_3 
      prop(n0+1+16) = 0.055764d0     ! n_1 
      prop(n0+1+17) = 1.018603d0     ! n_2 
      prop(n0+1+18) = 0.778150d0     ! n_3 
c     --------------------------
      ndyld=1+18
c
      case ( -7 )
      prop(n0+1   )=float(-7)     ! R.Hill (Hill90)   not pass
c          sample : J.JSTP v40(n457)pp145-149(1999)
      am    =   2.4d0
      sig0  = 180.0d0
      sig45 = 188.0d0/sig0
      sig90 = 184.0d0/sig0
      sigb  = 184.0d0/sig0
      sig0  =   1.0d0
c     am    =   2.0d0        ! mises model   pass
c     sig0  =   1.0d0
c     sig45 =   1.0d0/sig0
c     sig90 =   1.0d0/sig0
c     sigb  =   1.0d0/sig0
c     sig0  =   1.0d0
      a=0.25d0*((2.0d0*sigb/sig90)**am-(2.0d0*sigb/sig0 )**am)
      b=0.50d0*((2.0d0*sigb/sig0 )**am+(2.0d0*sigb/sig90)**am)
     &        - (2.0d0*sigb/sig45)**am
      c     = (2.0d0*sigb/sig45 )**am-1.0d0
      tau   = sigb/(c**(1.0d0/am) )
      prop(n0+1+1)=a        ! a
      prop(n0+1+2)=b        ! b
      prop(n0+1+3)=tau      ! tau
      prop(n0+1+4)=sigb     ! sigmab
      prop(n0+1+5)=am       ! order M
      ndyld=1+5
c
      case default
        write (6,*) 'error retype'
        goto 10
      end select
      n0=n0+ndyld
c
c
c
c                                 ---------------- isotropic hardening
      nih=5
   20 if ( maninp ) then
        write (6,1020)
        write (6,*   ) 'type of isotropic hardening'
        read  (5,*   ) nih
      endif
 1020 format( ' 0 : perfectly  sy0',/,
     &        ' 1 : linear     sy0+h*p',/,
     &        ' 2 : Swift      c*(e0+p)^en',/,
     &        ' 3 : Ludwick    sy0+c*p^en',/,
     &        ' 4 : Voce       sy0+q*(1-exp(-b*p))',/,
     &        ' 5 : Voce+lin.  sy0+q*(1-exp(-b*p)))+c',/,
     &        ' 6 : Voce+Swift ....')
c
      select case ( nih )
      case ( 0 )                      
      prop(n0+1  )=float(0)  ! perfectly plastic
      prop(n0+1+1)=1.0d0     !   sy0
      ndihd=2
c
      case ( 1 )
      prop(n0+1   )=float(1) ! linear  sy0+h*p
      prop(n0+1+1)=1.0d0     !   sy0
      prop(n0+1+2)=1.0d0     !   h
      ndihd=3
c
      case ( 2 )
      prop(n0+1  )=float(2)  ! Swift   c*(e0+p)^en
      prop(n0+1+1)=1.0d0     !   c
      prop(n0+1+2)=1.0d0     !   e0
      prop(n0+1+3)=1.0d0     !   en
c                                    * 151221 TUAT
      prop(n0+1+1)=541.0d0     !   c
      prop(n0+1+2)=0.0036d0    !   e0
      prop(n0+1+3)=0.249d0     !   en
      ndihd=4
c
      case ( 3 )
      prop(n0+1  )=float(3)  ! Ludwick sy0+c*p^en
      prop(n0+1+1)=1.0d0     !   sy0
      prop(n0+1+2)=1.0d0     !   c
      prop(n0+1+3)=1.0d0     !   en
      ndihd=4
c
      case ( 4 )
      prop(n0+1  )=float(4)  ! Voce    sy0+q*(1-exp(-b*p))
      prop(n0+1+1)=1.0d0     !   sy0
      prop(n0+1+2)=1.0d0     !   q
      prop(n0+1+3)=1.0d0     !   b
      ndihd=4
c
      case ( 5 )
      prop(n0+1  )=float(5)  ! Voce+lin. sy0+q*(1-exp(-b*p))+c*p
      prop(n0+1+1)=1.0d0     !   sy0
      prop(n0+1+2)=1.0d0     !   q
      prop(n0+1+3)=1.0d0     !   b
      prop(n0+1+4)=1.0d0     !   c
      ndihd=5
c*****160130 TUAT Nagai
      prop(n0+1+1)=315.1d0-189.1d0     !   sy0
      prop(n0+1+2)=        189.1d0     !   q
      prop(n0+1+3)=         13.1d0     !   b
      prop(n0+1+4)=220.7d0             !   c
      ndihd=5
c
      case ( 6 )
      prop(n0+1  )=float(6) ! Voce+Swift a*(sy0+q*(1-exp(-b*p))+(1-a)*(c*(e0+p)**n) 
      write (6,*) 'nih=6 set here'
      stop
c
      case default
        write (6,*) 'error retype'
        goto 20
      end select
      n0=n0+ndihd
c
c
c
c                                 ---------------- kinematic hardening
      nkh=0
   30 if ( maninp ) then   
        write (6,1030)
        write (6,*   ) 'type of kinematic hardening'
        read  (5,*   ) nkh
      endif
 1030 format( ' 0 : no kinematic hardening',/,
     &        ' 1 : Prager(1949)',/,
     &        ' 2 : Ziegler (1959)',/,
     &        ' 3 : Armstrong-Frederick (1966)',/,
     &        ' 4 : Chaboche (1979)',/,
     &        ' 5 : Yoshida-Uemori',/)
c
      select case ( nkh )
      case ( 0 )
      prop(n0+1  )=float(0)  ! no kinematic hardening
      ndkin=1
      npbs=0
c
      case ( 1 )
      prop(n0+1  )=float(1)  ! Prager(1949)
c                              dx(1)=(2/3)*c*{dpe}
      prop(n0+1+1)=10.0d0    !   constant c
      ndkin=2
      npbs=1
c
      case ( 2 )
      prop(n0+1  )=float(2)  ! Ziegler (1959)
c                              dx(1)=dp*c*{{s}-{X}}
      prop(n0+1+1)=1.0d0     !   constant c
      ndkin=2
      npbs=1
c
      case ( 3 )
      prop(n0+1  )=float(3)  ! Armstrong-Frederick (1966)
c                              dx(1)=(2/3)*c*{dpe}-dp*g*{X}
      prop(n0+1+1)=10.0d0    !   constant c
      prop(n0+1+2)= 0.1d0    !   constant gamma
      ndkin=3
      npbs=1
c
      case ( 4 )
      prop(n0+1  )=float(4)   ! Chaboche (1979)
c                              dx(j)=c(j)*(2/3)*{dpe}-dp*g(j)*{x(j)}
      prop(n0+1+1)= 80000.0d0 ! constant c(1)
      prop(n0+1+2)=   800.0d0 ! constant g(1)
      prop(n0+1+3)=300000.0d0 ! constant c(2)
      prop(n0+1+4)= 10000.0d0 ! constant g(2)
      prop(n0+1+5)=  7500.0d0 ! constant c(3)
      prop(n0+1+6)=     0.0d0 ! constant g(3)
      ndkin=7
      npbs=3
c
      case ( 5 )
c     prop(n0+1  )=float(5)   ! Yoshida-Uemori (*****)
      pc=1.0d0
      py=1.0d0
      pa=1.0d0
      pk=1.0d0
      pb=1.0d0
      prop(n0+1+1)= pc ! C
      prop(n0+1+2)= py ! Y
      prop(n0+1+3)= pa ! a
      prop(n0+1+4)= pk ! k
      prop(n0+1+5)= pb ! b
      ndkin=6
      npbs=2
c
      case default
        write (6,*) 'error retype'
        goto 30
      end select
      n0=n0+ndkin
c
c
      if ( n.gt.nprop ) then
        write (6,*) 'nprop error in jancae_set_prop'
        write (6,*) 'nprop=',nprop
        write (6,*) 'n    =',n
        do i=1,5
          write (6,*) 'prop(',i,')=',prop(i)
        enddo
        call jancae_exit ( 9000 )
      endif
c
      return
      end
c
c
c
c-------------------------------------------------------------
c     set material parameters
c
      subroutine jancae_file_prop ( prop,nprop,io )
c-------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension prop(nprop)
      character flname*32
c
      call jancae_clear1 ( prop,nprop )
c
      flname='prop_ummdp.txt'
      write (6,*) 'read material properties :',flname,io
      open  (io,file=flname,status='old')
      do i=1,5
        read  (io,*) prop(i)
      enddo
      ndela=prop(2)
      ndyld=prop(3)
      ndihd=prop(4)
      ndkin=prop(5)
      npmax=5+ndela+ndyld+ndihd+ndkin
      do i=6,npmax      
        read  (io,*) prop(i)
      enddo
      close (io)
c
      return
      end
c
c
c

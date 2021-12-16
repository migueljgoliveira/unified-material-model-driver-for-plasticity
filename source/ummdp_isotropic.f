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
c
c     ISOTROPIC HARDENING LAW
c
      subroutine ummdp_isotropic ( sy,dsydp,d2sydp2,nreq,p,prihd,ndihd )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: ndihd,nreq
      real*8 ,intent(in) :: p
      real*8 ,intent(in) :: prihd(ndihd)
c
      real*8,intent(out) :: sy,dsydp,d2sydp2
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
        if ( nreq >= 1 ) then
          dsydp = 0.0d0
          if ( nreq >= 2 ) then
            d2sydp2 = 0.0d0
          end if
        end if
c
      case ( 1 )                                                ! Linear
        sy0  = prihd(1+1)
        hard = prihd(1+2)
c
        sy = sy0 + hard*p
        if ( nreq >= 1 ) then
          dsydp = hard
          if ( nreq >= 2 ) then
            d2sydp2 = 0.0d0
          end if
        end if
c
      case ( 2 )                                                 ! Swift
        c  = prihd(1+1)
        e0 = prihd(1+2)
        en = prihd(1+3)
c
        sy = c * (e0+p)**en
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
c
        sy = sy0 + c*p**en
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
c
        sy = sy0 + q*(1.0d0-exp(-b*p))
        if ( nreq >= 1 ) then
          dsydp = q * b  *exp(-b*p)
          if ( nreq >= 2 ) then
            d2sydp2 = -q * b * b * exp(-b*p)
          end if
        end if
c
      case ( 5 )                                         ! Voce & Linear
        sy0 = prihd(1+1)
        q   = prihd(1+2)
        b   = prihd(1+3)
        c   = prihd(1+4)
c
        sy = sy0 + q*(1.0d0-exp(-b*p)) + c*p
        if ( nreq >= 1 ) then
          dsydp = q*b*exp(-b*p) + c
          if ( nreq >= 2 ) then
            d2sydp2 = -q*b*b*exp(-b*p)
          end if
        end if
c
      case ( 6 )                                          ! Voce & Swift
        a   = prihd(1+1)
        sy0 = prihd(1+2)
        q   = prihd(1+3)
        b   = prihd(1+4)
        c   = prihd(1+5)
        e0  = prihd(1+6)
        en  = prihd(1+7)
c
        sy = a*(sy0+q*(1.0d0-exp(-b*p))) + (1.0d0-a)*(c*(e0+p)**en)
        if ( nreq >= 1 ) then
          dsydp = a*(q*b*exp(-b*p)) +(1.0d0-a)*(en*c*(e0+p)**(en-1.0d0))
          if ( nreq >= 2 ) then
            d2sydp2 = a*(-q*b*b*exp(-b*p))
     1                + (1.0d0-a)*(en*c*(en-1.0d0)*(e0+p)**(en-2.0d0))
          end if
        end if
c
      case default
        write (6,*) 'hardening type error',ntihd
        call ummdp_exit ( 203 )
      end select
c
      return
      end subroutine ummdp_isotropic
c
c
c

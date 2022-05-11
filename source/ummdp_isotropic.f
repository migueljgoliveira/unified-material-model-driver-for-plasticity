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
c      5 : Voce & Linear
c      6 : Voce & Swift
c      7 : p-Model
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
      real*8 sy0,hard,c,e0,en,q,b,a,p1,emax
c-----------------------------------------------------------------------
c
      ntihd = nint(prihd(1))
      select case ( ntihd )
c
      case ( 0 )                                     ! Perfectly Plastic
        sy = prihd(1+1)
c
        if ( nreq >= 1 ) then
          dsydp = 0.0d0
        end if
c
        if ( nreq >= 2 ) then
          d2sydp2 = 0.0d0
        end if
c
      case ( 1 )                                                ! Linear
        sy0  = prihd(1+1)
        hard = prihd(1+2)
c
        sy = sy0 + hard*p
c
        if ( nreq >= 1 ) then
          dsydp = hard
        end if
c
        if ( nreq >= 2 ) then
          d2sydp2 = 0.0d0
        end if
c
      case ( 2 )                                                 ! Swift
        c  = prihd(1+1)
        e0 = prihd(1+2)
        en = prihd(1+3)
c
        sy = c * (e0+p)**en
c
        if ( nreq >= 1 ) then
          dsydp = en*c*(e0+p)**(en-1.0d0)
        end if
c
        if ( nreq >= 2 ) then
          d2sydp2 = en*c*(en-1.0d0)*(e0+p)**(en-2.0d0)
        end if
c
      case ( 3 )                                                ! Ludwik
        sy0 = prihd(1+1)
        c   = prihd(1+2)
        en  = prihd(1+3)
c
        sy = sy0 + c*p**en
c
        if ( nreq >= 1 ) then
          p1 = max(p,1.0d-16)
          dsydp = en*c*p1**(en-1.0d0)
        end if
c
        if ( nreq >= 2 ) then
          d2sydp2 = en*c*(en-1.0d0)*p**(en-2.0d0)
        end if
c
      case ( 4 )                                                  ! Voce
        sy0 = prihd(1+1)
        q   = prihd(1+2)
        b   = prihd(1+3)
c
        sy = sy0 + q*(1.0d0-exp(-b*p))
c
        if ( nreq >= 1 ) then
          dsydp = q*b*exp(-b*p)
        end if
c
        if ( nreq >= 2 ) then
          d2sydp2 = -q*b*b*exp(-b*p)
        end if
c
      case ( 5 )                                         ! Voce & Linear
        sy0 = prihd(1+1)
        q   = prihd(1+2)
        b   = prihd(1+3)
        c   = prihd(1+4)
c
        sy = sy0 + q*(1.0d0-exp(-b*p)) + c*p
c
        if ( nreq >= 1 ) then
          dsydp = q*b*exp(-b*p) + c
        end if
c
        if ( nreq >= 2 ) then
          d2sydp2 = -q*b*b*exp(-b*p)
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
c
        if ( nreq >= 1 ) then
          dsydp = a*(q*b*exp(-b*p)) 
     1            + (1.0d0-a)*(en*c*(e0+p)**(en-1.0d0))
        end if
c
        if ( nreq >= 2 ) then
          d2sydp2 = a*(-q*b*b*exp(-b*p))
     1              + (1.0d0-a)*(en*c*(en-1.0d0)*(e0+p)**(en-2.0d0))
        end if
c
      case ( 7 )                                               ! p-Model
        c    = prihd(1+1)
        e0   = prihd(1+2)
        en   = prihd(1+3)
        emax = prihd(1+4)
        b    = prihd(1+5)
c
        if ( p <= emax ) then
          sy = c*(e0+p)**en
        else
          q = (c*en*(e0+emax)**(en-1))/b
          sy = c*(e0+emax)**en + q*(1-exp(-b*(p-emax)))
        end if
c
        if ( nreq >= 1 ) then
          if ( p <= emax ) then
            dsydp = en*c*(e0+p)**(en-1.0d0)
          else
            dsydp = q*b*exp(-b*(p-emax))
          end if
        end if
c
        if ( nreq >= 2 ) then
          if ( p <= emax ) then
            d2sydp2 = en*c*(en-1.0d0)*(e0+p)**(en-2.0d0)
          else
            d2sydp2 = -q*b*b*exp(-b*(p-emax))
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

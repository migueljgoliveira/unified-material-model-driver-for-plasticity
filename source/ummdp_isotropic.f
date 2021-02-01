c***********************************************************************
c     JANCAE/UMMDp : Isotropic Hardening
c***********************************************************************
c
c      0 : Perfectly Plastic
c      1 : Linear
c      2 : Swift
c      3 : Ludwick
c      4 : Voce
c      5 : Voce + Linear
c      6 : Voce + Swift
c
c-----------------------------------------------------------------------
c     hardening curve
c
      subroutine jancae_hardencurve ( sy,dsydp,d2sydp2,
     &                                nreq,p,prihd,ndihd )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension prihd(ndihd)
c
      ntihd = nint(prihd(1))
      select case ( ntihd )
c
      case ( 0 )                                     ! Perfectly Plastic
        sy = prihd(1+1)
        if ( nreq . ge.1 ) then
          dsydp = 0.0
          if ( nreq .ge. 2 ) then
            d2sydp2 = 0.0
          endif
        endif
c
      case ( 1 )                                                ! Linear
        sy0 = prihd(1+1)
        hard = prihd(1+2)
        sy = sy0 + hard*p
        if ( nreq .ge. 1 ) then
          dsydp = hard
          if ( nreq .ge. 2 ) then
            d2sydp2 = 0.0
          endif
        endif
c
      case ( 2 )                                                 ! Swift
        c = prihd(1+1)
        e0 = prihd(1+2)
        en = prihd(1+3)
        sy = c*(e0+p)**en
        if ( nreq .ge. 1 ) then
          dsydp = en*c*(e0+p)**(en-1.0d0)
          if ( nreq .ge. 2 ) then
            d2sydp2 = en*c*(en-1.0d0)*(e0+p)**(en-2.0d0)
          endif
        endif
c
      case ( 3 )                                               ! Ludwick
        sy0 = prihd(1+1)
        c   = prihd(1+2)
        en  = prihd(1+3)
        sy = sy0+c*p**en
        if ( nreq .ge. 1 ) then
          dsydp = en*c*p**(en-1.0d0)
          if ( nreq .ge. 2 ) then
            d2sydp2 = en*c*(en-1.0d0)*p**(en-2.0d0)
          endif
        endif
c
      case ( 4 )                                                  ! Voce
        sy0 = prihd(1+1)
        q = prihd(1+2)
        b = prihd(1+3)
        sy = sy0+q*(1.0d0-exp(-b*p))
        if ( nreq .ge. 1 ) then
          dsydp = q*b*exp(-b*p)
          if ( nreq .ge. 2 ) then
            d2sydp2 = -q*b*b*exp(-b*p)
          endif
        endif
c
      case ( 5 )                                         ! Voce + Linear
        sy0 = prihd(1+1)
        q = prihd(1+2)
        b = prihd(1+3)
        c = prihd(1+4)
        sy = sy0+q*(1.0d0-exp(-b*p))+c*p
        if ( nreq .ge. 1 ) then
          dsydp = q*b*exp(-b*p)+c
          if ( nreq .ge. 2 ) then
            d2sydp2 = -q*b*b*exp(-b*p)
          endif
        endif
c
      case ( 6 )                                          ! Voce + Swift
        a = prihd(1+1)
        sy0 = prihd(1+2)
        q = prihd(1+3)
        b = prihd(1+4)
        c = prihd(1+5)
        e0 = prihd(1+6)
        en = prihd(1+7)
        sy = a*(sy0+q*(1.0d0-exp(-b*p))) + (1.0d0-a)*(c*(e0+p)**en)
        if ( nreq .ge. 1 ) then
          dsydp = a*(q*b*exp(-b*p)) + (1.0d0-a)*(en*c*(e0+p)**(en-1.0d0))
          if ( nreq .ge. 2 ) then
            d2sydp2 = a*(-q*b*b*exp(-b*p)) + 
     &                (1.0d0-a)*(en*c*(en-1.0d0)*(e0+p)**(en-2.0d0))
          endif
        endif
c
      case default
        write (6,*) 'hardening type error',ntihd
        call jancae_exit (9000)
c
      end select
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     print parameters for isotropic hardening info
c
      subroutine jancae_harden_print ( prihd,ndihd )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension prihd(ndihd)
c
      ntihd = nint(prihd(1))
      write (6,*)
      write (6,*) '*** Isotropic Hardening Law',ntihd
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
     1 c*(e0+p)^en)'
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
      end
c
c
c

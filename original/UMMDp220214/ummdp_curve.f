c************************************************************
c     JANCAE/UMMDp : Flow Curve Library 
c************************************************************
c
c
c------------------------------------------------------------
c     hardening curve 
c
      subroutine jancae_hardencurve ( sy,dsydp,d2sydp2,
     &                                nreq,p,prihd,ndihd )
c------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension prihd(ndihd)
c
      ntihd=nint(prihd(1))
      select case ( ntihd )
c                                       ** perfectly plastic
      case ( 0 )
        sy=prihd(1+1)
        if ( nreq.ge.1 ) then
          dsydp=0.0
          if ( nreq.ge.2 ) then
            d2sydp2=0.0
          endif
        endif
c                                        ** linear hardening
      case ( 1 )
        sy0 =prihd(1+1)
        hard=prihd(1+2)
        sy=sy0+hard*p
        if ( nreq.ge.1 ) then
          dsydp=hard
          if ( nreq.ge.2 ) then
            d2sydp2=0.0
          endif
        endif
c                                              ** Swift type
      case ( 2 )
        c =prihd(1+1)
        e0=prihd(1+2)
        en=prihd(1+3)
        sy=c*(e0+p)**en
        if ( nreq.ge.1 ) then
          dsydp=en*c*(e0+p)**(en-1.0d0)
          if ( nreq.ge.2 ) then
            d2sydp2=en*c*(en-1.0d0)*(e0+p)**(en-2.0d0)
          endif
        endif
c                                            ** Ludwick type
      case ( 3 ) 
        sy0=prihd(1+1)
        c  =prihd(1+2)
        en =prihd(1+3)
        sy=sy0+c*p**en
        if ( nreq.ge.1 ) then
          p1=max(p,1.0d-16)            ! 220214 updated
          dsydp=en*c*p1**(en-1.0d0)
          if ( nreq.ge.2 ) then
            d2sydp2=en*c*(en-1.0d0)*p1**(en-2.0d0)
          endif
        endif
c                                               ** Voce type
      case ( 4 )
        sy0=prihd(1+1)
        q  =prihd(1+2)
        b  =prihd(1+3)
        sy=sy0+q*(1.0d0-exp(-b*p))
        if ( nreq.ge.1 ) then
          dsydp=q*b*exp(-b*p)
          if ( nreq.ge.2 ) then
            d2sydp2=-q*b*b*exp(-b*p)
          endif
        endif
c                                       ** Voce + Linear type
      case ( 5 )
        sy0=prihd(1+1)
        q  =prihd(1+2)
        b  =prihd(1+3)
        c  =prihd(1+4)
        sy=sy0+q*(1.0d0-exp(-b*p))+c*p
        if ( nreq.ge.1 ) then
          dsydp=q*b*exp(-b*p)+c
          if ( nreq.ge.2 ) then
            d2sydp2=-q*b*b*exp(-b*p)
          endif
        endif
c                                        ** Voce + Swift type
      case ( 6 )
        a  =prihd(1+1)
        sy0=prihd(1+2)
        q  =prihd(1+3)
        b  =prihd(1+4)
        c  =prihd(1+5)
        e0 =prihd(1+6)
        en =prihd(1+7)
        sy=       a  * ( sy0+q*(1.0d0-exp(-b*p)) )+
     &     (1.0d0-a) * ( c*(e0+p)**en            )
        if ( nreq.ge.1 ) then
          dsydp=       a  * ( q*b*exp(-b*p)           )+
     &          (1.0d0-a) * ( en*c*(e0+p)**(en-1.0d0) )
          if ( nreq.ge.2 ) then
            d2sydp2=
     &           a  * ( -q*b*b*exp(-b*p)                   )+
     &    (1.0d0-a) * ( en*c*(en-1.0d0)*(e0+p)**(en-2.0d0) )
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
c------------------------------------------------------------
c     print parameters for isotropic hardening info
c
      subroutine jancae_harden_print ( prihd,ndihd )
c------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension prihd(ndihd)
c
      ntihd=nint(prihd(1))
      write (6,*) '*** isotropic hardening curve',ntihd
      select case ( ntihd )
      case ( 0 )
        write (6,*) 'Perfect plasticity'
        write (6,*) 'sy_const=',prihd(1+1)
      case ( 1 )
        write (6,*) 'Linear sy0+h*p'
        write (6,*) 'sy0=',prihd(1+1)
        write (6,*) 'h  =',prihd(1+2)
      case ( 2 )
        write (6,*) 'Swift c*(e0+p)^en'
        write (6,*) 'c =',prihd(1+1)
        write (6,*) 'e0=',prihd(1+2)
        write (6,*) 'en=',prihd(1+3)
      case ( 3 )
        write (6,*) 'Ludwick sy0+c*p^en'
        write (6,*) 'sy0=',prihd(1+1)
        write (6,*) 'c  =',prihd(1+2)
        write (6,*) 'en =',prihd(1+3)
      case ( 4 )
        write (6,*) 'Voce sy0+q*(1-exp(-b*p))'
        write (6,*) 'sy0=',prihd(1+1)
        write (6,*) 'q  =',prihd(1+2)
        write (6,*) 'b  =',prihd(1+3)
      case ( 5 )
        write (6,*) 'Voce+Linear sy0+q*(1-exp(-b*p))+c*p'
        write (6,*) 'sy0=',prihd(1+1)
        write (6,*) 'q  =',prihd(1+2)
        write (6,*) 'b  =',prihd(1+3)
        write (6,*) 'c  =',prihd(1+4)
      case ( 6 )
        write (6,*) 'Voce+Swift a *( sy0+q*(1-exp(-b*p)) )+'
        write (6,*) '        (1-a)*( c*(e0+p)^en         )'
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

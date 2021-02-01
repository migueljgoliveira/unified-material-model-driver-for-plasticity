c---------------------------------------------------------------(vegter)
c     Vegter(2006) yield function and its differentials
c     (International Journal of Plasticity 22(2006) P-557-580)
c     
c     This is plane stress yield function.
c
c     Version1.5 (Build20121224)
c        coded by Tomokage Inoue(Aisin AW Industries Co.,Ltd.)
c     Version1.5a 20150225
c        modified by H.Takizawa
c
      subroutine jancae_vegter ( s,se,dseds,d2seds2,nreq,
     &                           pryld,ndyld )
c
      implicit none
      integer, intent(in) :: nreq , ndyld
      real*8,  intent(in) :: s(3), pryld(ndyld)
      real*8,  intent(out) :: se, dseds(3), d2seds2(3,3)
c
      integer  :: nf
c
      nf=nint(pryld(2))-1
      call jancae_vegter_core ( s,se,dseds,d2seds2,nreq,
     &                          pryld,ndyld,nf )
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c
      subroutine jancae_vegter_core ( s,se,dseds,d2seds2,nreq,
     &                                pryld,ndyld,nf )
c-----------------------------------------------------------------------
c
      implicit none
c                                                    ---- input & output  
      integer, intent(in) :: nreq , ndyld , nf
      real*8, intent(in) :: s(3), pryld(ndyld)
      real*8, intent(out) :: se, dseds(3), d2seds2(3,3)
c
c                                                   ---- local variables
      real*8, parameter :: pi=3.141592653589793d0
      real*8, parameter :: TOL0=1.0d-8
      real*8, parameter :: TOL2=1.0d-2
      real*8, parameter :: TOL =1.0d-4
c
c        TOL0    : exception treatment tolerance of stress state 
c        TOL2    : se tolerance change f(1) to f(2)
c        TOL     : exception treatment tolerance of vsin2t
c
c ---------------------------------
      real*8 x(4),vsqrt,vcos2t,vsin2t,theta,theta_rv
      integer i, j, k, m, n,isflag,iareaflag
      real*8 phi_un(0:nf),phi_sh(0:nf),phi_ps(0:nf),f_bi0,
     &       omg(0:nf),r_bi0
      real*8 fun1,fun2,run,fsh1,fsh2,rsh,fps1,fps2,rps,fbi1,fbi2,rbi,
     &       fun1r,fun2r,runr,fps1r,fps2r,rpsr
      real*8 alfa,a(2),b(2),c(2),mm(2),nn(2),mu,f(2)
      real*8 beta,aa,bb,cc,dd
c
c                            over this line for equivalent stress
c                            ---------------------------------------
c
      real*8 dxds(3,3),dphidx(3),dfdmu(2),dadc(2),dcdc(2),dbdc(2)
      real*8 dndc(2),dmdc(2),dmdctmp,dndctmp,P(2),nnmm, dfdc(2)
c  -------------------------------------------
c
c                            over this line for 1st order differential
c                            -----------------------------------------
c
      real*8 d2adc2(2),d2bdc2(2),d2cdc2(2),d2ndc2(2),d2mdc2(2),
     &       vvtmp(0:6),vvtmp_rv(0:6),d2mdc2tmp,d2ndc2tmp,
     &       d2phidx2(3,3)
      integer ithetaflag
c
c                            over this line for 2nd order differential
c-----------------------------------------------------------------------
c                             ---- vegter's parameters from main routine
c     pryld(2) : nf           ! fourier term
c     pryld(3) : f_bi0        ! eq. bi-axial
c     pryld(4) : r_bi0        ! eq. bi-axial
c     i : 1~nf
c     n0=4+i*4
c     test angle  : 90/(nf)*i
c     pryld(n0+1) : phi_un(i) ! uniaxial
c     pryld(n0+2) : phi_sh(i) ! pure shear
c     pryld(n0+3) : phi_ps(i) ! plane strain
c     pryld(n0+4) : omg(   i) 
c
c     this program assumes 2 type of tests
c       type 1 : 3 tests (0,45,90 deg)
c                set :phi_un,phi_sh,phi_ps,omg (1 to 3)
c                     phi_un,phi_sh,phi_ps,omg (4 to 6) =0.0
c       type 2 : 6 tests (0,15,30,45,60,75,90 deg)
c                set :phi_un,phi_sh,phi_ps,omg (1 to 6)
c
c
c
c     2015.2.25 new storage format : H.Takizawa
      f_bi0=     pryld(3)
      r_bi0=     pryld(4)
      do i=0,nf
        phi_un(i)=pryld(4+i*4+1)
        phi_sh(i)=pryld(4+i*4+2)
        phi_ps(i)=pryld(4+i*4+3)
        omg(   i)=pryld(4+i*4+4)
      enddo
c
c
      se=0.0
      do i=1,3
        se=se+s(i)**2
      enddo
      if ( se.le.0.0 ) then
        se=0.0
        return
      endif

c                                       ---- calc x(i) from eq.(14)~(17)
c     x(i)    : Principal stress(i=1^2)
c             : cos2theta(i=3), sin2theta(i=4)
c     isflag  : 0 s(i,i=1,3)=0
c             : 1 not s(i)=0
c                                 ---- exception treatment if all s(i)=0
c
      if(abs(s(1)).le.TOL0.and.abs(s(2)).le.TOL0.and.
     &                                         abs(s(3)).le.TOL0)then
      isflag=0
      call jancae_clear1 ( x,4 )
      goto 100
c
      else
      isflag=1
      call jancae_clear1 ( x,4 )
c
c                           ---- exception treatment if s(1)=s(2),s(3)=0
c
      if(abs(s(1)-s(2)).le.TOL0.and.abs(s(3)).le.TOL0) then
        theta=0.0d0
        theta_rv=0.5d0*pi
        x(1)=0.5d0*(s(1)+s(2))
        x(2)=0.5d0*(s(1)+s(2))
        x(3)=cos(2.0d0*theta)
        x(4)=sin(2.0d0*theta)
        vcos2t=x(3)
        vsin2t=x(4)
c
c                                                    ---- normal process
      else
        vsqrt=sqrt((s(1)-s(2))**2+4.0d0*s(3)**2)
        x(1)=0.5d0*(s(1)+s(2)+vsqrt)
        x(2)=0.5d0*(s(1)+s(2)-vsqrt)
        x(3)=(s(1)-s(2))/vsqrt
        x(4)=2.0d0*s(3)/vsqrt
        vcos2t=x(3)
        vsin2t=x(4)
        theta=0.5d0*acos(vcos2t)
        theta_rv=0.5d0*pi-theta
        endif
      endif
c
c                        ---- calc fk(k=un,sh,ps,bi) , rk(k=un,sh,ps,bi)
c
c                                                          ! un=uniaxial
        fun1=0.0d0
        fun1r=0.0d0
        run=0.0d0
        runr=0.0d0
      do m=0,nf
        fun1=fun1+phi_un(m)*cos(2.0d0*dble(m)*theta)
        fun1r=fun1r+phi_un(m)*cos(2.0d0*dble(m)*theta_rv)
        run=run+omg(m)*cos(2.0d0*dble(m)*theta)
        runr=runr+omg(m)*cos(2.0d0*dble(m)*theta_rv)
      enddo
        fun2=0.0d0
        fun2r=0.0d0
c                                                        ! sh=pure shear
        fsh1=0.0d0
        fsh2=0.0d0
      do m=0,nf
        fsh1=fsh1+phi_sh(m)*cos(2.0d0*dble(m)*theta)
        fsh2=fsh2-phi_sh(m)*cos(2.0d0*dble(m)*theta_rv)
      enddo
        rsh=-1.0d0
c                                                      ! ps=plane strain
        fps1=0.0d0
        fps1r=0.0d0
        rps=0.0d0
        rpsr=0.0d0
      do m=0,nf
        fps1=fps1+phi_ps(m)*cos(2.0d0*dble(m)*theta)
        fps2=0.5d0*fps1
        fps1r=fps1r+phi_ps(m)*cos(2.0d0*dble(m)*theta_rv)
        fps2r=0.5d0*fps1r
      enddo
        rps=-0.0d0
        rpsr=0.0d0
c                                                      ! bi=equi-biaxial
        fbi1=f_bi0
        fbi2=fbi1
        rbi=((r_bi0+1.0d0)+(r_bi0-1.0d0)*vcos2t)/
     &           ((r_bi0+1.0d0)-(r_bi0-1.0d0)*vcos2t)
c
c                                 ---- case distribution by stress state
      if(x(1).ne.0.0d0)then
        alfa=x(2)/x(1)
      endif
      if(x(2).ne.0.0d0)then
        beta=x(1)/x(2)
      endif
c
c     iareaflag    :stress state flag(i=0~6)
c
      if(x(1).gt.0.0d0.and.alfa.lt.0.0d0.and.alfa.ge.fsh2/fsh1) then
        iareaflag=1
      else if(x(1).gt.0.0d0.and.alfa.ge.0.0d0
     &                           .and.alfa.lt.fps2/fps1) then
        iareaflag=2
      else if(x(1).gt.0.0d0.and.alfa.ge.fps2/fps1
     &                           .and.alfa.le.1.0d0) then
        iareaflag=3
c
      else if(x(1).lt.0.0d0.and.alfa.ge.1.0d0
     &                           .and.alfa.lt.fps1r/fps2r) then
        iareaflag=4
      else if(x(1).lt.0.0d0.and.beta.le.fps2r/fps1r
     &                                      .and.beta.gt.0.0d0) then
        iareaflag=5
      else if(x(1).ge.0.0d0.and.beta.le.0.0d0
     &                           .and.beta.gt.fsh1/fsh2) then
        iareaflag=6
c
      else
        go to 100
      endif
c
c                                       ---- calc. hingepoint b(i,i=1~2)
      select case ( iareaflag )
c                                            
      case ( 1 )                                           ! iareaflag=1
         a(1)=fsh1
         a(2)=fsh2
         c(1)=fun1
         c(2)=fun2
         nn(1)=1.0d0
         nn(2)=rsh
         mm(1)=1.0d0
         mm(2)=run
         call jancae_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
c                                            
      case ( 2 )                                           ! iareaflag=2
         a(1)=fun1
         a(2)=fun2
         c(1)=fps1
         c(2)=fps2
         nn(1)=1.0d0
         nn(2)=run
         mm(1)=1.0d0
         mm(2)=rps
         call jancae_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
c                                   
      case ( 3 )                                           ! iareaflag=3
         a(1)=fps1
         a(2)=fps2
         c(1)=fbi1
         c(2)=fbi2
         nn(1)=1.0d0
         nn(2)=rps
         mm(1)=1.0d0
         mm(2)=rbi
         call jancae_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
c                                            
      case ( 4 )                                           ! iareaflag=4
         a(1)=-fbi1
         a(2)=-fbi2
         c(1)=-fps2r
         c(2)=-fps1r
         nn(1)=-1.0d0
         nn(2)=-rbi
         mm(1)=-rpsr
         mm(2)=-1.0d0
         call jancae_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
c                                            
      case ( 5 )                                           ! iareaflag=5
         a(1)=-fps2r
         a(2)=-fps1r
         c(1)=-fun2r
         c(2)=-fun1r
         nn(1)=-rpsr
         nn(2)=-1.0d0
         mm(1)=-runr
         mm(2)=-1.0d0
         call jancae_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
c                                    
      case ( 6 )                                           ! iareaflag=6
         a(1)=-fun2r
         a(2)=-fun1r
         c(1)=fsh1
         c(2)=fsh2
         nn(1)=-runr
         nn(2)=-1.0d0
         mm(1)=1.0d0
         mm(2)=rsh
         call jancae_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
c
      case default
         write (6,*) 'iareaflag error :',iareaflag
         call jancae_exit (9000)
      end select

c                            ---- calc. fourier coefficient mu(0<=mu<=1)
      call jancae_vegter_calc_mu ( x,a,b,c,mu,iareaflag,s,theta
     &                                            ,aa,bb,cc,dd )
c
c                            ---- calc. normalized yield locus f(i)i=1~2
      call jancae_vegter_calc_fi ( x,a,b,c,mu,f )
      go to 200
c
c                                                 ---- equivalent stress
  100 continue
      se=0.0d0
      go to 300
  200 continue
      if(f(1).le.TOL2) then
         se=x(2)/f(2)
      else
         se=x(1)/f(1)
      endif
c
      go to 300
c
  300 continue
c
c                                            ---- 1st order differential
c
      if ( nreq.ge.1 ) then
c                        ---- set dadc,dcdc,dndc,dmdc for eq.(A.7)^(A.9)
c
      call jancae_clear1 ( dadc,2 )
      call jancae_clear1 ( dbdc,2 )
      call jancae_clear1 ( dcdc,2 )
      call jancae_clear1 ( dndc,2 )
      call jancae_clear1 ( dmdc,2 )
c
      select case ( iareaflag )
c                                            
      case ( 1 )                                           ! iareaflag=1
       if(abs(vsin2t).ge.TOL) then
        do m=0,nf
          dadc(1)=dadc(1)+phi_sh(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta)
     &                             /sin(2.0d0*theta)
          dadc(2)=dadc(2)+phi_sh(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta_rv)
     &                             /sin(2.0d0*theta)
        enddo
       else
          do m=0,nf
              dadc(1)=dadc(1)+phi_sh(m)*dble(m)**2
              dadc(2)=dadc(2)+phi_sh(m)*dble(m)**2
          enddo
       endif
c
       if(abs(vsin2t).ge.TOL) then
         do m=0,nf
          dcdc(1)=dcdc(1)+phi_un(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta)
     &                             /sin(2.0d0*theta)
         enddo
        else
          do m=0,nf
              dcdc(1)=dcdc(1)+phi_un(m)*dble(m)**2
          enddo
        endif
          dcdc(2)=0.0d0
c
          dndc(1)=0.0d0
          dndc(2)=0.0d0
c
          dmdc(1)=0.0d0
       if(abs(vsin2t).ge.TOL) then
         do m=0,nf
          dmdc(2)=dmdc(2)+omg(m)*dble(m)
     &                          *sin(2.0d0*dble(m)*theta)
     &                          /sin(2.0d0*theta)
         enddo
        else
          do m=0,nf
              dmdc(2)=dmdc(2)+omg(m)*dble(m)**2
          enddo
        endif
      call jancae_vegter_calc_dbdc ( a,b,c,dadc,dbdc,dcdc,mm,nn,
     &                               dndc,dmdc,iareaflag,P,nnmm)
c
c                                            
      case ( 2 )                                           ! iareaflag=2
       if(abs(vsin2t).ge.TOL) then
         do m=0,nf
          dadc(1)=dadc(1)+phi_un(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta)
     &                             /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dadc(1)=dadc(1)+phi_un(m)*dble(m)**2
         enddo
       endif
          dadc(2)=0.0d0
c
       if(abs(vsin2t).ge.TOL) then
         do m=0,nf
          dcdc(1)=dcdc(1)+phi_ps(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta)
     &                             /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dcdc(1)=dcdc(1)+phi_ps(m)*dble(m)**2
         enddo
       endif
          dcdc(2)=0.5d0*dcdc(1)
c
          dndc(1)=0.0d0
       if(abs(vsin2t).ge.TOL) then
        do m=0,nf
          dndc(2)=dndc(2)+omg(m)*dble(m)
     &                           *sin(2.0d0*dble(m)*theta)
     &                           /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dndc(2)=dndc(2)+omg(m)*dble(m)**2
         enddo
       endif
c
          dmdc(1)=0.0d0
          dmdc(2)=0.0d0
      call jancae_vegter_calc_dbdc ( a,b,c,dadc,dbdc,dcdc,mm,nn,
     &                               dndc,dmdc,iareaflag,P,nnmm)
c
c                                           
      case ( 3 )                                           ! iareaflag=3
       if(abs(vsin2t).ge.TOL) then
        do m=0,nf
          dadc(1)=dadc(1)+phi_ps(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta)
     &                             /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dadc(1)=dadc(1)+phi_ps(m)*dble(m)**2
         enddo
       endif
          dadc(2)=0.5d0*dadc(1)
c
          dcdc(1)=0.0d0
          dcdc(2)=0.0d0
c
          dndc(1)=0.0d0
          dndc(2)=0.0d0
c
          dmdc(1)=0.0d0
          dmdctmp=r_bi0+1.0d0-(r_bi0-1.0d0)*vcos2t
          dmdc(2)=2.0d0*(r_bi0*r_bi0-1.0d0)/(dmdctmp*dmdctmp)
      call jancae_vegter_calc_dbdc ( a,b,c,dadc,dbdc,dcdc,mm,nn,
     &                               dndc,dmdc,iareaflag,P,nnmm)
c
c                                            
      case ( 4 )                                           ! iareaflag=4
          dadc(1)=0.0d0
          dadc(2)=0.0d0
c
       if(abs(vsin2t).ge.TOL) then
        do m=0,nf
          dcdc(2)=dcdc(2)+phi_ps(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta_rv)
     &                             /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dcdc(2)=dcdc(2)+phi_ps(m)*dble(m)**2
         enddo
       endif
          dcdc(1)=0.5d0*dcdc(2)
c
          dndc(1)=0.0d0
          dndctmp=r_bi0+1.0d0-(r_bi0-1.0d0)*vcos2t
          dndc(2)=-2.0d0*(r_bi0*r_bi0-1.0d0)/(dndctmp*dndctmp)
c
          dmdc(1)=0.0d0
          dmdc(2)=0.0d0
      call jancae_vegter_calc_dbdc ( a,b,c,dadc,dbdc,dcdc,mm,nn,
     &                               dndc,dmdc,iareaflag,P,nnmm)
c
c                                            
      case ( 5 )                                           ! iareaflag=5
       if(abs(vsin2t).ge.TOL) then
        do m=0,nf
          dadc(2)=dadc(2)+phi_ps(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta_rv)
     &                             /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dadc(2)=dadc(2)+phi_ps(m)*dble(m)**2
         enddo
       endif
          dadc(1)=0.5d0*dadc(2)
c
          dcdc(1)=0.0d0
       if(abs(vsin2t).ge.TOL) then
        do m=0,nf
          dcdc(2)=dcdc(2)+phi_un(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta_rv)
     &                             /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dcdc(2)=dcdc(2)+phi_un(m)*dble(m)**2
         enddo
       endif
c
          dndc(1)=0.0d0
          dndc(2)=0.0d0
c
       if(abs(vsin2t).ge.TOL) then
        do m=0,nf
          dmdc(1)=dmdc(1)+omg(m)*dble(m)
     &                          *sin(2.0d0*dble(m)*theta_rv)
     &                          /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dmdc(1)=dmdc(1)+omg(m)*dble(m)**2
         enddo
       endif
          dmdc(2)=0.0d0
      call jancae_vegter_calc_dbdc ( a,b,c,dadc,dbdc,dcdc,mm,nn,
     &                               dndc,dmdc,iareaflag,P,nnmm)
c
c                                            
      case ( 6 )                                           ! iareaflag=6
          dadc(1)=0.0d0
       if(abs(vsin2t).ge.TOL) then
        do m=0,nf
          dadc(2)=dadc(2)+phi_un(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta_rv)
     &                             /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dadc(2)=dadc(2)+phi_un(m)*dble(m)**2
         enddo
       endif
c
       if(abs(vsin2t).ge.TOL) then
        do m=0,nf
          dcdc(1)=dcdc(1)+phi_sh(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta)
     &                             /sin(2.0d0*theta)
          dcdc(2)=dcdc(2)+phi_sh(m)*dble(m)
     &                             *sin(2.0d0*dble(m)*theta_rv)
     &                             /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dcdc(1)=dcdc(1)+phi_sh(m)*dble(m)**2
          dcdc(2)=dcdc(2)+phi_sh(m)*dble(m)**2
         enddo
       endif
c
       if(abs(vsin2t).ge.TOL) then
        do m=0,nf
          dndc(1)=dndc(1)+omg(m)*dble(m)
     &                          *sin(2.0d0*dble(m)*theta_rv)
     &                          /sin(2.0d0*theta)
         enddo
        else
         do m=0,nf
          dndc(1)=dndc(1)+omg(m)*dble(m)**2
         enddo
       endif
          dndc(2)=0.0d0
c
          dmdc(1)=0.0d0
          dmdc(2)=0.0d0
      call jancae_vegter_calc_dbdc ( a,b,c,dadc,dbdc,dcdc,mm,nn,
     &                               dndc,dmdc,iareaflag,P,nnmm)
c
      case default
        write (6,*) 'iareaflag error(dseds) :',iareaflag
        call jancae_exit (9000)
      end select
c
c
c                             ---- calc. dphidx(i) (i=1~3)  eq.(21)~(23)
      call jancae_vegter_calc_dphidx ( dphidx,se,a,b,c,dadc,dbdc,dcdc,
     &                      mm,nn,dndc,dmdc,f,iareaflag,mu,x,dfdmu,dfdc)
c
c
c                                   ---- calc. dseds(i) (i=1~3)  eq.(20)
      call jancae_vegter_calc_dseds ( dseds,x,dphidx,vcos2t,vsin2t,
     &                                iareaflag,isflag,dxds)
c
      endif
c
c
c                                            ---- 2nd order differential
      if ( nreq.ge.2 ) then
c                     ---- set d2adc2,d2cdc2,d2ndc2,d2mdc2 for d2bdc2(2)
c
      call jancae_clear1 ( d2adc2,2 )
      call jancae_clear1 ( d2bdc2,2 )
      call jancae_clear1 ( d2cdc2,2 )
      call jancae_clear1 ( d2ndc2,2 )
      call jancae_clear1 ( d2mdc2,2 )
c
c                     ---- define exception treatment condition of theta
c                        ---- if theta<=0.002865deg then apply exception
c
      if(abs(vsin2t).le.TOL) then
         ithetaflag=1
      else 
         ithetaflag=0
      endif
c
      if(ithetaflag.eq.1) then
        vvtmp(0)=0.0d0
        vvtmp(1)=0.0d0
        vvtmp(2)=2.0d0
        vvtmp(3)=8.0d0
        vvtmp(4)=20.0d0
        vvtmp(5)=40.0d0
        vvtmp(6)=70.0d0
c
        vvtmp_rv(0)=0.0d0
        vvtmp_rv(1)=0.0d0
        vvtmp_rv(2)=-2.0d0
        vvtmp_rv(3)=8.0d0
        vvtmp_rv(4)=-20.0d0
        vvtmp_rv(5)=40.0d0
        vvtmp_rv(6)=-70.0d0
c
      else
       do m=0,nf
        vvtmp(m)=cos(2.0d0*theta)*sin(2.0d0*dble(m)*theta)/
     &          (sin(2.0d0*theta)**3)-dble(m)*cos(2.0d0*dble(m)*theta)/
     &                                           (sin(2.0d0*theta)**2)
        vvtmp_rv(m)=cos(2.0d0*theta)*sin(2.0d0*dble(m)*theta_rv)/
     &       (sin(2.0d0*theta)**3)+dble(m)*cos(2.0d0*dble(m)*theta_rv)/
     &                                           (sin(2.0d0*theta)**2)
       enddo
c
      endif
c
      select case ( iareaflag )
c                                            
      case ( 1 )                                           ! iareaflag=1
         do m=0,nf
          d2adc2(1)=d2adc2(1)+phi_sh(m)*dble(m)*vvtmp(m)
          d2adc2(2)=d2adc2(2)+phi_sh(m)*dble(m)*vvtmp_rv(m)
         enddo
c
         do m=0,nf
          d2cdc2(1)=d2cdc2(1)+phi_un(m)*dble(m)*vvtmp(m)
         enddo
          d2cdc2(2)=0.0d0
c
          d2ndc2(1)=0.0d0
          d2ndc2(2)=0.0d0
c
          d2mdc2(1)=0.0d0
         do m=0,nf
          d2mdc2(2)=d2mdc2(2)+omg(m)*dble(m)*vvtmp(m)
         enddo
      call jancae_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,nn,dndc,
     &     dmdc,iareaflag,d2adc2,d2bdc2,d2cdc2,d2ndc2,d2mdc2,P,nnmm,s)
c
c                                            
      case ( 2 )                                           ! iareaflag=2
         do m=0,nf
          d2adc2(1)=d2adc2(1)+phi_un(m)*dble(m)*vvtmp(m)
         enddo
          d2adc2(2)=0.0d0
c
         do m=0,nf
          d2cdc2(1)=d2cdc2(1)+phi_ps(m)*dble(m)*vvtmp(m)
         enddo
          d2cdc2(2)=0.5d0*d2cdc2(1)
c
          d2ndc2(1)=0.0d0
         do m=0,nf
          d2ndc2(2)=d2ndc2(2)+omg(m)*dble(m)*vvtmp(m)
         enddo
c
          d2mdc2(1)=0.0d0
          d2mdc2(2)=0.0d0
      call jancae_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,nn,dndc,
     &     dmdc,iareaflag,d2adc2,d2bdc2,d2cdc2,d2ndc2,d2mdc2,P,nnmm,s)
c
c                                            
      case ( 3 )                                           ! iareaflag=3
         do m=0,nf  
          d2adc2(1)=d2adc2(1)+phi_ps(m)*dble(m)*vvtmp(m)
         enddo
          d2adc2(2)=0.5d0*d2adc2(1)
c
          d2cdc2(1)=0.0d0
          d2cdc2(2)=0.0d0
c
          d2ndc2(1)=0.0d0
          d2ndc2(2)=0.0d0
c
          d2mdc2(1)=0.0d0
             d2mdc2tmp=r_bi0+1.0d0-(r_bi0-1.0d0)*vcos2t
          d2mdc2(2)=4.0d0*(r_bi0**2-1.0d0)*(r_bi0-1.0d0)/
     &                                              (d2mdc2tmp**3)
      call jancae_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,nn,dndc,
     &     dmdc,iareaflag,d2adc2,d2bdc2,d2cdc2,d2ndc2,d2mdc2,P,nnmm,s)
c
c                                            
      case ( 4 )                                          !  iareaflag=4
          d2adc2(1)=0.0d0
          d2adc2(2)=0.0d0
c
         do m=0,nf
          d2cdc2(2)=d2cdc2(2)+phi_ps(m)*dble(m)*vvtmp_rv(m)
         enddo
          d2cdc2(1)=0.5d0*d2cdc2(2)
c
          d2ndc2(1)=0.0d0
             d2ndc2tmp=r_bi0+1.0d0-(r_bi0-1.0d0)*vcos2t
          d2ndc2(2)=-4.0d0*(r_bi0**2-1.0d0)*(r_bi0-1.0d0)/
     &                                              (d2ndc2tmp**3)
c
          d2mdc2(1)=0.0d0
          d2mdc2(2)=0.0d0
      call jancae_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,nn,dndc,
     &     dmdc,iareaflag,d2adc2,d2bdc2,d2cdc2,d2ndc2,d2mdc2,P,nnmm,s)
c
c                                           
      case ( 5 )                                           ! iareaflag=5
         do m=0,nf
          d2adc2(2)=d2adc2(2)+phi_ps(m)*dble(m)*vvtmp_rv(m)
         enddo
          d2adc2(1)=0.5d0*d2adc2(2)
c
          d2cdc2(1)=0.0d0
         do m=0,nf
          d2cdc2(2)=d2cdc2(2)+phi_un(m)*dble(m)*vvtmp_rv(m)
         enddo
c
          d2ndc2(1)=0.0d0
          d2ndc2(2)=0.0d0
c
         do m=0,nf
          d2mdc2(1)=d2mdc2(1)+omg(m)*dble(m)*vvtmp_rv(m)
         enddo
          d2mdc2(2)=0.0d0
      call jancae_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,nn,dndc,
     &     dmdc,iareaflag,d2adc2,d2bdc2,d2cdc2,d2ndc2,d2mdc2,P,nnmm,s)
c
c                                            
      case ( 6 )                                           ! iareaflag=6
          d2adc2(1)=0.0d0
         do m=0,nf
          d2adc2(2)=d2adc2(2)+phi_un(m)*dble(m)*vvtmp_rv(m)
         enddo
c
         do m=0,nf
          d2cdc2(1)=d2cdc2(1)+phi_sh(m)*dble(m)*vvtmp(m)
          d2cdc2(2)=d2cdc2(2)+phi_sh(m)*dble(m)*vvtmp_rv(m)
         enddo
cn
         do m=0,nf
          d2ndc2(1)=d2ndc2(1)+omg(m)*dble(m)*vvtmp_rv(m)
         enddo
          d2ndc2(2)=0.0d0
c
          d2mdc2(1)=0.0d0
          d2mdc2(2)=0.0d0
      call jancae_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,nn,dndc,
     &     dmdc,iareaflag,d2adc2,d2bdc2,d2cdc2,d2ndc2,d2mdc2,P,nnmm,s)
c
      case default
        write (6,*) 'iareaflag error(d2seds2) :',iareaflag
        call jancae_exit (9000)
      end select
c
c
c                                     ---- calc. d2phidx2(k,l) (k,l=1~3)
      call jancae_vegter_calc_d2phidx2 (d2phidx2,se,a,b,c,dadc,dbdc,
     &             dcdc,mm,nn,dndc,dmdc,f,iareaflag,mu,x,d2adc2,d2bdc2,
     &            d2cdc2,d2ndc2,d2mdc2,dfdmu,dfdc,s,aa,bb,cc,dd,dphidx)
c
c
c                                      ---- calc. d2seds2(i,j) (i,j=1~3)
      call jancae_vegter_calc_d2seds2 (d2seds2,d2phidx2,se,a,b,c,mu,x,
     &         vcos2t,vsin2t,iareaflag,dxds,dphidx,isflag,s,dseds,
     &                                 pryld,ndyld)
c
      endif
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c                                      under this line is branch_routine
c---------------------------------------------------------------(vegter)
c     calc. hingepoint b(i,i=1~2)
c
      subroutine jancae_vegter_hingepoint ( a,b,c,mm,nn,iareaflag,s )
c-----------------------------------------------------------------------
      implicit none
c
      real*8, intent(in) :: a(2),c(2),mm(2),nn(2),s(3)
      real*8, intent(out) :: b(2)
      real*8  bb, b1u,b2u
      real*8, parameter :: TOL=1.0d-8
      integer, intent(in) ::iareaflag
c
      b1u=mm(2)*(nn(1)*a(1)+nn(2)*a(2))-nn(2)*(mm(1)*c(1)+mm(2)*c(2))
      b2u=nn(1)*(mm(1)*c(1)+mm(2)*c(2))-mm(1)*(nn(1)*a(1)+nn(2)*a(2))
c
      bb=nn(1)*mm(2)-mm(1)*nn(2)
      if(abs(bb).le.TOL) then
         write (6,*) 'hingepoint singular error! '
         call jancae_exit (9000)
      endif
c
      b(1)=b1u/bb
      b(2)=b2u/bb
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     calc. fourier coefficient mu(0<=mu<=1)
c
      subroutine jancae_vegter_calc_mu ( x,a,b,c,mu,iareaflag,s,theta,
     &                                                     aa,bb,cc,dd)
c-----------------------------------------------------------------------
      implicit none
c
      real*8, intent(in) :: x(4),a(2),b(2),c(2),s(3),theta 
      real*8, intent(out) :: mu
      real*8  aa,bb,cc,dd,xx(2)
      real*8, parameter :: TOL1=1.0d-8
      real*8, parameter :: TOL2=1.0d-8
      integer, intent(in) ::iareaflag
      integer imuflag
c
      aa=x(2)*(a(1)+c(1)-2.0d0*b(1))-x(1)*(a(2)+c(2)-2.0d0*b(2))
      bb=2.0d0*x(2)*(b(1)-a(1))-2.0d0*x(1)*(b(2)-a(2))
      cc=x(2)*a(1)-x(1)*a(2)
c
      if(abs(aa).le.TOL1) then
         write (6,*) 'calc. mu singular error! ',abs(aa),iareaflag
         call jancae_exit (9000)
      endif
c
      dd=bb*bb-4.0d0*aa*cc
      if(dd.ge.0.0d0) then
        xx(1)=0.5d0*(-bb+sign(sqrt(dd),-bb))/aa
        xx(2)=cc/(aa*xx(1))
c
      else
         write (6,*) 'negative dd ! ',dd,iareaflag
         call jancae_exit (9000)
      endif
c
      if(xx(1).ge.0.0d0.and.xx(1).le.1.0000005d0) then
          mu=xx(1)
           imuflag=1
        else if(xx(2).ge.0.0d0.and.xx(2).le.1.0000005d0) then
          mu=xx(2)
           imuflag=2
        else if(abs(xx(1)).le.TOL2.or.abs(xx(2)).le.TOL2)then
          mu=0.0d0
        else
         write (6,*) 'can not find mu ! solve error ',iareaflag,xx(1)
     & ,xx(2)
         call jancae_exit (9000)
      endif
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     calc. normalized yield locus f(i)i=1~2
c
      subroutine jancae_vegter_calc_fi ( x,a,b,c,mu,f )
c-----------------------------------------------------------------------
      implicit none
c
      real*8, intent(in) :: x(4),a(2),b(2),c(2),mu
      real*8, intent(out) :: f(2)
      integer i
c
      do i=1,2
        f(i)=a(i)+2.0d0*mu*(b(i)-a(i))+mu*mu*(a(i)+c(i)-2.0d0*b(i))
      enddo
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     calc. dbdc(i) (i=1~2) eq.(A.7)
c
      subroutine jancae_vegter_calc_dbdc ( a,b,c,dadc,dbdc,
     &                 dcdc, mm,nn,dndc,dmdc,iareaflag,P,nnmm)
c-----------------------------------------------------------------------
      implicit none
c
      real*8, intent(in) :: a(2),b(2),c(2),dadc(2),dcdc(2),mm(2),
     &                      nn(2),dndc(2),dmdc(2)
      real*8, intent(out) :: dbdc(2),P(2),nnmm
      integer, intent(in) :: iareaflag
      real*8, parameter :: TOL=1.0d-8
      real*8 nminv
      integer i
c
      P(1)=nn(1)*dadc(1)+dndc(1)*(a(1)-b(1))+nn(2)*dadc(2)+
     &                                        dndc(2)*(a(2)-b(2))
c
      P(2)=mm(1)*dcdc(1)+dmdc(1)*(c(1)-b(1))+mm(2)*dcdc(2)+
     &                                        dmdc(2)*(c(2)-b(2))
c
      nnmm=nn(1)*mm(2)-mm(1)*nn(2)
         if(abs(nnmm).lt.TOL) then
            write (6,*) 'nnmm too small! ',nnmm
            call jancae_exit (9000)
         endif
      nminv=1.0d0/nnmm
c
      dbdc(1)=nminv*(P(1)*mm(2)-P(2)*nn(2))
      dbdc(2)=nminv*(P(2)*nn(1)-P(1)*mm(1))
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     calc. dphidx(i) (i=1~3)  eq.(21)~(23)
c
      subroutine jancae_vegter_calc_dphidx ( dphidx,se,a,b,c,dadc,
     &            dbdc,dcdc,mm,nn,dndc,dmdc,f,iareaflag,mu,x,dfdmu,dfdc)
c-----------------------------------------------------------------------
      implicit none
c
      real*8, intent(in) :: se,a(2),b(2),c(2),dadc(2),dbdc(2),dcdc(2),
     &                      mm(2),nn(2),dndc(2),dmdc(2),f(2),mu,x(4)
      real*8, intent(out) :: dphidx(3),dfdmu(2),dfdc(2)
      integer, intent(in) :: iareaflag
      real*8, parameter :: TOL1=0.996515d0       ! =44.9deg
      real*8, parameter :: TOL2=1.003497d0       ! =45.1deg
      real*8, parameter :: TOL3=1.0d-8
      real*8 dphidxtmp(3),dphidxcoe,dphidxcinv,
     &       tmp_u,tmp_b,vtan
      integer i
c
c                                    ---- calc. dfdc(i) (i=1~2)  eq.(23)
      do i=1,2
        dfdc(i)=dadc(i)+2.0d0*mu*(dbdc(i)-dadc(i))+mu*mu*(dadc(i)
     &                                       +dcdc(i)-2.0d0*dbdc(i))
      enddo
c
c                                   ---- calc. dfdmu(i) (i=1~2)  eq.(22)
      do i=1,2
        dfdmu(i)=2.0d0*(b(i)-a(i))+2.0d0*mu*(a(i)+c(i)-2.0d0*b(i))
      enddo
c
c                            ---- calc. dphidx(i) (i=1~3)  eq.(21),(C.1)
      dphidxcinv=f(1)*dfdmu(2)-f(2)*dfdmu(1)
         if(abs(dphidxcinv).lt.TOL3) then
            write (6,*) 'eq.(21) too small! ',dphidxcinv
            call jancae_exit (9000)
         endif
             dphidxcoe=1.0d0/dphidxcinv
c
c                            ---- if condition to avoid singular eq.(20)
c                                            ---- apply 44.9 to 45.1 deg
c
      if(iareaflag.eq.3.or.iareaflag.eq.4) then
      vtan=x(2)/x(1)
      end if
c
      if(iareaflag.eq.4.and.vtan.ge.TOL1.and.vtan.le.TOL2) then
c
      tmp_u=1.0d0*(2.0d0*(1.0d0-mu)*dbdc(2)+mu*dcdc(2))*dfdmu(1)
     &     -1.0d0*(2.0d0*(1.0d0-mu)*dbdc(1)+mu*dcdc(1))*dfdmu(2)
c
      tmp_b=2.0d0*(1.0d0-mu)*(b(1)-b(2))+mu*(c(1)-c(2))

      dphidxtmp(1)=dfdmu(2)
      dphidxtmp(2)=-dfdmu(1)
      dphidxtmp(3)=tmp_u/tmp_b
c
      else if(iareaflag.eq.3.and.vtan.ge.TOL1.and.vtan.le.TOL2) then
c
      tmp_u=1.0d0*(2.0d0*mu*dbdc(2)+(1.0d0-mu)*dadc(2))*dfdmu(1)
     &     -1.0d0*(2.0d0*mu*dbdc(1)+(1.0d0-mu)*dadc(1))*dfdmu(2)
c
      tmp_b=2.0d0*mu*(b(1)-b(2))+(1.0d0-mu)*(a(1)-a(2))

      dphidxtmp(1)=dfdmu(2)
      dphidxtmp(2)=-dfdmu(1)
      dphidxtmp(3)=tmp_u/tmp_b
c
      else
c
      dphidxtmp(1)=dfdmu(2)
      dphidxtmp(2)=-dfdmu(1)
      dphidxtmp(3)=se*(dfdc(2)*dfdmu(1)-dfdc(1)*dfdmu(2))
c
      endif
      do i=1,3
        dphidx(i)=dphidxcoe*dphidxtmp(i)
      enddo
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     calc. dseds(i) (i=1~3)  eq.(20)
c
      subroutine jancae_vegter_calc_dseds ( dseds,x,dphidx,vcos2t,
     &                                     vsin2t,iareaflag,isflag,dxds)
c-----------------------------------------------------------------------
      implicit none
c
      real*8, intent(in) :: x(4),dphidx(3),vcos2t,vsin2t
      real*8, intent(out) :: dseds(3),dxds(3,3)
      integer, intent(in) :: iareaflag,isflag
      real*8, parameter :: TOL1=0.996515d0       ! =44.9deg
      real*8, parameter :: TOL2=1.003497d0       ! =45.1deg
      real*8 dxds_t(3,3),vtan
      integer i,j
c
c                  ---- set linear transformation matrix dxds(3,3) eq.18
c
c                            ---- if condition to avoid singular eq.(20)
c                                            ---- apply 44.9 to 45.1 deg
c
      if(iareaflag.eq.3.or.iareaflag.eq.4) then
      vtan=x(2)/x(1)
      end if
c
      if(iareaflag.eq.3.and.vtan.ge.TOL1.and.vtan.le.TOL2) then
        dxds(1,1)=0.5d0*(1.0d0+vcos2t)
        dxds(2,1)=0.5d0*(1.0d0-vcos2t)
        dxds(3,1)=vsin2t*vsin2t
        dxds(1,2)=0.5d0*(1.0d0-vcos2t)
        dxds(2,2)=0.5d0*(1.0d0+vcos2t)
        dxds(3,2)=-vsin2t*vsin2t
        dxds(1,3)=vsin2t
        dxds(2,3)=-vsin2t
        dxds(3,3)=-2.0d0*vsin2t*vcos2t
        dxds_t=transpose(dxds)
c
      else if(iareaflag.eq.4.and.vtan.ge.TOL1.and.vtan.le.TOL2) then
        dxds(1,1)=0.5d0*(1.0d0+vcos2t)
        dxds(2,1)=0.5d0*(1.0d0-vcos2t)
        dxds(3,1)=vsin2t*vsin2t
        dxds(1,2)=0.5d0*(1.0d0-vcos2t)
        dxds(2,2)=0.5d0*(1.0d0+vcos2t)
        dxds(3,2)=-vsin2t*vsin2t
        dxds(1,3)=vsin2t
        dxds(2,3)=-vsin2t
        dxds(3,3)=-2.0d0*vsin2t*vcos2t
        dxds_t=transpose(dxds)
c
      else
        dxds(1,1)=0.5d0*(1.0d0+vcos2t)
        dxds(2,1)=0.5d0*(1.0d0-vcos2t)
        dxds(3,1)=vsin2t*vsin2t/(x(1)-x(2))
        dxds(1,2)=0.5d0*(1.0d0-vcos2t)
        dxds(2,2)=0.5d0*(1.0d0+vcos2t)
        dxds(3,2)=-vsin2t*vsin2t/(x(1)-x(2))
        dxds(1,3)=vsin2t
        dxds(2,3)=-vsin2t
        dxds(3,3)=-2.0d0*vsin2t*vcos2t/(x(1)-x(2))
        dxds_t=transpose(dxds)
      endif
c
c                                   ---- calc. dseds(i) (1=1~3)  eq.(20)
      call jancae_clear1 ( dseds,3 )
c
      do i=1,3
        do j=1,3
          dseds(i)=dseds(i)+dxds_t(i,j)*dphidx(j)
        enddo
      enddo
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     calc. d2bdc2(i) (i=1~2)
c
      subroutine jancae_vegter_calc_d2bdc2 ( a,b,c,dadc,dbdc,dcdc,mm,
     $      nn,dndc,dmdc,iareaflag,d2adc2,d2bdc2,d2cdc2,d2ndc2,d2mdc2,
     &      P,nnmm,s)
c-----------------------------------------------------------------------
      implicit none
c
      real*8, intent(in) :: a(2),b(2),c(2),dadc(2),dcdc(2),mm(2),
     &                      nn(2),dndc(2),dmdc(2),dbdc(2),
     &                      d2adc2(2),d2cdc2(2),d2ndc2(2),d2mdc2(2),
     &                      P(2),nnmm,s(3)
      real*8, intent(out) :: d2bdc2(2)
      integer, intent(in) :: iareaflag
      real*8  dnnmmdc,dPdc(2),dp1m2p2n2,dp2n1p1m1
      integer i
c
      dnnmmdc=dndc(1)*mm(2)+nn(1)*dmdc(2)-dmdc(1)*nn(2)-mm(1)*dndc(2)
c
      dPdc(1)=dndc(1)*dadc(1)+nn(1)*d2adc2(1)+d2ndc2(1)*(a(1)-b(1))
     &                                   +dndc(1)*(dadc(1)-dbdc(1))
     &       +dndc(2)*dadc(2)+nn(2)*d2adc2(2)+d2ndc2(2)*(a(2)-b(2))
     &                                   +dndc(2)*(dadc(2)-dbdc(2))
c
      dPdc(2)=dmdc(1)*dcdc(1)+mm(1)*d2cdc2(1)+d2mdc2(1)*(c(1)-b(1))
     &                                   +dmdc(1)*(dcdc(1)-dbdc(1))
     &       +dmdc(2)*dcdc(2)+mm(2)*d2cdc2(2)+d2mdc2(2)*(c(2)-b(2))
     &                                   +dmdc(2)*(dcdc(2)-dbdc(2))
c
      dp1m2p2n2=dPdc(1)*mm(2)+P(1)*dmdc(2)-dPdc(2)*nn(2)-P(2)*dndc(2)
      dp2n1p1m1=dPdc(2)*nn(1)+P(2)*dndc(1)-dPdc(1)*mm(1)-P(1)*dmdc(1)
c
      d2bdc2(1)=-1.0d0*dnnmmdc*(P(1)*mm(2)-P(2)*nn(2))/(nnmm*nnmm)
     &         +dp1m2p2n2/nnmm
      d2bdc2(2)=-1.0d0*dnnmmdc*(P(2)*nn(1)-P(1)*mm(1))/(nnmm*nnmm)
     &         +dp2n1p1m1/nnmm
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     calc. d2phidx2(k,l) (k,l=1~3)
c
      subroutine jancae_vegter_calc_d2phidx2 (d2phidx2,se,a,b,c,dadc,
     &              dbdc,dcdc,mm,nn,dndc,dmdc,f,iareaflag,mu,x,d2adc2,
     &      d2bdc2,d2cdc2,d2ndc2,d2mdc2,dfdmu,dfdc,s,aa,bb,cc,dd,dphidx)
c-----------------------------------------------------------------------
      implicit none
c
      real*8, intent(in) :: se,a(2),b(2),c(2),dadc(2),dbdc(2),dcdc(2),
     &                      mm(2),nn(2),dndc(2),dmdc(2),f(2),mu,x(4),
     &                      d2adc2(2),d2cdc2(2),d2ndc2(2),d2mdc2(2),
     &                              d2bdc2(2),dfdmu(2),dfdc(2),s(3),
     &                              aa,bb,cc,dd,dphidx(3)
      real*8, intent(out) :: d2phidx2(3,3)
      integer, intent(in) :: iareaflag
c     real*8, parameter :: TOL1=0.996515d0       ! =44.9deg
c     real*8, parameter :: TOL2=1.003497d0       ! =45.1deg
c
      integer i,j
      real*8 daadx(3),dbbdx(3),dccdx(3),ddddx(3),dmudx(3),
     &       d2fdmu2(2),d2fdmudc(2),d2fdcdmu(2),d2fdc2(2)
      real*8 vcommon,vtmp1,vtmp2,vtmp3,vtmp4
      real*8 va,vc,dvadx(3),dsedx(3),dvcdx(3)
c
c                                           ---- calc.  dmudx(i) (i=1~3)
      daadx(1)=-a(2)-c(2)+2.0d0*b(2)
      dbbdx(1)=-2.0d0*(b(2)-a(2))
      dccdx(1)=-a(2)
c
      daadx(2)=a(1)+c(1)-2.0d0*b(1)
      dbbdx(2)=2.0d0*(b(1)-a(1))
      dccdx(2)=a(1)
c
      daadx(3)=x(2)*(dadc(1)+dcdc(1)-2.0d0*dbdc(1))
     &        -x(1)*(dadc(2)+dcdc(2)-2.0d0*dbdc(2))
      dbbdx(3)=2.0d0*x(2)*
     &        (dbdc(1)-dadc(1))-2.0d0*x(1)*(dbdc(2)-dadc(2))
      dccdx(3)=x(2)*dadc(1)-x(1)*dadc(2)
c
      do i=1,3
       dmudx(i)=0.5d0*daadx(i)*(bb+sqrt(dd))/(aa*aa)
     &         +0.5d0*(-dbbdx(i)-0.5d0/(sqrt(dd))*(2.0d0*bb*dbbdx(i)
     &           -4.0d0*daadx(i)*cc-4.0d0*aa*dccdx(i)))/aa
      enddo
c
c                                         ---- calc.  d2fdmu2(i) (i=1~2)
      do i=1,2
         d2fdmu2(i)=2.0d0*(a(i)+c(i)-2.0d0*b(i))
      enddo
c
c                                        ---- calc.  d2fdmudc(i) (i=1~2)
      do i=1,2
        d2fdmudc(i)=2.0d0*(dbdc(i)-dadc(i))+2.0d0*mu*(dadc(i)+dcdc(i)
     &                                             -2.0d0*dbdc(i))
      enddo
c
c                                        ---- calc.  d2fdcdmu(i) (i=1~2)
      do i=1,2
        d2fdcdmu(i)=2.0d0*(dbdc(i)-dadc(i))+2.0d0*mu*(dadc(i)+dcdc(i)
     &                                             -2.0d0*dbdc(i))
      enddo
c
c                                          ---- calc.  d2fdc2(i) (i=1~2)
      do i=1,2
        d2fdc2(i)=d2adc2(i)+2.0d0*mu*(d2bdc2(i)-d2adc2(i))
     &            +mu*mu*(d2adc2(i)+d2cdc2(i)-2.0d0*d2bdc2(i))
      enddo
c
c                                                 ---- for d2phidx2(k,l)
c
      vcommon=1.0d0/(f(1)*dfdmu(2)-f(2)*dfdmu(1))
      vtmp1=dfdc(1)*dfdmu(2)+f(1)*d2fdmudc(2)
     &                           -dfdc(2)*dfdmu(1)-f(2)*d2fdmudc(1)
      vtmp2=dfdmu(1)*dfdmu(2)+f(1)*d2fdmu2(2)
     &                           -dfdmu(2)*dfdmu(1)-f(2)*d2fdmu2(1)
      vtmp3=d2fdcdmu(2)*dfdmu(1)+dfdc(2)*d2fdmu2(1)
     &                     -d2fdcdmu(1)*dfdmu(2)-dfdc(1)*d2fdmu2(2)
      vtmp4=d2fdc2(2)*dfdmu(1)+dfdc(2)*d2fdmudc(1)
     &                      -d2fdc2(1)*dfdmu(2)-dfdc(1)*d2fdmudc(2)
c
      va=vcommon
      vc=dfdc(2)*dfdmu(1)-dfdc(1)*dfdmu(2)
c
      do i=1,2
        dvadx(i)=-vtmp2*vcommon*vcommon*dmudx(i)
      enddo
        dvadx(3)=-vtmp1*vcommon*vcommon
     &           -vtmp2*vcommon*vcommon*dmudx(3)
c
      do i=1,3
        dsedx(i)=dphidx(i)
      enddo
c
      do i=1,2
        dvcdx(i)=vtmp3*dmudx(i)
      enddo
        dvcdx(3)=vtmp4+vtmp3*dmudx(3)
c
c                                    ---- calc.  d2phidx2(i,j) (i,j=1~3)
      do j=1,2
        d2phidx2(1,j)=(-dfdmu(2)*vtmp2*vcommon*vcommon
     &                             +d2fdmu2(2)*vcommon)*dmudx(j)
      enddo
        d2phidx2(1,3)=-dfdmu(2)*vtmp1*vcommon*vcommon
     &                +d2fdmudc(2)*vcommon+(-dfdmu(2)*vtmp2
     &                  *vcommon*vcommon+d2fdmu2(2)*vcommon)*dmudx(3)
c
      do j=1,2
        d2phidx2(2,j)=(dfdmu(1)*vtmp2*vcommon*vcommon
     &                             -d2fdmu2(1)*vcommon)*dmudx(j)
      enddo
        d2phidx2(2,3)=dfdmu(1)*vtmp1*vcommon*vcommon
     &                -d2fdmudc(1)*vcommon+(dfdmu(1)*vtmp2
     &                  *vcommon*vcommon-d2fdmu2(1)*vcommon)*dmudx(3)
c
      do j=1,3
        d2phidx2(3,j)=dvadx(j)*se*vc+va*dsedx(j)*vc+va*se*dvcdx(j)
      enddo
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     calc. d2seds2(i,j) (i,j=1~3)
c
      subroutine jancae_vegter_calc_d2seds2 (d2seds2,d2phidx2,se,a,b,c,
     &      mu,x,vcos2t,vsin2t,iareaflag,dxds,dphidx,isflag,s,dseds,
     &      pryld,ndyld)
c-----------------------------------------------------------------------
      implicit none
c
      integer, intent(in) :: iareaflag,isflag,ndyld
      real*8, intent(in) :: d2phidx2(3,3),se,a(2),b(2),c(2),mu,x(4),
     &            vcos2t,vsin2t,dxds(3,3),dphidx(3),s(3),dseds(3),
     &            pryld(ndyld)
      real*8, intent(out) :: d2seds2(3,3)
      real*8, parameter :: TOL1 = 0.996515d0       ! =44.9deg
      real*8, parameter :: TOL2 = 1.003497d0       ! =45.1deg
      real*8, parameter :: TOL3a= 1.0d-6
      real*8, parameter :: TOL3b= 1.0d0-TOL3a
      real*8, parameter :: TOL4 = 1.0d-7
      real*8, parameter :: TOL4a= 1.0d0-TOL4       !thera<0.012812deg
      real*8, parameter :: TOL4b=-1.0d0+TOL4       !thera>89.98719deg
      real*8 vtan,d2xds2(3,3,3),vx1x2
      integer i,j,k,l,iflag
c
c                             ---- if condition to apply numerical diff.
c
      if(iareaflag.eq.3.or.iareaflag.eq.4) then
      vtan=x(2)/x(1)
      end if
c
      if(iareaflag.eq.4.and.vtan.ge.TOL1.and.vtan.le.TOL2) then
      iflag=1
c
      else if(iareaflag.eq.3.and.vtan.ge.TOL1.and.vtan.le.TOL2) then
      iflag=2
c
      else if(abs(mu).le.TOL3a.or.abs(mu).ge.TOL3b) then
      iflag=3
c
      else if(vcos2t.ge.TOL4a.or.vcos2t.le.TOL4b) then
      iflag=4
c
      else
      vx1x2=x(1)-x(2)
      iflag=0
      endif
c
c                                    ---- set d2xds2(k,i,j)  (k,i,j=1~3)
c
      if (iflag.eq.0) then
c
      d2xds2(1,1,1)=0.5d0*(1.0d0-vcos2t*vcos2t)/vx1x2
      d2xds2(1,1,2)=0.5d0*(vcos2t*vcos2t-1.0d0)/vx1x2
      d2xds2(1,1,3)=-vcos2t*vsin2t/vx1x2
c
      d2xds2(1,2,1)=0.5d0*(vcos2t*vcos2t-1.0d0)/vx1x2
      d2xds2(1,2,2)=0.5d0*(1-vcos2t*vcos2t)/vx1x2
      d2xds2(1,2,3)=vcos2t*vsin2t/vx1x2
c
      d2xds2(1,3,1)=-vcos2t*vsin2t/vx1x2
      d2xds2(1,3,2)=vcos2t*vsin2t/vx1x2
      d2xds2(1,3,3)=-2.0d0*(vsin2t*vsin2t-1.0d0)/vx1x2
c
      d2xds2(2,1,1)=0.5d0*(vcos2t*vcos2t-1.0d0)/vx1x2
      d2xds2(2,1,2)=0.5d0*(1.0d0-vcos2t*vcos2t)/vx1x2
      d2xds2(2,1,3)=vcos2t*vsin2t/vx1x2
c
      d2xds2(2,2,1)=0.5d0*(1.0d0-vcos2t*vcos2t)/vx1x2
      d2xds2(2,2,2)=0.5d0*(vcos2t*vcos2t-1.0d0)/vx1x2
      d2xds2(2,2,3)=-vcos2t*vsin2t/vx1x2
c
      d2xds2(2,3,1)=vcos2t*vsin2t/vx1x2
      d2xds2(2,3,2)=-vcos2t*vsin2t/vx1x2
      d2xds2(2,3,3)=2.0d0*(vsin2t*vsin2t-1.0d0)/vx1x2
c
      d2xds2(3,1,1)=-3.0d0*vcos2t*vsin2t*vsin2t/(vx1x2*vx1x2)
      d2xds2(3,1,2)=3.0d0*vcos2t*vsin2t*vsin2t/(vx1x2*vx1x2)
      d2xds2(3,1,3)=2.0d0*vsin2t*(2.0d0-3.0d0*vsin2t*vsin2t)/
     &                                           (vx1x2*vx1x2)
c
      d2xds2(3,2,1)=3.0d0*vcos2t*vsin2t*vsin2t/(vx1x2*vx1x2)
      d2xds2(3,2,2)=-3.0d0*vcos2t*vsin2t*vsin2t/(vx1x2*vx1x2)
      d2xds2(3,2,3)=-2.0d0*vsin2t*(2.0d0-3.0d0*vsin2t*vsin2t)/
     &                                           (vx1x2*vx1x2)
c
      d2xds2(3,3,1)=2.0d0*vsin2t*(3.0d0*vcos2t*vcos2t-1.0d0)/
     &                                           (vx1x2*vx1x2)
      d2xds2(3,3,2)=-2.0d0*vsin2t*(3.0d0*vcos2t*vcos2t-1.0d0)/
     &                                           (vx1x2*vx1x2)
      d2xds2(3,3,3)=4.0d0*vcos2t*(3.0d0*vsin2t*vsin2t-1.0d0)/
     &                                           (vx1x2*vx1x2)
      end if
c
c                                      ---- calc. d2seds2(i,j) (i,j=1~3)
c
      call jancae_clear2 ( d2seds2,3,3 )
c
      if (iflag.ne.0) then
      call jancae_vegter_d2seds2n(d2seds2,s,dseds,pryld,ndyld,se)
c
      else
      do i=1,3
        do j=1,3
          do k=1,3
            do l=1,3
              d2seds2(i,j)=d2seds2(i,j)
     &                    +d2phidx2(k,l)*dxds(l,j)*dxds(k,i)
            enddo
              d2seds2(i,j)=d2seds2(i,j)
     &                    +dphidx(k)*d2xds2(k,i,j)
          enddo
        enddo
      enddo
      endif
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     numerical differential for 2nd order differentials
c
      subroutine jancae_vegter_d2seds2n (d2seds2,s,dseds,
     &                                   pryld,ndyld,se)
c-----------------------------------------------------------------------
      implicit none
c
      integer, intent(in)  :: ndyld
      real*8, intent(out) :: d2seds2(3,3)
      real*8, intent(in)  :: dseds(3),pryld(ndyld),se,s(3)
      real*8, parameter :: delta=1.0d-3
      real*8 s0(3),sea,seb,a,b,seba,seaa,sebb,seab,se0,ss(3)
      integer j,k
c
      s0(:) = s(:)
      ss(:) = s(:)
        do j=1,3
          do k=1,3
           if ( j.eq.k ) then
             se0=se
             ss(j)=s0(j)-delta
             call jancae_vegter_yieldfunc(3,ss,sea,dseds,d2seds2,0,
     &                                    pryld,ndyld)
             ss(j)=s0(j)+delta
             call jancae_vegter_yieldfunc(3,ss,seb,dseds,d2seds2,0,
     &                                    pryld,ndyld)
             ss(j)=s0(j)
             a=(se0-sea)/delta
             b=(seb-se0)/delta
             d2seds2(j,k)=(b-a)/delta
           else
             ss(j)=s0(j)-delta
             ss(k)=s0(k)-delta
             call jancae_vegter_yieldfunc(3,ss,seaa,dseds,d2seds2,0,
     &                                    pryld,ndyld)
             ss(j)=s0(j)+delta
             ss(k)=s0(k)-delta
             call jancae_vegter_yieldfunc(3,ss,seba,dseds,d2seds2,0,
     &                                    pryld,ndyld)
             ss(j)=s0(j)-delta
             ss(k)=s0(k)+delta
             call jancae_vegter_yieldfunc(3,ss,seab,dseds,d2seds2,0,
     &                                    pryld,ndyld)
             ss(j)=s0(j)+delta
             ss(k)=s0(k)+delta
             call jancae_vegter_yieldfunc(3,ss,sebb,dseds,d2seds2,0,
     &                                    pryld,ndyld)
             ss(j)=s0(j)
             ss(k)=s0(k)
             a=(seba-seaa)/(2.0d0*delta)
             b=(sebb-seab)/(2.0d0*delta)
             d2seds2(j,k)=(b-a)/(2.0d0*delta)
           endif
          end do
         end do
c
      return
      end
c
c
c
c---------------------------------------------------------------(vegter)
c     calc. equivalent stress for d2seds2n
c
      subroutine jancae_vegter_yieldfunc (nttl,s,se,dseds,d2seds2,
     &                                    nreq,pryld,ndyld)
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: nttl,nreq,ndyld
      real*8, intent(in) :: pryld(ndyld),s(3)
      real*8, intent(out)  :: se,d2seds2(3,3),dseds(3)
c
      call jancae_vegter ( s,se,dseds,d2seds2,nreq,pryld,ndyld )
c
      return
      end
c
c
c
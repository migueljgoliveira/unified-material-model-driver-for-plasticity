c***********************************************************************
c
c     UMMDp : Uncoupled Rupture Criteria
c
c**********************************************************************
c
c      0 : No Rupture Criterion
c
c      1 : Equivalent Plastic Strain
c      2 : Cockroft and Latham
c      3 : Rice and Tracey
c      4 : Ayada
c      5 : Brozzo
c      6 : Forming Limit Diagram (only plane-stress)
c
c-----------------------------------------------------------------------
c     calculated rupture criteria
c
      subroutine jancae_rupture ( ntens,sdv,nsdv,uvar2,uvar1,nuvarm,
     &                            jrcd,jmac,jmatyp,matlayo,laccfla,
     &                            nt,ndrup,prrup )
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension UVAR1(NUVARM),jmac(*),jmatyp(*)
      dimension sdv(nsdv),uvar2(nuvarm),prrup(ndrup)
      real*8 lim,wlimnorm
c-----------------------------------------------------------------------
c
c      prrup(1) : criteria id
c      prrup(2) : flag to terminate analysis if limit is reached
c      prrup(3) : rupture limit
c
c 																		       ---- rupture criteria limit
      lim = prrup(3)
c                                           ---- select rupture criteria
      ntrup = nint(prrup(1))
      select case ( ntrup )
c
      case ( 0 )                                  ! No Rupture Criterion
        return
c
      case ( 1 )                             ! Equivalent Plastic Strain
        call jancae_rup_eqstrain ( sdv,nsdv,uvar2,uvar1,nuvarm,
     &                             nt,lim,wlimnorm )
c
      case ( 2 )                                   ! Cockroft and Latham
        call jancae_rup_cockroft ( sdv,nsdv,uvar2,uvar1,nuvarm,
     &                             jrcd,jmac,jmatyp,matlayo,laccfla,
     &                             nt,lim,wlimnorm )
c
      case ( 3 )                                       ! Rice and Tracey
        call jancae_rup_rice ( sdv,nsdv,uvar2,uvar1,nuvarm,
     &                         jrcd,jmac,jmatyp,matlayo,laccfla,
     &                         nt,lim,wlimnorm )
c
      case ( 4 )                                                 ! Ayada
        call jancae_rup_ayada ( sdv,nsdv,uvar2,uvar1,nuvarm,
     &                          jrcd,jmac,jmatyp,matlayo,laccfla,
     &                          nt,lim,wlimnorm )
c
      case ( 5 )                                                ! Brozzo
        call jancae_rup_brozzo ( sdv,nsdv,uvar2,uvar1,nuvarm,
     &                           jrcd,jmac,jmatyp,matlayo,laccfla,
     &                           nt,lim,wlimnorm )
c
      case ( 6 )                                 ! Forming Limit Diagram
        call jancae_rup_fld ( ntens,uvar2,uvar1,nuvarm,
     &                        jrcd,jmac,jmatyp,matlayo,laccfla,
     &                        nt,lim,wlimnorm )
c
      case default
        write (6,*) 'error in jancae_rupture'
        write (6,*) 'ntrup error :',ntrup
        call jancae_exit ( 9000 )
      end select
c
c                    ---- terminate analysis if rupture limit is reached
      end = nint(prrup(2))
      if ( end == 1 ) then
        if ( wlimnorm >= 1.0d0 ) then 
          write (6,*) 'analysis terminated by rupture criterion'
          write (6,*) 'stop in uvrm.'
          call jancae_exit( 10000 )
        end if
      end if
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     Equivalent Plastic Strain
c
      subroutine jancae_rup_eqstrain ( sdv,nsdv,uvar2,uvar1,nuvarm,
     &                                 nt,lim,wlimnorm )
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension sdv(nsdv),uvar2(nuvarm),uvar1(nuvarm)
      real*8 lim,peeq
c-----------------------------------------------------------------------
c
c     nuvarm : 2
c
c     uvar(2+nt+1) : equivalent plastic strain
c     uvar(2+nt+2) : rupture parameter normalised by critical value
c
c                                            ---- get uvar before update
      wlimnorm  = uvar1(2+nt+2)
c
c                                              ---- get sdv after update
      peeq = sdv(1)
c
c                                                       ---- update uvar
      uvar2(2+nt+1) = peeq
      uvar2(2+nt+2) = peeq / lim
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     Cockroft and Latham
c
      subroutine jancae_rup_cockroft ( sdv,nsdv,uvar2,uvar1,nuvarm,
     &                                 jrcd,jmac,jmatyp,matlayo,laccfla,
     &                                 nt,lim,wlimnorm )
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension UVAR1(NUVARM),jmac(*),jmatyp(*)
      dimension array(15),jarray(15)
      dimension sdv(nsdv),uvar2(nuvarm)
      character*3 flgray(15)
      real*8 lim,wlimnorm
      real*8 se1,peeq1,maxsp1,maxsp1se1,wlim1
      real*8 se2,peeq2,maxsp2,maxsp2se2,wlim2
c-----------------------------------------------------------------------
c
c     nuvarm : 4
c
c     uvar(1)      : equivalent stress
c     uvar(2+nt+1) : equivalent plastic strain
c     uvar(2+nt+2) : maximum principal stress
c     uvar(2+nt+3) : rupture parameter
c     uvar(2+nt+4) : rupture parameter normalised by critical value
c
c                                            ---- get uvar before update
      se1    = uvar1(1)
      peeq1  = uvar1(2+nt+1)
      maxsp1 = uvar1(2+nt+2)
      wlim1  = uvar1(2+nt+3)
c
      wlimnorm  = uvar1(2+nt+4)
c
c                                     ---- get sdv and uvar after update
      se2 = uvar2(1)
      peeq2 = sdv(1)
c
c                                 ---- get principal stress after update
      call getvrm ('SP',array,jarray,flgray,jrcd,jmac,jmatyp,
     &                  matlayo,laccfla )
      if ( jrcd /= 0 ) then
        write (6,*) 'request error in uvarm for sp'
        write (6,*) 'stop in uvrm.'
        call jancae_exit ( 9000 )
      end if
      maxsp2 = array(3)
c
c                                                 ---- rupture criterion
      maxsp1se1 = 0.0d0
      maxsp2se2 = 0.0d0
      if ( se1 > 0.0d0 ) maxsp1se1 = maxsp1 / se1
      if ( se2 > 0.0d0 ) maxsp2se2 = maxsp2 / se2
c
      wlim2 = wlim1 + (maxsp2se2+maxsp1se1)*(peeq2-peeq1)/2.0d0
c
c                                                       ---- update uvar
      uvar2(2+nt+1) = peeq2
      uvar2(2+nt+2) = maxsp2
      uvar2(2+nt+3) = wlim2
      uvar2(2+nt+4) = wlim2/lim
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     Rice and Tracey
c
      subroutine jancae_rup_rice ( sdv,nsdv,uvar2,uvar1,nuvarm,
     &                             jrcd,jmac,jmatyp,matlayo,laccfla,
     &                             nt,lim,wlimnorm )
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension UVAR1(NUVARM),jmac(*),jmatyp(*)
      dimension array(15),jarray(15)
      dimension sdv(nsdv),uvar2(nuvarm)
      character*3 flgray(15)
      real*8 lim,wlimnorm
      real*8 se1,peeq1,shyd1,shyd1se1,wlim1
      real*8 se2,peeq2,shyd2,shyd2se2,wlim2
c-----------------------------------------------------------------------
c
c     nuvarm : 4
c
c     uvar(1)      : equivalent stress
c     uvar(2+nt+1) : equivalent plastic strain
c     uvar(2+nt+2) : hydrostatic stress
c     uvar(2+nt+3) : rupture parameter
c     uvar(2+nt+4) : rupture parameter normalised by critical value
c
c                                            ---- get uvar before update
      se1   = uvar1(1)
      peeq1 = uvar1(2+nt+1)
      shyd1 = uvar1(2+nt+2)
      wlim1 = uvar1(2+nt+3)
c
      wlimnorm  = uvar1(2+nt+4)
c
c                                     ---- get sdv and uvar after update
      se2 = uvar2(1)
      peeq2 = sdv(1)
c
c                               ---- get hydrostatic stress after update
      call getvrm ( 'SINV',array,jarray,flgray,jrcd,jmac,jmatyp,
     &                     matlayo,laccfla )
      if ( jrcd /= 0 ) then
        write (6,*) 'request error in uvarm for sinv'
        write (6,*) 'stop in uvrm.'
        call jancae_exit ( 9000 )
      end if
      shyd2 = -array(3)
c
c                                                 ---- rupture criterion
      shyd1se1 = 0.0d0
      shyd2se2 = 0.0d0
      if ( se1 > 0.0d0 ) shyd1se1 = exp(1.5d0*shyd1/se1)
      if ( se2 > 0.0d0 ) shyd2se2 = exp(1.5d0*shyd2/se2)
c
      wlim2 = wlim1 + (shyd1se1+shyd2se2)*(peeq2-peeq1)/2.0d0
c
c                                                       ---- update uvar
      uvar2(2+nt+1) = peeq2
      uvar2(2+nt+2) = shyd2
      uvar2(2+nt+3) = wlim2
      uvar2(2+nt+4) = wlim2/lim
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     Ayada
c
      subroutine jancae_rup_ayada ( sdv,nsdv,uvar2,uvar1,nuvarm,
     &                              jrcd,jmac,jmatyp,matlayo,laccfla,
     &                              nt,lim,wlimnorm )
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension UVAR1(NUVARM),jmac(*),jmatyp(*)
      dimension array(15),jarray(15)
      dimension sdv(nsdv),uvar2(nuvarm)
      character*3 flgray(15)
      real*8 lim,wlimnorm
      real*8 se1,peeq1,shyd1,shyd1se1,wlim1
      real*8 se2,peeq2,shyd2,shyd2se2,wlim2
c-----------------------------------------------------------------------
c
c     nuvarm : 4
c
c     uvar(1)      : equivalent stress
c     uvar(2+nt+1) : equivalent plastic strain
c     uvar(2+nt+2) : hydrostatic stress
c     uvar(2+nt+3) : rupture parameter
c     uvar(2+nt+4) : rupture parameter normalised by critical value
c
c                                            ---- get uvar before update
      se1   = uvar1(1)
      peeq1 = uvar1(2+nt+1)
      shyd1 = uvar1(2+nt+2)
      wlim1 = uvar1(2+nt+3)
c
      wlimnorm  = uvar1(2+nt+4)
c
c                                     ---- get sdv and uvar after update
      se2 = uvar2(1)
      peeq2 = sdv(1)
c
c                               ---- get hydrostatic stress after update
      call getvrm ( 'SINV',array,jarray,flgray,jrcd,jmac,jmatyp,
     &                     matlayo,laccfla )
      if ( jrcd /= 0 ) then
        write (6,*) 'request error in uvarm for sinv'
        write (6,*) 'stop in uvrm.'
        call jancae_exit ( 9000 )
      end if
      shyd2 = -array(3)
c
c                                                 ---- rupture criterion
      shyd1se1 = 0.0d0
      shyd2se2 = 0.0d0
      if ( se1 > 0.0d0 ) shyd1se1 = shyd1 / se1
      if ( se2 > 0.0d0 ) shyd2se2 = shyd2 / se2
c
      wlim2 = wlim1 + (shyd1se1+shyd2se2)*(peeq2-peeq1)/2.0d0
c
c                                                       ---- update uvar
      uvar2(2+nt+1) = peeq2
      uvar2(2+nt+2) = shyd2
      uvar2(2+nt+3) = wlim2
      uvar2(2+nt+4) = wlim2/lim
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     Brozzo
c
      subroutine jancae_rup_brozzo ( sdv,nsdv,uvar2,uvar1,nuvarm,
     &                               jrcd,jmac,jmatyp,matlayo,laccfla,
     &                               nt,lim,wlimnorm )
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension UVAR1(NUVARM),jmac(*),jmatyp(*)
      dimension array(15),jarray(15)
      dimension sdv(nsdv),uvar2(nuvarm)
      character*3 flgray(15)
      real*8 lim,wlimnorm
      real*8 se1,peeq1,shyd1,maxsp1,maxsp1shyd1,wlim1
      real*8 se2,peeq2,shyd2,maxsp2,maxsp2shyd2,wlim2
c-----------------------------------------------------------------------
c
c     nuvarm : 5
c
c     uvar(2+nt+1) : equivalent plastic strain
c     uvar(2+nt+2) : maximum principal stress
c     uvar(2+nt+3) : hydrostatic stress
c     uvar(2+nt+4) : rupture parameter
c     uvar(2+nt+5) : rupture parameter normalised by critical value
c
c                                            ---- get uvar before update
      peeq1  = uvar1(2+nt+1)
      maxsp1 = uvar1(2+nt+2)
      shyd1  = uvar1(2+nt+3)
      wlim1  = uvar1(2+nt+4)
c
      wlimnorm  = uvar1(2+nt+5)
c
c                                              ---- get sdv after update
      peeq2 = sdv(1)
c
c                                 ---- get principal stress after update
      call getvrm ('SP',array,jarray,flgray,jrcd,jmac,jmatyp,
     &                  matlayo,laccfla )
      if ( jrcd /= 0 ) then
        write (6,*) 'request error in uvarm for sp'
        write (6,*) 'stop in uvrm.'
        call jancae_exit ( 9000 )
      end if
      maxsp2 = array(3)
c
c                               ---- get hydrostatic stress after update
      call getvrm ('SINV',array,jarray,flgray,jrcd,jmac,jmatyp,
     &                    matlayo,laccfla )
      if ( jrcd /= 0 ) then
        write (6,*) 'request error in uvarm for sinv'
        write (6,*) 'stop in uvrm.'
        call jancae_exit ( 9000 )
      end if
      shyd2 = -array(3)
c
c                                                 ---- rupture criterion
      maxsp1shyd1 = 0.0d0
      maxsp2shyd2 = 0.0d0
      if ( shyd1 > 0.0d0 ) then
        maxsp1shyd1 = (2.0d0/3.0d0) * maxsp1 / (maxsp1-shyd1)
      end if
      if ( shyd2 > 0.0d0 ) then 
        maxsp2shyd2 = (2.0d0/3.0d0) * maxsp2 / (maxsp2-shyd2)
      end if
c
      wlim2 = wlim1 + (maxsp2shyd2+maxsp1shyd1)*(peeq2-peeq1)/2.0d0
c
c                                                       ---- update uvar
      uvar2(2+nt+1) = peeq2
      uvar2(2+nt+2) = maxsp2
      uvar2(2+nt+3) = shyd2
      uvar2(2+nt+4) = wlim2
      uvar2(2+nt+5) = wlim2/lim
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     Forming Limit Diagram (FLD)
c       . only plane-stress formulation
c
      subroutine jancae_rup_fld ( ntens,uvar2,uvar1,nuvarm,
     &                            jrcd,jmac,jmatyp,matlayo,laccfla,
     &                            nt,lim,wlimnorm )
c
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      parameter (mxflc=50)
c
      dimension jmac(*),jmatyp(*)
      dimension array(15),jarray(15)
      dimension uvar2(nuvarm),uvar1(nuvarm)
      character*3 flgray(15)
      real*8 lim,wlimnorm
      real*8 e1,e2,e1fld,wlim
      real*8 le(3,3),es(3),ev(3,3)
      dimension dum1(0),dum2(mxflc,0)
      dimension fld1(mxflc),fld2(mxflc)
c-----------------------------------------------------------------------
c
c     nuvarm : 5
c
c     uvar(2+nt+1) : maximum principal strain
c     uvar(2+nt+2) : minimum principal strain
c     uvar(2+nt+3) : projection of major principal strain on FLD
c     uvar(2+nt+4) : rupture parameter
c     uvar(2+nt+5) : rupture parameter normalised by critical value
c
c                                            ---- get uvar before update
      wlimnorm  = uvar1(2+nt+6)
c      
c                                 ---- get principal strain after update
      if ( ntens == 3 ) then
        call getvrm ( 'LEP',array,jarray,flgray,jrcd,jmac,jmatyp,
     &                      matlayo,laccfla )
        if ( jrcd /= 0 ) then
          write (6,*) 'request error in uvarm for lep'
          write (6,*) 'stop in uvrm.'
          call jancae_exit ( 9000 )
        end if
        e2 = array(1)
        e1 = array(2)
c                               ---- get logarithmic strain after update
      else
        write (6,*) 'request error in uvarm for fld, only plane-stress'
        write (6,*) 'stop in uvrm.'
        call jancae_exit ( 9000 )
!         call getvrm ( 'LE',array,jarray,flgray,jrcd,jmac,jmatyp,
!      &                     matlayo,laccfla )
!         if ( jrcd /= 0 ) then
!           write (6,*) 'request error in uvarm for le'
!           write (6,*) 'stop in uvrm.'
!           call jancae_exit ( 9000 )
!         end if
! c                                            ---- assemble strain tensor
!         call jancae_clear2 ( le,3,3 )
!         le(1,1) = array(1)
!         le(2,2) = array(2)
!         le(3,3) = array(3)
!         le(1,2) = array(4)/2
!         le(1,3) = array(5)/2
!         le(2,3) = array(6)/2
!         le(2,1) = le(1,2)
!         le(3,1) = le(1,3)
!         le(3,2) = le(2,3)
! c                           ---- strain tensor eigen- values and vectors
!         call jancae_clear1 ( es,3 )
!         call jancae_clear2 ( ev,3,3 )
!         call jancae_eigen_sym3 ( es,ev,le )
!         e2 = es(2)
!         e1 = es(1)
      end if
c
c                                     ---- activate fld table collection
      call settablecollection ( 'FLD',jerror )

      if ( jerror /= 0 ) then
        write (6,*) 'request error in uvarm for table collection fld'
        write (6,*) 'stop in uvrm.'
        call jancae_exit ( 9000 )
      end if
c	                     		                        ---- get fld E1 values
      call getpropertytable ( 'FLD1',dum1,dum1,dum1,nfld,fld1,dum2,0,
     &                               jerror )
      if ( jerror /= 0 ) then
        write (6,*) 'request error in uvarm for property table fld1'
        write (6,*) 'stop in uvrm.'
        call jancae_exit ( 9000 )
      end if
c                                                 ---- get fld E2 values
      call getpropertytable ( 'FLD2',dum1,dum1,dum1,nfld,fld2,dum2,0,
     &                               jerror )
      if ( jerror /= 0 ) then
        write (6,*) 'request error in uvarm for property table fld2'
        write (6,*) 'stop in uvrm.'
        call jancae_exit ( 9000 )
      end if
c
c                         ---- linear extra/inter -polation of E1 on FLD
      n = nfld
c                              --- linear extrapolation on the left side
      if ( e2 < fld2(1) ) then
        e1fld = fld1(2) + ( (e2-fld2(2)) / (fld2(1)-fld2(2)) )
     &          * ( fld1(1) - fld1(2) )
c
c                             --- linear extrapolation on the right side
      else if ( e2 > fld2(n) ) then
        e1fld = fld1(n-1) + ( (e2-fld2(n-1)) / (fld2(n)-fld2(n-1)) ) 
     &          * ( fld1(n) - fld1(n-1) )
c
c                                  --- linear interpolation inside range
      else
        k = 0
        do i = 1,n-1
          if ( ( e2 >= fld2(i) ) .and. ( e2 <= fld2(i+1) ) ) then
            k = i
          end if
        end do
        e1fld = fld1(k) + ( fld1(k+1) - fld1(k) )
     &          * ( (e2-fld2(k)) / (fld2(k+1)-fld2(k)) )
      end if
c
c                                                 ---- rupture criterion
      wlim = e1/e1fld
c
c                                                       ---- update uvar
      uvar2(2+nt+1) = e1
      uvar2(2+nt+2) = e2
      uvar2(2+nt+3) = e1fld
      uvar2(2+nt+4) = wlim
      uvar2(2+nt+5) = wlim/lim
c
      return
      end
c
c
c
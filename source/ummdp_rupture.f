c***********************************************************************
c     JANCAE/UMMDp : Rupture Criteria
c**********************************************************************
c
c      0 : No Rupture Criterion
c
c      1 : Equivalent Plastic Strain
c      2 : Cockroft and Latham (CL)
c      3 : Rice and Tracey (RT)
c      4 : Ayada
c      5 : Brozzo
c
c-----------------------------------------------------------------------
c     calc. rupture criterion
c
      subroutine jancae_rupture ( sdv,nsdv,uvar2,uvar1,nuvarm,jrcd,jmac,
     1                            jmatyp,matlayo,laccfla,nt,
     2                            ndrup,prrup )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension UVAR1(NUVARM),JMAC(*),JMATYP(*)
      dimension sdv(nsdv),uvar2(nuvarm),prrup(ndrup)
c
      ntrup = prrup(1)
c
      select case ( ntrup )
      case ( 0 )                                  ! No Rupture Criterion
        return
c
      case ( 1 )                             ! Equivalent Plastic Strain
        call jancae_rup_eqstrain ( sdv,nsdv,uvar2,nuvarm,nt,
     1                             ndrup,prrup )
c
      case ( 2 )                                   ! Cockroft and Latham
        call jancae_rup_cockroft ( sdv,nsdv,uvar2,uvar1,nuvarm,jrcd,
     1                             jmac,jmatyp,matlayo,laccfla,nt,
     2                             ndrup,prrup )
c
      case ( 3 )                                       ! Rice and Tracey
        call jancae_rup_rice ( sdv,nsdv,uvar2,uvar1,nuvarm,jrcd,
     1                         jmac,jmatyp,matlayo,laccfla,nt,
     2                         ndrup,prrup )
c
      case ( 4 )                                                 ! Ayada
        call jancae_rup_ayada ( sdv,nsdv,uvar2,uvar1,nuvarm,jrcd,
     1                          jmac,jmatyp,matlayo,laccfla,nt,
     2                          ndrup,prrup )
     c
      case ( 5 )                                                ! Brozzo
        call jancae_rup_brozzo ( sdv,nsdv,uvar2,uvar1,nuvarm,jrcd,
     1                           jmac,jmatyp,matlayo,laccfla,nt,
     2                           ndrup,prrup )
c
      case default
        write (6,*) 'error in jancae_rupture'
        write (6,*) 'ntrup error :',ntrup
        call jancae_exit (9000)
      end select
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     print parameters for uncoupled rupture criterion
c
      subroutine jancae_rupture_print ( prrup,ndrup )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension prrup(ndrup)
c
      ntrup = nint(prrup(1))
      write (6,*)
      write (6,*) '*** Uncoupled Rupture Criterion',ntrup
      select case ( ntrup )
      case ( 0 )
        write (6,*) 'No Uncoupled Rupture Criterion'
      case ( 1 )
        write (6,*) 'Eq. Plastic Strain'
        write (6,*) 'W=int[dp]'
      case ( 2 )
        write (6,*) 'Cockroft and Latham'
        write (6,*) 'W=int[(sp1/se)*dp]'
      case ( 3 )
        write (6,*) 'Rice and Tracey'
        write (6,*) 'W=int[exp(1.5*sh/se)*dp]'
      case ( 4 )
        write (6,*) 'Ayada'
        write (6,*) 'W=int[(sh/se)*dp]'
      case ( 5 )
        write (6,*) 'Brozzo'
        write (6,*) 'W=int[(2/3)*(sp1/(sp1-se))*dp]'
      end select
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     Equivalent Plastic Strain
c
      subroutine jancae_rup_eqstrain ( sdv,nsdv,uvar2,nuvarm,nt,
     1                                 ndrup,prrup)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension sdv(nsdv),uvar2(nuvarm),prrup(ndrup)
      real*8 peeq
c
c     uvar(2+nt+1) : eq. plastic strain
c     prrup(1)     : id
c
c                                              ---- get sdv after update
      peeq = sdv(1)
c
c                                                       ---- update uvar
      uvar2(2+nt+1) = peeq
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     Cockroft and Latham (CL)
c
      subroutine jancae_rup_cockroft ( sdv,nsdv,uvar2,uvar1,nuvarm,jrcd,
     1                                 jmac,jmatyp,matlayo,laccfla,nt,
     2                                 ndrup,prrup )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension UVAR1(NUVARM),JMAC(*),JMATYP(*)
      dimension ARRAY(15),JARRAY(15),prrup(ndrup)
      dimension sdv(nsdv),uvar2(nuvarm)
      character*3 FLGRAY(15)
      real*8 se1,peeq1,maxsp1,maxsp1se1,wlim1
      real*8 se2,peeq2,maxsp2,maxsp2se2,wlim2
c
c     uvar(2+nt+1) : equivalent plastic strain
c     uvar(2+nt+2) : maximum principal stress
c     uvar(2+nt+3) : rupture criterion limit
c
c     prrup(1) : id
c
c                                            ---- get uvar before update
      se1    = uvar1(1)
      peeq1  = uvar1(2+nt+1)
      maxsp1 = uvar1(2+nt+2)
      wlim1  = uvar1(2+nt+3)
c
c                                     ---- get sdv and uvar after update
      se2 = uvar2(1)
      peeq2 = sdv(1)
c
c                                 ---- get principal stress after update
      CALL GETVRM('SP',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1            MATLAYO,LACCFLA)
      if ( JRCD.ne.0 ) then
        write (6,*) 'request error in uvarm for sp'
        write (6,*) 'stop in uvrm.'
        call jancae_exit ( 9000 )
      endif
      maxsp2 = array(3)
c
c                                                 ---- rupture criterion
      maxsp1se1 = 0.0d0
      maxsp2se2 = 0.0d0
      if ( se1 .gt. 0.0d0 ) then
        maxsp1se1 = maxsp1/se1
      endif
      if ( se2 .gt. 0.0d0 ) then
        maxsp2se2 = maxsp2/se2
      endif
c
      wlim2 = wlim1 + (maxsp2se2+maxsp1se1)*(peeq2-peeq1)/2.0d0
c                                                       ---- update uvar
      uvar2(2+nt+1) = peeq2
      uvar2(2+nt+2) = maxsp2
      uvar2(2+nt+3) = wlim2
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     Rice and Tracey (RT)
c
      subroutine jancae_rup_rice ( sdv,nsdv,uvar2,uvar1,nuvarm,jrcd,
     1                             jmac,jmatyp,matlayo,laccfla,nt,
     2                             ndrup,prrup )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension UVAR1(NUVARM),JMAC(*),JMATYP(*)
      dimension ARRAY(15),JARRAY(15),prrup(ndrup)
      dimension sdv(nsdv),uvar2(nuvarm)
      character*3 FLGRAY(15)
      real*8 se1,peeq1,shyd1,shyd1se1,wlim1
      real*8 se2,peeq2,shyd2,shyd2se2,wlim2
c
c     uvar(2+nt+1) : eq. plastic strain
c     uvar(2+nt+2) : hydrostatic stress
c     uvar(2+nt+3) : rupture criterion limit
c
c     prrup(1) : id
c
c                                            ---- get uvar before update
      se1   = uvar1(1)
      peeq1 = uvar1(2+nt+1)
      shyd1 = uvar1(2+nt+2)
      wlim1 = uvar1(2+nt+3)
c
c                                     ---- get sdv and uvar after update
      se2 = uvar2(1)
      peeq2 = sdv(1)
c
c                               ---- get hydrostatic stress after update
      CALL GETVRM('SINV',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1            MATLAYO,LACCFLA)
      if ( JRCD .ne. 0 ) then
        write (6,*) 'request error in uvarm for sinv'
        write (6,*) 'stop in uvrm.'
        call jancae_exit ( 9000 )
      endif
      shyd2 = -array(3)
c
c                                                 ---- rupture criterion
      shyd1se1 = 0.0d0
      shyd2se2 = 0.0d0
      if ( se1 .gt. 0.0d0 ) then
        shyd1se1 = exp(1.5d0*shyd1/se1)
      endif
      if ( se2 .gt. 0.0d0 ) then
        shyd2se2 = exp(1.5d0*shyd2/se2)
      endif
c
      wlim2 = wlim1 + (shyd1se1+shyd2se2)*(peeq2-peeq1)/2.0d0
c
c                                                       ---- update uvar
      uvar2(2+nt+1) = peeq2
      uvar2(2+nt+2) = shyd2
      uvar2(2+nt+3) = wlim2
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     Ayada
c
      subroutine jancae_rup_ayada ( sdv,nsdv,uvar2,uvar1,nuvarm,jrcd,
     1                              jmac,jmatyp,matlayo,laccfla,nt,
     2                              ndrup,prrup )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension UVAR1(NUVARM),JMAC(*),JMATYP(*)
      dimension ARRAY(15),JARRAY(15),prrup(ndrup)
      dimension sdv(nsdv),uvar2(nuvarm)
      character*3 FLGRAY(15)
      real*8 se1,peeq1,shyd1,shyd1se1,wlim1
      real*8 se2,peeq2,shyd2,shyd2se2,wlim2
c
c     uvar(2+nt+1) : eq. plastic strain
c     uvar(2+nt+2) : hydrostatic stress
c     uvar(2+nt+3) : rupture criterion limit
c
c     prrup(1) : id
c
c                                            ---- get uvar before update
      se1   = uvar1(1)
      peeq1 = uvar1(2+nt+1)
      shyd1 = uvar1(2+nt+2)
      wlim1 = uvar1(2+nt+3)
c
c                                     ---- get sdv and uvar after update
      se2 = uvar2(1)
      peeq2 = sdv(1)
c
c                               ---- get hydrostatic stress after update
      CALL GETVRM('SINV',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1            MATLAYO,LACCFLA)
      if ( JRCD.ne.0 ) then
        write (6,*) 'request error in uvarm for sinv'
        write (6,*) 'stop in uvrm.'
        call jancae_exit ( 9000 )
      endif
      shyd2 = -array(3)
c
c                                                 ---- rupture criterion
      shyd1se1 = 0.0d0
      shyd2se2 = 0.0d0
      if ( se1 .gt. 0.0d0 ) then
        shyd1se1 = hyd1/se1
      endif
      if ( se2 .gt. 0.0d0 ) then
        shyd2se2 = shyd2/se2
      endif
      wlim2 = wlim1 + (shyd1se1+shyd2se2)*(peeq2-peeq1)/2.d0
c
c                                                       ---- update uvar
      uvar2(2+nt+1) = peeq2
      uvar2(2+nt+2) = shyd2
      uvar2(2+nt+3) = wlim2
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     Brozzo
c
      subroutine jancae_rup_brozzo ( sdv,nsdv,uvar2,uvar1,nuvarm,jrcd,
     1                               jmac,jmatyp,matlayo,laccfla,nt,
     2                               ndrup,prrup )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension UVAR1(NUVARM),JMAC(*),JMATYP(*)
      dimension ARRAY(15),JARRAY(15),prrup(ndrup)
      dimension sdv(nsdv),uvar2(nuvarm)
      character*3 FLGRAY(15)
      real*8 se1,peeq1,maxsp1,maxsp1se1,wlim1
      real*8 se2,peeq2,maxsp2,maxsp2se2,wlim2
c
c     uvar(2+nt+1) : eq. plastic strain
c     uvar(2+nt+2) : max. principal stress
c     uvar(2+nt+3) : rupture criterion limit
c
c     prrup(1) : id
c
c                                            ---- get uvar before update
      se1    = uvar1(1)
      peeq1  = uvar1(2+nt+1)
      maxsp1 = uvar1(2+nt+2)
      wlim1  = uvar1(2+nt+3)
c
c                                     ---- get sdv and uvar after update
      se2 = uvar2(1)
      peeq2 = sdv(1)
c
c                                 ---- get principal stress after update
      CALL GETVRM('SP',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1            MATLAYO,LACCFLA)
      if ( JRCD.ne.0 ) then
        write (6,*) 'request error in uvarm for sp'
        write (6,*) 'stop in uvrm.'
        call jancae_exit ( 9000 )
      endif
      maxsp2 = array(3)
c
c                                                 ---- rupture criterion
      maxsp1se1 = 0.0d0
      maxsp2se2 = 0.0d0
      if ( se1 .gt. 0.0d0 ) then
        maxsp1se1 = (2.0d0/3.0d0) * maxsp1 / (maxsp1-se1)
      endif
      if ( se2 .gt. 0.0d0 ) then
        maxsp2se2 = (2.0d0/3.0d0) * maxsp2 / (maxsp2-se2)
      endif
c
      wlim2 = wlim1 + (maxsp2se2+maxsp1se1)*(peeq2-peeq1)/2.0d0
c
c                                                       ---- update uvar
      uvar2(2+nt+1) = peeq2
      uvar2(2+nt+2) = maxsp2
      uvar2(2+nt+3) = wlim2
c
      return
      end
c
c
c
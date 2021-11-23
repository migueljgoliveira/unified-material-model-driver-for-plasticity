c***********************************************************************
c
c     UMMDp : Yield Criteria
c
c***********************************************************************
c
c      0 : von Mises (1913)
c
c     3D
c       1 : Hill 1948
c       2 : Yld2004-18p
c       3 : CPB 2006
c       4 : Karafillis-Boyce 1993
c       5 : Hu 2005
c       6 : Yoshida 2011
c
c     2D
c      -1 : Gotoh Biquadratic
c      -2 : Yld2000-2d
c      -3 : Vegter
c      -4 : BBC 2005
c      -5 : Yld89
c      -6 : BBC 2008
c      -7 : Hill 1990
c
************************************************************************
c     CALCULATE YIELD FUNCTION AND DERIVATIVES
c
      subroutine jancae_yfunc ( se,cdseds,cd2seds2,nreq,
     &                          cs,nttl,nnrm,nshr,
     &                          pryld,ndyld )
c-----------------------------------------------------------------------
      implicit none
c
      integer nreq,nttl,nnrm,nshr,ndyld
      real*8 se
			real*8 cdseds(nttl),cs(nttl),pryld(ndyld)
			real*8 cd2seds2(nttl,nttl)
c       
      integer i,j
      integer ntyld
      integer indx(6)
			real*8 ss
	    real*8 s(6),dseds(6)
			real*8 d2seds2(6,6)
c-----------------------------------------------------------------------
c
      ntyld = nint(pryld(1))
c
      if ( ntyld < 0 ) then
        if ( ( nnrm /= 2 ) .or. ( nshr /= 1 ) ) then
          write (6,*) 'error in jancae_yfunc'
          write (6,*) 'ntyld<0 for plane stress'
          write (6,*) 'nnrm,nshr,ntyld:',nnrm,nshr,ntyld
          call jancae_exit (9000)
        end if
        goto 100
      end if
c
      ss = 0.0
      do i = 1,nttl
        ss = ss + cs(i)**2
      end do
      if ( ( ss <= 0.0 ) .and. ( nreq == 0 ) ) then
        se = 0.0
        return
      end if
c
c                                                ---- 3D yield functions
c
c                                        ---- set index to s(i) to cs(i)
      do i = 1,6
        indx(i) = 0
      end do
      if ( nnrm == 3 ) then
        do i = 1,nttl
          indx(i) = i
        end do
      else if ( nnrm == 2 ) then
        indx(1) = 1
        indx(2) = 2
        indx(3) = 0
        do i = 1,nshr
          indx(3+i) = 2 + i
        end do
      end if
c                                                          ---- set s(i)
      call jancae_clear1 ( s,6 )
      do i = 1,6
        if ( indx(i) /= 0 ) then
          s(i) = cs(indx(i))
        end if
      end do
c
      select case ( ntyld )
      case ( 0 )                                             ! von Mises
        call jancae_mises ( s,se,dseds,d2seds2,nreq )
c
      case ( 1 )                                             ! Hill 1948
        call jancae_hill1948 ( s,se,dseds,d2seds2,nreq,
     &                         pryld,ndyld )
c
      case ( 2 )                                           ! Yld2004-18p
        call jancae_yld2004_18p ( s,se,dseds,d2seds2,nreq,
     &                            pryld,ndyld )
c
      case ( 3 )                                              ! CPB 2006
        call jancae_cpb2006 ( s,se,dseds,d2seds2,nreq,
     &                        pryld,ndyld )
c
      case ( 4 )                                 ! Karafillis-Boyce 1993
        call jancae_KarafillisBoyce ( s,se,dseds,d2seds2,nreq,
     &                                pryld,ndyld )
c
      case ( 5 )                                               ! Hu 2005
        call jancae_hu2005 ( s,se,dseds,d2seds2,nreq,
     &                       pryld,ndyld )
c
      case ( 6 )                                          ! Yoshida 2011
        call jancae_yoshida2011 ( s,se,dseds,d2seds2,nreq,
     &                            pryld,ndyld )
c
      case default
        write (6,*) 'error in jancae_yfunc'
        write (6,*) 'ntyld error :',ntyld
        call jancae_exit (9000)
      end select
c
c                                                        ---- set dse/ds
      if ( nreq >= 1 ) then
        do i = 1,6
          if ( indx(i) /= 0 ) cdseds(indx(i)) = dseds(i)
        end do
      end if
c                                                      ---- set d2se/ds2
      if ( nreq >= 2 ) then
        do i = 1,6
          if ( indx(i) /= 0 ) then
            do j = 1,6
              if ( indx(j) /= 0 ) then
                cd2seds2(indx(i),indx(j)) = d2seds2(i,j)
              end if
            end do
          end if
        end do
      end if
c
      return
c
c
  100 continue
c                                       ---- plane stress yield criteria
c
      select case ( ntyld )
      case ( -1 )                                    ! Gotoh Biquadratic
        call jancae_gotoh ( cs,se,cdseds,cd2seds2,nreq,
     &                      pryld,ndyld )
c
      case ( -2 )                                           ! Yld2000-2d
        call jancae_yld2000 ( cs,se,cdseds,cd2seds2,nreq,
     &                        pryld,ndyld )
c
      case ( -3 )                                               ! Vegter
        call jancae_vegter ( cs,se,cdseds,cd2seds2,nreq,
     &                       pryld,ndyld )
c
      case ( -4 )                                             ! BBC 2005
        call jancae_bbc2005 ( cs,se,cdseds,cd2seds2,nreq,
     &                        pryld,ndyld )
c
      case ( -5 )                                                ! Yld89
        call jancae_yld89 ( cs,se,cdseds,cd2seds2,nreq,
     &                      pryld,ndyld )
c
      case ( -6 )                                             ! BBC 2008
        call jancae_bbc2008 ( cs,se,cdseds,cd2seds2,nreq,
     &                        pryld,ndyld )
c
      case ( -7 )                                            ! Hill 1990
        call jancae_hill1990 ( cs,se,cdseds,cd2seds2,nreq,
     &                         pryld,ndyld )
c
      case default
        write (6,*) 'error in jancae_yfunc'
        write (6,*) 'ntyld error :',ntyld
        call jancae_exit (9000)
      end select
c
      return
      end subroutine jancae_yfunc
c
c
c
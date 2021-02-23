c***********************************************************************
c     JANCAE/UMMDp : Yield Functions 
c***********************************************************************
c-----------------------------------------------------------------------
c     yield function and its dfferentials
c
      subroutine jancae_yfunc ( se,cdseds,cd2seds2,nreq,
     &                          cs,nttl,nnrm,nshr,
     &                          pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension cs(nttl),cdseds(nttl),cd2seds2(nttl,nttl),
     &          pryld(ndyld)
      dimension s(6),dseds(6),d2seds2(6,6),indx(6)
c
      ntyld = nint(pryld(1))
c
c     0 : von Mises isotropic (1913)
c
c     1 : Hill quadratic (1948)
c     2 : Barlat yld2004 (2005)
c     3 : Cazacu (2006)
c     4 : Karafillis-Boyce (1993)
c     5 : Hu (2005)
c     6 : Yohsida (2011)
c
c    -1 : Gotoh biquadratic (1978)
c    -2 : Barlat YLD2000-2d (2000)
c    -3 : Vegter
c    -4 : Banabic BBC2005
c    -5 : Barlat YLD89
c    -6 : Banabic BBC2008
c    -7 : Hill 1990
c
      if ( ntyld .lt. 0 ) then
        if ( ( nnrm .ne. 2 ) .or. ( nshr .ne. 1 ) ) then
          write (6,*) 'error in jancae_yfunc'
          write (6,*) 'ntyld<0 for plane stress'
          write (6,*) 'nnrm,nshr,ntyld:',nnrm,nshr,ntyld
          call jancae_exit (9000)
        endif
        goto 100
      endif
c
      ss = 0.0
      do i = 1,nttl
        ss = ss + cs(i)**2
      enddo
      if ( ( ss .le. 0.0 ) .and. ( nreq .eq. 0 ) ) then
        se = 0.0
        return
      endif
c
c                                      *** 3D yield functions
c
c                                  set index to s(i) to cs(i)
      do i = 1,6
        indx(i) = 0
      enddo
      if ( nnrm .eq. 3 ) then
        do i = 1,nttl
          indx(i) = i
        enddo
      else if ( nnrm .eq. 2 ) then
        indx(1) = 1
        indx(2) = 2
        indx(3) = 0
        do i = 1,nshr
          indx(3+i) = 2 + i
        enddo
      endif
c                                                    set s(i)
      call jancae_clear1 ( s,6 )
      do i = 1,6
        if ( indx(i) .ne. 0 ) then
          s(i) = cs(indx(i))
        endif
      enddo
c
      select case ( ntyld )
      case ( 0 )                           									 ! von Mises
        call jancae_mises ( s,se,dseds,d2seds2,nreq )
      case ( 1 )                     												 ! Hill 1948
        call jancae_hill_1948 ( s,se,dseds,d2seds2,nreq,
     &                          pryld,ndyld )
      case ( 2 )                      										 ! Yld2004-18p
        call jancae_yld2004_18p ( s,se,dseds,d2seds2,nreq,
     &                            pryld,ndyld )
      case ( 3 )                                              ! CPB 2006
        call jancae_cazacu2006 ( s,se,dseds,d2seds2,nreq,
     &                           pryld,ndyld )
      case ( 4 )                    				     ! Karafillis-Boyce 1993
        call jancae_KarafillisBoyce ( s,se,dseds,d2seds2,nreq,
     &                                pryld,ndyld )
      case ( 5 )                                               ! Hu 2005
        call jancae_hu_2005 ( s,se,dseds,d2seds2,nreq,
     &                       pryld,ndyld )
      case ( 6 )                                          ! Yoshida 2011
        call jancae_yoshida_2011 ( s,se,dseds,d2seds2,nreq,
     &                             pryld,ndyld )
c
      case default
        write (6,*) 'error in jancae_yfunc'
        write (6,*) 'ntyld error :',ntyld
        call jancae_exit (9000)
      end select
c
c                                                         --- set dse/ds
      if ( nreq .ge. 1 ) then
        do i = 1,6
          if ( indx(i) .ne. 0 ) cdseds(indx(i)) = dseds(i)
        enddo
      endif
c                                                       --- set d2se/ds2
      if ( nreq .ge. 2 ) then
        do i = 1,6
          if ( indx(i) .ne. 0 ) then
            do j = 1,6
              if ( indx(j) .ne. 0 ) then
                cd2seds2(indx(i),indx(j)) = d2seds2(i,j)
              endif
            enddo
          endif
        enddo
      endif
c
      return
c
c
  100 continue
c                                       *** plane stress yield functions
c
      select case ( ntyld )
      case ( -1 )                                    ! Gotoh Biquadratic
        call jancae_gotoh ( cs,se,cdseds,cd2seds2,nreq,
     &                      pryld,ndyld )
      case ( -2 )                 							            ! Yld2000-2d
        call jancae_yld2000 ( cs,se,cdseds,cd2seds2,nreq,
     &                        pryld,ndyld )
      case ( -3 )                 					                    ! Vegter
        call jancae_vegter ( cs,se,cdseds,cd2seds2,nreq,
     &                       pryld,ndyld )
      case ( -4 )                                             ! BBC 2005
        call jancae_bbc2005 ( cs,se,cdseds,cd2seds2,nreq,
     &                        pryld,ndyld )
      case ( -5 )                                                ! Yld89
        call jancae_yld89 ( cs,se,cdseds,cd2seds2,nreq,
     &                      pryld,ndyld )
      case ( -6 )                                             ! BBC 2008
        call jancae_bbc2008 ( cs,se,cdseds,cd2seds2,nreq,
     &                        pryld,ndyld )
      case ( -7 )                                            ! Hill 1990
        call jancae_hill90  ( cs,se,cdseds,cd2seds2,nreq,
     &                        pryld,ndyld )

      case default
        write (6,*) 'error in jancae_yfunc'
        write (6,*) 'ntyld error :',ntyld
        call jancae_exit (9000)
      end select
c
      return
      end
c
c
c-----------------------------------------------------------------------
c     print type and parameters for yield functions
c
      subroutine jancae_yfunc_print ( pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension pryld(ndyld)
c
      ntyld = pryld(1)
      write (6,*) '*** yield function',ntyld
      select case ( ntyld )
      case ( 0 )
        write (6,*) 'von Mises'
      case ( 1 )
        write (6,*) 'Hill 1948'
        write (6,*) 'F=',pryld(2)
        write (6,*) 'G=',pryld(3)
        write (6,*) 'H=',pryld(4)
        write (6,*) 'L=',pryld(5)
        write (6,*) 'M=',pryld(6)
        write (6,*) 'N=',pryld(7)
      case ( 2 )
        write (6,*) 'Yld2004-18p'
        n0 = 1
        do i = 1,18
           n0 = n0+1
           write (6,*) 'a(',i,')=',pryld(n0)
        enddo
        write (6,*) 'M=',pryld(1+18+1)
      case ( 3 )
        write (6,*) 'CPB 2006'
        n0 = 1
        do i = 1,3
          do j = 1,3
            n0 = n0 + 1 ; write (6,*) 'c(',i,',',j,')=',pryld(n0)
          enddo
        enddo
        do i = 4,6
          n0 = n0 + 1 ; write (6,*) 'c(',i,',',i,')=',pryld(n0)
        enddo
        n0 = n0 + 1 ; write (6,*) 'a =',pryld(n0)
        n0 = n0 + 1 ; write (6,*) 'ck=',pryld(n0)
      case ( 4 )
        write (6,*) 'Karafillis-Boyce 1993'
        n0 = 1
        do i = 1,6
          do j = i,6
            n0 = n0 + 1
            write (6,*) 'L(',i,',',j,') =',pryld(n0)
          end do
        end do
        n0 = n0 + 1 ; write (6,*) 'k =',pryld(n0)
        n0 = n0 + 1 ; write (6,*) 'c =',pryld(n0)
      case ( 5 )
        write (6,*) 'Hu 2005'
        n0 = 1
        do i = 1,5
          n0 = n0 + 1
          write (6,*) 'X(',i,')=',pryld(n0)
        enddo
        n0 = n0 + 1 ; write (6,*) 'X(',7,')=',pryld(n0)
        do i = 1,3
          n0 = n0 + 1
          write (6,*) 'C(',i,')=',pryld(n0)
        enddo
      case ( 6 )
        write (6,*) 'Yoshida 2011'
        n0 = 1
        do i = 1,16
          n0 = n0 + 1
          write (6,*) 'c(',i,')=',pryld(n0)
        enddo
c
      case ( -1 )
        write (6,*) 'Gotoh Biquadratic'
        do i = 1,9
          write (6,*) 'A(',i,')=',pryld(i+1)
        enddo
      case ( -2 )
        write (6,*) 'Yld2000-2d'
        do i = 1,8
          write (6,*) 'a(',i,')=',pryld(i+1)
        enddo
        write (6,*) 'M=',pryld(9+1)
      case ( -3 )
        write (6,*) 'Vegter'
        write (6,*) 'nf=',nint(pryld(2))
        write (6,*) 'f_bi0=',pryld(3)
        write (6,*) 'r_bi0=',pryld(4)
        do i = 0,nint(pryld(2))
          write (6,*) 'test angle=',90.0d0*float(i)/pryld(2)
          write (6,*) 'phi_un(',i,')=',pryld(4+i*4+1)
          write (6,*) 'phi_sh(',i,')=',pryld(4+i*4+2)
          write (6,*) 'phi_ps(',i,')=',pryld(4+i*4+3)
          write (6,*) 'omg(   ',i,')=',pryld(4+i*4+4)
        enddo
c       do i=1,7
c         write (6,*) 'phi_un(',i-1,')=',pryld(1+i   )
c         write (6,*) 'phi_sh(',i-1,')=',pryld(1+i+ 7)
c         write (6,*) 'phi_ps(',i-1,')=',pryld(1+i+14)
c         write (6,*) 'omg   (',i-1,')=',pryld(1+i+23)
c       enddo
c       write (6,*)   'f_bi0=',pryld(1+22)
c       write (6,*)   'r_bi0=',pryld(1+23)
c       write (6,*)   'nf   =',nint(pryld(1+31))
      case ( -4 )
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
      case ( -5 )
        write (6,*) 'Yld89'
        write (6,*) 'order M=',pryld(1+1)
        write (6,*) 'a      =',pryld(1+2)
        write (6,*) 'h      =',pryld(1+3)
        write (6,*) 'p      =',pryld(1+4)
      case ( -6 )
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
        enddo
      case ( -7 )
        write (6,*) 'Hill 1990'
        write (6,*) 'a   =',pryld(1+1)
        write (6,*) 'b   =',pryld(1+2)
        write (6,*) 'tau =',pryld(1+3)
        write (6,*) 'sigb=',pryld(1+4)
        write (6,*) 'M   =',pryld(1+5)
      end select
c
      return
      end
c
c

************************************************************************
*
*     YIELD FUNCTIONS
*
************************************************************************
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     YIELD FUNCTION
c
      subroutine ummdp_yield ( se,cdseds,cd2seds2,nreq,cs,nttl,nnrm,
     1                         nshr,pryld,ndyld )
c
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
          write (6,*) 'error in ummdp_yield'
          write (6,*) 'ntyld<0 for plane stress'
          write (6,*) 'nnrm,nshr,ntyld:',nnrm,nshr,ntyld
          call ummdp_exit ( 304 )
        end if
        goto 100
      end if
c
      ss = 0.0d0
      do i = 1,nttl
        ss = ss + cs(i)**2
      end do
      if ( (ss <= 0.0d0) .and. (nreq == 0) ) then
        se = 0.0d0
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
      s = 0.0d0
      do i = 1,6
        if ( indx(i) /= 0 ) then
          s(i) = cs(indx(i))
        end if
      end do
c
      select case ( ntyld )
      case ( 0 )
        call ummdp_yield_mises ( s,se,dseds,d2seds2,nreq )
      case ( 1 )
        call ummdp_yield_hill1948 ( s,se,dseds,d2seds2,nreq,pryld,
     1                              ndyld )
      case ( 2 )
        call ummdp_yield_yld2004 ( s,se,dseds,d2seds2,nreq,pryld,ndyld )
     1
      case ( 3 )
        call ummdp_yield_cpb2006 ( s,se,dseds,d2seds2,nreq,pryld,ndyld )
      case ( 4 )
        call ummdp_yield_karafillisboyce ( s,se,dseds,d2seds2,nreq,
     1                                     pryld,ndyld )
      case ( 5 )
        call ummdp_yield_hu2005 ( s,se,dseds,d2seds2,nreq,pryld,ndyld )
      case ( 6 )
        call ummdp_yield_yoshida2011 ( s,se,dseds,d2seds2,nreq,pryld,
     1                                 ndyld )
c
      case default
        write (6,*) 'error in ummdp_yield'
        write (6,*) 'ntyld error :',ntyld
        call ummdp_exit ( 202 )
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
      case ( -1 )
        call ummdp_yield_gotoh ( cs,se,cdseds,cd2seds2,nreq,pryld,
     1                           ndyld )
      case ( -2 )
        call ummdp_yield_yld2000 ( cs,se,cdseds,cd2seds2,nreq,pryld,
     1                             ndyld )
      case ( -3 )
        call ummdp_yield_vegter ( cs,se,cdseds,cd2seds2,nreq,pryld,
     1                            ndyld )
      case ( -4 )
        call ummdp_yield_bbc2005 ( cs,se,cdseds,cd2seds2,nreq,pryld,
     1                             ndyld )
      case ( -5 )
        call ummdp_yield_yld89 ( cs,se,cdseds,cd2seds2,nreq, pryld,
     1                           ndyld )
      case ( -6 )
        call ummdp_yield_bbc2008 ( cs,se,cdseds,cd2seds2,nreq,pryld,
     1                             ndyld )
      case ( -7 )
        call ummdp_yield_hill1990 ( cs,se,cdseds,cd2seds2,nreq,pryld,
     1                              ndyld )
      case default
        write (6,*) 'error in ummdp_yield'
        write (6,*) 'ntyld error :',ntyld
        call ummdp_exit ( 202 )
      end select
c
      return
      end subroutine ummdp_yield
c
c
c
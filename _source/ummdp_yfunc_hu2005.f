c---------------------------------------------------------------(hu2005)
c     WeiLong Hu (2005) yield function and its dfferentials
c
c     "An orthotropic yield criterion in a 3-D general stress state"
c      IJP,v.21(2005), pp.1771-1796. )
c
      subroutine jancae_hu_2005 ( s,se,dseds,d2seds2,nreq,
     &                            pryld,ndyld )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (maxa=100)
      dimension s(6),dseds(6),d2seds2(6,6),pryld(ndyld)
      dimension a(maxa),ipow(maxa,3)
c
c
      nd0=2
c
      n=0
      do it=0,nd0
        n=n+(nd0-it)*2+1
      enddo
      nterms=n
      if ( maxa.lt.nterms ) then
        write (6,*) 'increase maxa :',maxa,nterms
        call jancae_exit ( 9000 )
      endif
c
      n=0
      ipow=0
      do it=0,nd0
        ndmax=nd0*2-it*2
        do jy=0,ndmax
          jx=ndmax-jy
          n=n+1
          ipow(n,1)=jx
          ipow(n,2)=jy
          ipow(n,3)=it
        enddo
      enddo
c
      a     = 0.0
      a( 1)= 1.0d0*pryld(1+ 1)    ! X1 
      a( 2)= 1.0d0*pryld(1+ 2)    ! X2
      a( 3)= 1.0d0*pryld(1+ 3)    ! X3
      a( 4)= 1.0d0*pryld(1+ 4)    ! X4
      a( 5)= 1.0d0*pryld(1+ 5)    ! X5
      a( 6)= 1.0d0*pryld(1+ 7)    ! C1 <-
      a( 7)=-1.0d0*pryld(1+ 9)    ! C3 <-
      a( 8)= 1.0d0*pryld(1+ 8)    ! C2 <-
      a( 9)= 1.0d0*pryld(1+ 6)    ! X7
c
      call jancae_hy_polytype ( s,se,dseds,d2seds2,nreq,
     &                          nd0,a,ipow,maxa,nterms )
c
      return
      end
c
c
c
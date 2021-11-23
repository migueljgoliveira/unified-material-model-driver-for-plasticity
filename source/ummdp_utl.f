c***********************************************************************
c
c     UMMDp: Utility Subroutines
c
c***********************************************************************
c
c     jancae_clear1 ( a,n )
c       clear 1st order vector
c
c     jancae_clear2 ( a,n,m )
c       clear 2nd order matrix
c
c     jancae_clear3 ( a,n,m,l )
c       clear 3rd order tensor
c
c     jancae_setunitm ( a,n )
c       set unit 2nd order matrix
c
c     jancae_print1 ( text,a,n )
c       print vector with text
c
c     jancae_print2 ( text,a,n,m )
c       print matrix with text
c
c     jancae_mv (v,a,u,nv,nu)
c       mutiply matrix and vector
c
c     jancae_mm (a,b,c,na1,na2,nbc)
c       mutiply matrix and matrix
c
c     jancae_vvs ( s,u,v,n )
c       calculate scalar product of vectors 
c
c     jancae_minv ( b,a,n,d )
c       calculate inverse matrix using lu decomposition
c
c         jancae_ludcmp( a,n,indx,d )
c           lu decomposition
c         jancae_lubksb(a,n,indx,b)
c           lu backward substitution
c         jancae_minv2 ( b,a,deta )
c           calculate inverse matrix 2x2
c         jancae_minv3 ( b,a,deta )
c           calculate inverse matrix 3x3
c
c     jancae_eigen_sym3 ( es,ev,a )
c       calculate eigenvalues and eigenvectors by jacobi method
c
c     jancae_printinfo
c       print informations for debug (info)
c
c     jancae_printinout
c       print informations for debug (input/output)
c
c
c
************************************************************************
c     CLEAR 1st ORDER VECTOR A(N)
c
      subroutine jancae_clear1 ( a,n )
c-----------------------------------------------------------------------
      implicit none
c
      integer n
      real*8 a(n)
c
      integer i
c-----------------------------------------------------------------------
c
      do i = 1,n
        a(i) = 0.0d0
      end do
c
      return
      end subroutine jancae_clear1
c
c
c
************************************************************************
c     CLEAR 2ND ORDER MATRIX
c
      subroutine jancae_clear2 ( a,n,m )
c-----------------------------------------------------------------------
      implicit none
c
      integer n,m
      real*8 a(n,m)
c
      integer i,j
c-----------------------------------------------------------------------
c
      do i = 1,n
        do j = 1,m
          a(i,j) = 0.0d0
        end do
      end do
c
      return
      end subroutine jancae_clear2
c
c
c
************************************************************************
c     CLEAR 3RD ORDER MATRIX
c
      subroutine jancae_clear3 ( a,n,m,l )
c-----------------------------------------------------------------------
      implicit none
c
      integer n,m,l
      real*8 a(n,m,l)
c
      integer i,j,k
c-----------------------------------------------------------------------
c
      do i = 1,n
        do j = 1,m
          do k = 1,l
            a(i,j,k) = 0.0d0
          end do
        end do
      end do
c
      return
      end subroutine jancae_clear3
c
c
c
************************************************************************
c     SET UNIT 2ND ORDER MATRIX
c
      subroutine jancae_setunitm ( a,n )
c-----------------------------------------------------------------------
      implicit none
c
      integer n
      real*8 a(n,n)
c
      integer i
c-----------------------------------------------------------------------
c
      call jancae_clear2 ( a,n,n )
      do i = 1,n
        a(i,i) = 1.0d0
      end do
c
      return
      end subroutine jancae_setunitm
c
c
c
************************************************************************
c     PRINT VECTOR WITH TEXT
c
      subroutine jancae_print1 ( text,a,n )
c-----------------------------------------------------------------------
      implicit none
c
      integer n
      real*8 a(n)
      character*32 text
c
      integer i
c-----------------------------------------------------------------------
c
      write (6,*) text
      write (6,9000) (a(i),i=1,n)
 9000 format (6e16.8)
c
      return
      end subroutine jancae_print1
c
c
c
************************************************************************
c     PRINT MATRIX WITH TEXT
c
      subroutine jancae_print2 ( text,a,n,m )
c-----------------------------------------------------------------------
      implicit none
c
      integer n,m
      real*8 a(n,m)
      character*32 text
c
      integer i,j
c-----------------------------------------------------------------------
      write (6,*) text
      do i = 1,n
        write (6,9000) (a(i,j),j=1,m)
      end do
 9000 format (6e16.8)
c
      return
      end subroutine jancae_print2
c
c
c
************************************************************************
c     MULTIPLY MATRIX AND VECTOR
c
      subroutine jancae_mv ( v,a,u,nv,nu )
c-----------------------------------------------------------------------
      implicit none
c
      integer nv,nu
      real*8 v(nv),u(nu)
      real*8 a(nv,nu)
c
      integer i,j
c-----------------------------------------------------------------------
c
      call jancae_clear1 ( v,nv )
      do i = 1,nv
        do j = 1,nu
          v(i) = v(i) + a(i,j)*u(j)
        end do
      end do
c
      return
      end subroutine jancae_mv
c
c
c
************************************************************************
c     MULTIPLY MATRIX AND MATRIX
c     
      subroutine jancae_mm ( a,b,c,na1,na2,nbc )
c-----------------------------------------------------------------------
      implicit none
c
      integer na1,na2,nbc
      real*8 a(na1,na2),b(na1,nbc),c(nbc,na2)
c
      integer i,j,k
c-----------------------------------------------------------------------
c
      call jancae_clear2 ( a,na1,na2 )
      do i = 1,na1
        do j = 1,na2
          do k = 1,nbc
            a(i,j) = a(i,j) + b(i,k)*c(k,j)
          end do
        end do
      end do
c
      return
      end subroutine jancae_mm
c
c
c
************************************************************************
c     CALCULATE SCALAR PRODUCT OF VECTORS
c
      subroutine jancae_vvs ( s,u,v,n )
c-----------------------------------------------------------------------
      implicit none
c
      integer n
      real*8 s
      real*8 v(n),u(n)
c
      integer i
c-----------------------------------------------------------------------
c
      s = 0.0
      do i = 1,n
        s = s + u(i)*v(i)
      end do
c
      return
      end subroutine jancae_vvs
c
c
c
************************************************************************
c     CALCULATE INVERSE MATRIX USING LU DECOMPOSITION
c
c     Ref.: http://astr-www.kj.yamagata-u.ac.jp/~shibata/kbg/
c
      subroutine jancae_minv ( b,a,n,d )
c-----------------------------------------------------------------------
      implicit none
c
      integer n
      real*8 d
      real*8 a(n,n),b(n,n)
c
      integer i,j,k
      real*8 eps,anorm,ani
      real*8 indx(n),y(n),c(n,n),aorg(n,n)
      character*32 text
      logical check
c-----------------------------------------------------------------------
c
      check = .false.
      eps = 1.0d-36
c
      do i = 1,n
        do j = 1,n
          aorg(i,j) = a(i,j)
        end do
      end do
c
      anorm = 0.0
      do i = 1,n
        do j = 1,n
          if ( anorm < abs(a(i,j)) ) anorm = abs(a(i,j))
        end do
      end do
      do i = 1,n
        do j = 1,n
          a(i,j) = a(i,j) / anorm
        end do
      end do
c
      if ( n == 2 ) then
        call jancae_minv2 ( b,a,d,eps )
        goto 100
      else if ( n == 3 ) then
        call jancae_minv3 ( b,a,d,eps )
        goto 100
      end if
c
      call jancae_ludcmp ( a,n,indx,d,eps )
c                                                 ---- check determinant
      if ( abs(d) <= eps ) then
         write (6,*) 'determinant det[a] error',d
         write (6,*) 'stop in minv'
         call jancae_exit ( 9000 ) 
      end if
c                                                            ---- B=A^-1
      do j = 1,n
        call jancae_clear1 ( y,n )
        y(j) = 1.0d0
        call jancae_lubksb ( a,n,indx,y,eps )
        do i = 1,n
          b(i,j) = y(i)
        end do
      end do
c
  100 continue
      ani = 1.0d0/anorm
      do i = 1,n
        do j = 1,n
          a(i,j) = aorg(i,j)
          b(i,j) = b(i,j) * ani
        end do
      end do
c                                                             ---- check
      if ( check ) then
        write (6,*) 'check inverse matrix',n
        text = 'original matrix [A]'
        call jancae_print2 ( text,a,n,n )
        text = 'inversed matrix [A]^-1'
        call jancae_print2 ( text,b,n,n )
        call jancae_clear2 ( c,n,n )
        do i = 1,n
          do j = 1,n
            do k = 1,n
              c(i,j) = c(i,j) + b(i,k)*a(k,j)
            end do
          end do
        end do
        text = '[A]^-1*[A]=[I] ?'
        call jancae_print2 ( text,c,n,n )
      end if
c
      return
      end subroutine jancae_minv
c
c
c
************************************************************************
c     LU DECOMPOSITION
c
      subroutine jancae_ludcmp ( a,n,indx,d,eps )
c
c-----------------------------------------------------------------------
      implicit none
c
			integer n
			integer indx(n)
			real*8 d,eps
      real*8 a(n,n)
c
			integer i,j,k,imax
			real*8 aamax,sum,dum,ajj
			real*8 vtemp(n)
      character*32 text
c-----------------------------------------------------------------------
c
      d = 1.0d0
      do i = 1,n
        aamax = 0.0d0
        do j = 1,n
          if ( abs(a(i,j)) > aamax ) aamax = abs(a(i,j))
        end do
        if ( aamax <= eps ) then
          write (6,*) 'singular matrix in jancae_ludcmp'
          text = 'matrix detail'
          call jancae_print2 ( text,a,n,n )
          call jancae_exit ( 9000 )
        end if
        vtemp(i) = 1.0d0 / aamax
      end do
c
      do j = 1,n
        do i = 1,j-1
          sum = a(i,j)
          do k = 1,i-1
            sum = sum - a(i,k)*a(k,j)
          end do
          a(i,j) = sum
        end do
        aamax = 0.0d0
        do i = j,n
          sum = a(i,j)
          do k = 1,j-1
            sum = sum - a(i,k)*a(k,j)
          end do
          a(i,j) = sum
          dum = vtemp(i)*abs(sum)
          if ( dum >= aamax ) then
            imax = i
            aamax = dum
          end if
        end do
        if ( j /= imax ) then
          do k = 1,n
            dum = a(imax,k)
            a(imax,k) = a(j,k)
            a(j,k) = dum
          end do
          d = -d
          vtemp(imax) = vtemp(j)
        end if
        indx(j) = imax
c       if ( abs(a(i,j)) <= eps ) a(i,j) = eps     !2010.07.02 c.out
        if ( j /= n ) then
          ajj = a(j,j)                               !2010.07.02 add
          if ( abs(ajj) <= eps ) ajj = eps         !2010.07.02 add
          dum = 1.0d0 / ajj                          !2010.07.02 mod
          do i = j+1,n
            a(i,j) = a(i,j) * dum
          end do
        end if
      end do
c                                                 ---- get the det. of A
      do j = 1,n
        d = d*a(j,j)
      end do
c
      return
      end subroutine jancae_ludcmp
c
c
c
************************************************************************
c     LU BACKWARD SUBSTITUTION
c
      subroutine jancae_lubksb ( a,n,indx,b,eps )
c-----------------------------------------------------------------------
      implicit none
c
			integer n
			integer indx(n)
			real*8 eps
			real*8 a(n,n),b(n)
c
			integer i,j,ii,ll
			real*8 sum
c-----------------------------------------------------------------------
c
      ii = 0
      do i = 1,n
        ll = indx(i)
        sum = b(ll)
        b(ll) = b(i)
        if ( ii /= 0 ) then
          do j = ii,i-1
            sum = sum - a(i,j)*b(j)
          end do
        else if ( abs(sum) >= eps ) then
          ii = i
        end if
        b(i) = sum
      end do
      do i = n,1,-1
        sum = b(i)
        do j = i+1,n
          sum = sum - a(i,j)*b(j)
        end do
        b(i) = sum / a(i,i)
      end do
c
      return
      end subroutine jancae_lubksb
c
c
c
************************************************************************
c     CALCULATE INVERSE MATRIX 2x2 
c
      subroutine jancae_minv2 ( b,a,deta,eps )
c-----------------------------------------------------------------------
      implicit none
c
			real*8 deta,eps
			real*8 b(2,2),a(2,2)
c
			real*8 detai
c-----------------------------------------------------------------------
c
      deta = a(1,1)*a(2,2) - a(1,2)*a(2,1)
      if ( abs(deta) <= eps ) then
         write (6,*) 'determinant det[a] error',deta
         write (6,*) 'stop in minv2'
         call jancae_exit ( 9000 )
      end if
c
      detai = 1.0d0 / deta
      b(1,1) = a(2,2) * detai
      b(1,2) = -a(1,2) * detai
      b(2,1) = -a(2,1) * detai
      b(2,2) = a(1,1) * detai
c
      return
      end subroutine jancae_minv2
c
c
c
************************************************************************
c     CALCULATE INVERSE MATRIX 3x3
c
      subroutine jancae_minv3 ( b,a,deta,eps )
c-----------------------------------------------------------------------
      implicit none
c
			real*8 deta,eps
			real*8 b(3,3),a(3,3)
c
			real*8 detai
c-----------------------------------------------------------------------
c
      deta = a(1,1) * (a(2,2)*a(3,3) - a(2,3)*a(3,2)) +
     &       a(1,2) * (a(2,3)*a(3,1) - a(2,1)*a(3,3)) +
     &       a(1,3) * (a(2,1)*a(3,2) - a(2,2)*a(3,1))
      if ( abs(deta) <= eps ) then
         write (6,*) 'determinant det[a] error',deta
         write (6,*) 'stop in minv3'
         call jancae_exit ( 9000 )
      end if
c
      detai = 1.0d0 / deta
      b(1,1) = ( a(2,2)*a(3,3) - a(2,3)*a(3,2) ) * detai
      b(1,2) = ( a(1,3)*a(3,2) - a(1,2)*a(3,3) ) * detai
      b(1,3) = ( a(1,2)*a(2,3) - a(1,3)*a(2,2) ) * detai
      b(2,1) = ( a(2,3)*a(3,1) - a(2,1)*a(3,3) ) * detai
      b(2,2) = ( a(1,1)*a(3,3) - a(1,3)*a(3,1) ) * detai
      b(2,3) = ( a(1,3)*a(2,1) - a(1,1)*a(2,3) ) * detai
      b(3,1) = ( a(2,1)*a(3,2) - a(2,2)*a(3,1) ) * detai
      b(3,2) = ( a(1,2)*a(3,1) - a(1,1)*a(3,2) ) * detai
      b(3,3) = ( a(1,1)*a(2,2) - a(1,2)*a(2,1) ) * detai
c
      return
      end subroutine jancae_minv3
c
c
c
************************************************************************
c     PRINT INFORMATIONS FOR DEBUG (INFO)
c
      subroutine jancae_printinfo ( inc,nnrm,nshr )
c-----------------------------------------------------------------------
      implicit none
c
			integer inc,nnrm,nshr
c
      integer nttl,nerr
c
			integer ne,ip,lay
      common /jancae1/ne,ip,lay
c-----------------------------------------------------------------------
c
      nttl = nnrm + nshr
c
      write (6,*) '----- JANCAE.UMMDp Debug Info -----'
      write (6,*) 'increment=',inc
      write (6,*) 'elem,ip,lay=',ne,ip,lay
      write (6,*) 'nttl,nnrm,nshr=',nttl,nnrm,nshr
      nerr = 0
      if ( nnrm == 3 ) then
        if ( nshr == 3 ) then
          write (6,*) '3d solid element'
        else if ( nshr == 1 ) then
          write (6,*) 'plane strain or axi-sym solid element'
        else
          nerr = nerr + 1
        end if
      else if ( nnrm == 2 ) then
        if ( nshr == 1 ) then
          write (6,*) 'plane stress or thin shell element'
        else if ( nshr == 3 ) then
          write (6,*) 'thick shell element'
        else
          nerr = nerr + 1
        end if
      else
        nerr = nerr + 1
      end if
      if ( nerr /= 0 ) then
        write (6,*) 'no supported element type',nnrm,nshr
        call jancae_exit ( 9000 )
      end if
c
      return
      end subroutine jancae_printinfo
c
c
c
************************************************************************
c     PRINT INFORMATIONS FOR DEBUG (INPUT/OUTPUT)
c
      subroutine jancae_printinout ( io,s,de,d,nttl,stv,nstv )
c-----------------------------------------------------------------------
      implicit none
c
			integer io,nttl,nstv
			real*8 s(nttl),stv(nstv),de(nttl),d(nttl,nttl)
c
      character*32 text
c-----------------------------------------------------------------------
c
      if ( io == 0 ) then
        text = 'initial stresses'
      else
        text = 'updated stresses'
      end if
      call jancae_print1 ( text,s,nttl )
c
      if ( io == 0 ) then
        text = 'initial internal state var.'
      else
        text = 'updated internal state var.'
      end if
      call jancae_print1 ( text,stv,nstv )
c
      if ( io == 0 ) then
        text = 'driving strain increment'
        call jancae_print1 ( text,de,nttl )
      else
        text = 'tangent modulus matrix'
        call jancae_print2 ( text,d,nttl,nttl )
      end if
c
      return
      end subroutine jancae_printinout
c
c
c
************************************************************************
c     CALCULATE EIGENVALUES AND EIGENVECTORS BY JACOBI METHOD
c
c     Ref.: http://www.flagshyp.com/
c
c     input
c       a(3,3)  : symmetric matrix to be analyzed
c
c     output
c       es(i)   : i-th eigenvalue
c       ev(i,3) : normalized eigenvector for i-th eigenvalue
c
      subroutine jancae_eigen_sym3 ( es,ev,a )
c-----------------------------------------------------------------------
      implicit none
c
      real*8 es(3)
			real*8 ev(3,3),a(3,3)
c
			integer msweep,i,j,is,ip,ir,iq
			real*8 eps,ax,er,sum,od,hd,theta,t,c,s,tau,h,g
      real*8 w(3,3),prc(3,3)
c-----------------------------------------------------------------------
c
      msweep = 100
      eps = 1.0d-8
c
c                                                       ---- preparation
      ax = 0.0d0
      er = 0.0d0
      do i = 1,3
        do j = 1,3
          if ( abs(a(i,j)) > ax ) ax = abs(a(i,j))
          er = er + abs(a(i,j)-a(j,i))
        end do
      end do
      if ( er/ax > eps ) then
        write (6,*) 'a is not symmetric'
        write (6,*) 'stop in jancae_eigen_sym3'
        call jancae_exit ( 9000 )
      end if
      do i = 1,3
        do j = 1,3
          w(i,j) = a(i,j) / ax
        end do
      end do
c                                    ---- initialise prc to the identity
      do i = 1,3
        do j = 1,3
          prc(i,j) = 0.0d0
        end do
        prc(i,i) = 1.0d0
        es(i) = w(i,i)
      end do
c                                                   ---- starts sweeping
      do is = 1,msweep
c
        sum = 0.0d0
        do ip = 1,2
          do iq = ip+1,3
            sum = sum + abs( w(ip,iq) )
          end do
        end do
c       write (6,*) 'ite',is,sum,eps
c            ---- if the sum of off-diagonal terms is zero evaluates the
c                                                     esches and returns
        if ( abs(sum) < eps ) then
          do i = 1,3
            do j = 1,3
              ev(i,j) = prc(j,i)
            end do
            es(i) = es(i)*ax
          end do
          return
        end if
c                             ---- performs the sweep in three rotations
c                                         ---- one per off diagonal term
        do ip = 1,2
          do iq = ip+1,3
            od = 100.0d0 * abs( w(ip,iq) )
            if ( abs(od) > eps ) then
              hd = es(iq) - es(ip)
c                                      ---- evaluates the rotation angle
              theta = 0.5d0 * hd / w(ip,iq)
              t = 1.0d0/(abs(theta) + sqrt(1.0d0+theta**2))
              if ( theta < 0.0d0 ) t = -t
c                                   ---- re-evaluates the diagonal terms
              c = 1.0d0 / sqrt(1.0d0+t**2)
              s = t * c
              tau = s / (1.0d0+c)
              h = t * w(ip,iq)
              es(ip) = es(ip) - h
              es(iq) = es(iq) + h
c                     ---- re-evaluates the remaining off-diagonal terms
              ir = 6 - ip - iq
              g = w( min(ir,ip),max(ir,ip) )
              h = w( min(ir,iq),max(ir,iq) )
              w( min(ir,ip),max(ir,ip) ) = g - s*(h+g*tau)
              w( min(ir,iq),max(ir,iq) ) = h + s*(g-h*tau)
c                                          ---- rotates the eigenvectors
              do ir = 1,3
                g = prc(ir,ip)
                h = prc(ir,iq)
                prc(ir,ip) = g - s*(h+g*tau)
                prc(ir,iq) = h + s*(g-h*tau)
              end do
            end if
            w(ip,iq) = 0.0d0
          end do
        end do
      end do
c                              ---- if convergence is not achieved stops
      write (6,*) 'did not converge in eigen calculation.'
      write (6,*) 'msweep=',msweep
      write (6,*) 'eps=',eps
      write (6,*) 'sum=',sum
      write (6,*) 'stop in jancae_eigen_sym3'
      call jancae_exit ( 9000 )
c
      return
      end subroutine jancae_eigen_sym3
c
c
c
************************************************************************
c     CHECKING EXISTENCE OF FILE NAMES 'FLNAME'
c
      logical function jancae_file_exist ( flname )
c
c-----------------------------------------------------------------------
      implicit none
c
      character*16 flname
c
			integer nio
c-----------------------------------------------------------------------
c
      nio = 616
      open ( nio,file=flname,status='old',err=10 )
c
      close ( nio,            status='keep' )
      jancae_file_exist = .true.
      return
c
   10 jancae_file_exist = .false.
      return
c
      end function jancae_file_exist
c
c
c

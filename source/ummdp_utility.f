************************************************************************
*
*     UTILITY SUBROUTINES
*
************************************************************************
c
c     ummdp_utility_clear1 ( a,n )
c       clear 1st order vector
c
c     ummdp_utility_clear2 ( a,n,m )
c       clear 2nd order matrix
c
c     ummdp_utility_clear3 ( a,n,m,l )
c       clear 3rd order tensor
c
c     ummdp_utility_setunitm ( a,n )
c       set unit 2nd order matrix
c
c     ummdp_utility_print1 ( text,a,n )
c       print vector with text
c
c     ummdp_utility_print2 ( text,a,n,m )
c       print matrix with text
c
c     ummdp_utility_mv (v,a,u,nv,nu)
c       mutiply matrix and vector
c
c     ummdp_utility_mm (a,b,c,na1,na2,nbc)
c       mutiply matrix and matrix
c
c     ummdp_utility_vvs ( s,u,v,n )
c       calculate scalar product of vectors 
c
c     ummdp_utility_minv ( b,a,n,d )
c       calculate inverse matrix using lu decomposition
c
c         ummdp_utility_ludcmp( a,n,indx,d )
c           lu decomposition
c         ummdp_utility_lubksb(a,n,indx,b)
c           lu backward substitution
c         ummdp_utility_minv2 ( b,a,deta )
c           calculate inverse matrix 2x2
c         ummdp_utility_minv3 ( b,a,deta )
c           calculate inverse matrix 3x3
c
c     ummdp_utility_eigen_sym3 ( es,ev,a )
c       calculate eigenvalues and eigenvectors by jacobi method
c
c     ummdp_utility_file_exist ( flname )
c       checking existence of files
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     CLEAR 1st ORDER VECTOR A(N)
c
      subroutine ummdp_utility_clear1 ( a,n )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: n
c
      real*8,intent(inout) :: a(n)
c
      integer i
c-----------------------------------------------------------------------
c
      do i = 1,n
        a(i) = 0.0d0
      end do
c
      return
      end subroutine ummdp_utility_clear1
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     CLEAR 2ND ORDER MATRIX
c
      subroutine ummdp_utility_clear2 ( a,n,m )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: n,m
c
      real*8,intent(inout) :: a(n,m)
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
      end subroutine ummdp_utility_clear2
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     CLEAR 3RD ORDER MATRIX
c
      subroutine ummdp_utility_clear3 ( a,n,m,l )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: n,m,l
c
      real*8,intent(inout) ::  a(n,m,l)
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
      end subroutine ummdp_utility_clear3
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     SET UNIT 2ND ORDER MATRIX
c
      subroutine ummdp_utility_setunitm ( a,n )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: n
c
      real*8,intent(out) :: a(n,n)
c
      integer i
c-----------------------------------------------------------------------
c
      call ummdp_utility_clear2 ( a,n,n )
      do i = 1,n
        a(i,i) = 1.0d0
      end do
c
      return
      end subroutine ummdp_utility_setunitm
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     PRINT VECTOR WITH TEXT
c
      subroutine ummdp_utility_print1 ( text,a,n,tab )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer      ,intent(in) :: n,tab
      real*8       ,intent(in) :: a(n)
      character*100,intent(in) :: text
c
      integer i
      character*20 fmt
c-----------------------------------------------------------------------
c
      write(fmt,'(A,I2,A)') '(/',12+tab,'xA,A)'
      write (6,fmt) '. ',text
c
      write(fmt,'(A,I2,A)') '(',14+tab,'x6E20.12)'
      ! write (6,'(14x6E16.8)') (a(i),i=1,n)
      write (6,fmt) (a(i),i=1,n)
      
c
      return
      end subroutine ummdp_utility_print1
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     PRINT MATRIX WITH TEXT
c
      subroutine ummdp_utility_print2 ( text,a,n,m,tab )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer      ,intent(in) :: n,m,tab
      real*8       ,intent(in) :: a(n,m)
      character*100,intent(in) :: text
c
      integer i,j
      character*20 fmt
c-----------------------------------------------------------------------
c
      write(fmt,'(A,I2,A)') '(/',12+tab,'xA,A)'
      write (6,fmt) '. ',text
c
      write(fmt,'(A,I2,A)') '(',14+tab,'x6E20.12)'
      do i = 1,n
        ! write (6,'(14x6E16.8)') (a(i,j),j=1,m)
        write (6,fmt) (a(i,j),j=1,m)
      end do
c
      return
      end subroutine ummdp_utility_print2
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     PRINT REAL WITH TEXT
c
      subroutine ummdp_utility_print3 ( text,a,tab )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer      ,intent(in) :: tab
      real*8       ,intent(in) :: a
      character*100,intent(in) :: text
c
      character*20 fmt
c-----------------------------------------------------------------------
c
      write(fmt,'(A,I2,A)') '(/',12+tab,'xA,A)'
      write (6,fmt) '. ',text
c
      write(fmt,'(A,I2,A)') '(',14+tab,'x6E20.12)'
      ! write (6,'(14x6E16.8)') a
      write (6,fmt) a
c
      return
      end subroutine ummdp_utility_print3
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     MULTIPLY MATRIX AND VECTOR
c
      subroutine ummdp_utility_mv ( v,a,u,nv,nu )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: nv,nu
      real*8 ,intent(in) :: u(nu)
      real*8 ,intent(in) :: a(nv,nu)
c
      real*8,intent(out) :: v(nv)
c
      integer i,j
c-----------------------------------------------------------------------
c
      call ummdp_utility_clear1 ( v,nv )
      do i = 1,nv
        do j = 1,nu
          v(i) = v(i) + a(i,j)*u(j)
        end do
      end do
c
      return
      end subroutine ummdp_utility_mv
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     MATRIX PRODUCT OF TWO MATRICES
c     
      subroutine ummdp_utility_mm ( a,b,c,na1,na2,nbc )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: na1,na2,nbc
      real*8 ,intent(in) :: b(na1,nbc),c(nbc,na2)
c
      real*8,intent(out) :: a(na1,na2)
c
      integer i,j,k
c-----------------------------------------------------------------------
c
      call ummdp_utility_clear2 ( a,na1,na2 )
      do i = 1,na1
        do j = 1,na2
          do k = 1,nbc
            a(i,j) = a(i,j) + b(i,k)*c(k,j)
          end do
        end do
      end do
c
      return
      end subroutine ummdp_utility_mm
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     SCALAR PRODUCT OF VECTORS
c
      subroutine ummdp_utility_vvs ( s,u,v,n )
c
c-----------------------------------------------------------------------
      implicit none
c
      integer,intent(in) :: n
      real*8 ,intent(in) :: v(n),u(n)
c
      real*8,intent(out) :: s
c
      integer i
c-----------------------------------------------------------------------
c
      s = 0.0d0
      do i = 1,n
        s = s + u(i)*v(i)
      end do
c
      return
      end subroutine ummdp_utility_vvs
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     CALCULATE INVERSE MATRIX USING LU DECOMPOSITION
c
c       Ref.: http://astr-www.kj.yamagata-u.ac.jp/~shibata/kbg/
c
      subroutine ummdp_utility_minv ( b,a,n,d )
c
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
        call ummdp_utility_minv2 ( b,a,d,eps )
        goto 100
      else if ( n == 3 ) then
        call ummdp_utility_minv3 ( b,a,d,eps )
        goto 100
      end if
c
      call ummdp_utility_ludcmp ( a,n,indx,d,eps )
c                                                 ---- check determinant
      if ( abs(d) <= eps ) then
         write (6,*) 'determinant det[a] error',d
         write (6,*) 'stop in minv'
         call ummdp_exit ( 9000 ) 
      end if
c                                                            ---- B=A^-1
      do j = 1,n
        call ummdp_utility_clear1 ( y,n )
        y(j) = 1.0d0
        call ummdp_utility_lubksb ( a,n,indx,y,eps )
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
        call ummdp_utility_print2 ( text,a,n,n )
        text = 'inversed matrix [A]^-1'
        call ummdp_utility_print2 ( text,b,n,n )
        call ummdp_utility_clear2 ( c,n,n )
        do i = 1,n
          do j = 1,n
            do k = 1,n
              c(i,j) = c(i,j) + b(i,k)*a(k,j)
            end do
          end do
        end do
        text = '[A]^-1*[A]=[I] ?'
        call ummdp_utility_print2 ( text,c,n,n )
      end if
c
      return
      end subroutine ummdp_utility_minv
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     LU DECOMPOSITION
c
      subroutine ummdp_utility_ludcmp ( a,n,indx,d,eps )
c
c-----------------------------------------------------------------------
      implicit none
c
			integer n
			integer indx(n)
			real*8 d,eps
      real*8 a(n,n)
c
			integer i,j,k
			integer imax
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
          write (6,*) 'singular matrix in ummdp_ludcmp'
          text = 'matrix detail'
          call ummdp_utility_print2 ( text,a,n,n )
          call ummdp_exit ( 9000 )
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
      end subroutine ummdp_utility_ludcmp
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     LU BACKWARD SUBSTITUTION
c
      subroutine ummdp_utility_lubksb ( a,n,indx,b,eps )
c-----------------------------------------------------------------------
      implicit none
c
			integer n
			integer indx(n)
			real*8 eps
			real*8 a(n,n),b(n)
c
			integer i,j
			integer ii,ll
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
      end subroutine ummdp_utility_lubksb
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     CALCULATE INVERSE MATRIX 2x2 
c
      subroutine ummdp_utility_minv2 ( b,a,deta,eps )
c
c-----------------------------------------------------------------------
      implicit none
c
			real*8,intent(in) :: eps
			real*8,intent(in) :: a(2,2)
c
      real*8,intent(out) :: deta
      real*8,intent(out) :: b(2,2)
c
			real*8 detai
c-----------------------------------------------------------------------
c
      deta = a(1,1)*a(2,2) - a(1,2)*a(2,1)
      if ( abs(deta) <= eps ) then
         write (6,*) 'determinant det[a] error',deta
         write (6,*) 'stop in minv2'
         call ummdp_exit ( 9000 )
      end if
c
      detai = 1.0d0 / deta
      b(1,1) =  a(2,2) * detai
      b(1,2) = -a(1,2) * detai
      b(2,1) = -a(2,1) * detai
      b(2,2) =  a(1,1) * detai
c
      return
      end subroutine ummdp_utility_minv2
c
c
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     CALCULATE INVERSE MATRIX 3x3
c
      subroutine ummdp_utility_minv3 ( b,a,deta,eps )
c-----------------------------------------------------------------------
      implicit none
c
			real*8,intent(in) :: eps
      real*8,intent(in) :: a(3,3)
c
      real*8,intent(out) :: deta 
      real*8,intent(out) :: b(3,3)
c
			real*8 detai
c-----------------------------------------------------------------------
c
      deta = a(1,1) * (a(2,2)*a(3,3) - a(2,3)*a(3,2))
     1       + a(1,2) * (a(2,3)*a(3,1) - a(2,1)*a(3,3))
     2       + a(1,3) * (a(2,1)*a(3,2) - a(2,2)*a(3,1))
      if ( abs(deta) <= eps ) then
         write (6,*) 'determinant det[a] error',deta
         write (6,*) 'stop in minv3'
         call ummdp_exit ( 9000 )
      end if
c
      detai = 1.0d0 / deta
      b(1,1) = (a(2,2)*a(3,3) - a(2,3)*a(3,2)) * detai
      b(1,2) = (a(1,3)*a(3,2) - a(1,2)*a(3,3)) * detai
      b(1,3) = (a(1,2)*a(2,3) - a(1,3)*a(2,2)) * detai
      b(2,1) = (a(2,3)*a(3,1) - a(2,1)*a(3,3)) * detai
      b(2,2) = (a(1,1)*a(3,3) - a(1,3)*a(3,1)) * detai
      b(2,3) = (a(1,3)*a(2,1) - a(1,1)*a(2,3)) * detai
      b(3,1) = (a(2,1)*a(3,2) - a(2,2)*a(3,1)) * detai
      b(3,2) = (a(1,2)*a(3,1) - a(1,1)*a(3,2)) * detai
      b(3,3) = (a(1,1)*a(2,2) - a(1,2)*a(2,1)) * detai
c
      return
      end subroutine ummdp_utility_minv3
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
      subroutine ummdp_utility_eigen_sym3 ( es,ev,a )
c-----------------------------------------------------------------------
      implicit none
c
      real*8 es(3)
			real*8 ev(3,3),a(3,3)
c
			integer i,j,is,ip,iq,ir
			integer msweep
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
        write (6,*) 'stop in ummdp_eigen_sym3'
        call ummdp_exit ( 9000 )
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
      write (6,*) 'stop in ummdp_eigen_sym3'
      call ummdp_exit ( 9000 )
c
      return
      end subroutine ummdp_utility_eigen_sym3
c
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c     CHECKING EXISTENCE OF FILE NAMES 'FLNAME'
c
      logical function ummdp_utility_file_exist ( flname )
c
c-----------------------------------------------------------------------
      implicit none
c
      character*16,intent(in) :: flname
c
			integer nio
c-----------------------------------------------------------------------
c
      nio = 616
      open  ( nio,file=flname,status='old',err=10 )
c
      close ( nio,status='keep' )
      ummdp_utility_file_exist = .true.
      return
c
   10 ummdp_utility_file_exist = .false.
      return
c
      end function ummdp_utility_file_exist
c
c
c

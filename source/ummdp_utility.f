************************************************************************
*
*     UTILITY SUBROUTINES
*
************************************************************************
c
c     ummdp_utility_setunitm
c       set unit 2nd order matrix
c
c     ummdp_utility_print1
c       print vector with text
c
c     ummdp_utility_print2
c       print matrix with text
c
c     ummdp_utility_mv
c       mutiply matrix and vector
c
c     ummdp_utility_mm
c       mutiply matrix and matrix
c
c     ummdp_utility_vvs
c       calculate scalar product of vectors 
c
c     ummdp_utility_minv
c       calculate inverse matrix using lu decomposition
c
c         ummdp_utility_ludcmp
c           lu decomposition
c         ummdp_utility_lubks
c           lu backward substitution
c         ummdp_utility_minv2
c           calculate inverse matrix 2x2
c         ummdp_utility_minv3
c           calculate inverse matrix 3x3
c
c     ummdp_utility_file_exist
c       checking existence of files
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
      a = 0.0d0 
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
      character*100 fmt
c-----------------------------------------------------------------------
c
      write(fmt,'(A,I2,A)') '(/',12+tab,'xA,A)'
      write (6,fmt) '. ',text
c
      write(fmt,'(A,I2,A)') '(',14+tab,'x6E23.15)'
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
      write(fmt,'(A,I2,A)') '(',14+tab,'x6E23.15)'
      do i = 1,n
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
      write(fmt,'(A,I2,A)') '(',14+tab,'x6E23.15)'
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
      v = 0.0d0
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
      a = 0.0d0
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
         call ummdp_exit ( 401 )
      end if
c                                                            ---- B=A^-1
      do j = 1,n
        y = 0.0d0
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
        c = 0.0d0
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
c
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
          call ummdp_exit ( 402 )
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
c
c     LU BACKWARD SUBSTITUTION
c
      subroutine ummdp_utility_lubksb ( a,n,indx,b,eps )
c
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
         call ummdp_exit ( 401 )
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
c
c     CALCULATE INVERSE MATRIX 3x3
c
      subroutine ummdp_utility_minv3 ( b,a,deta,eps )
c
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
         call ummdp_exit ( 401 )
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

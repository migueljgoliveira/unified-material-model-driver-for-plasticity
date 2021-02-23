c***********************************************************************
c
c     UMMDp : Utility Subroutines
c
c***********************************************************************
c     jancae_clear1 ( a,n )
c       clear 1st order tensor (vector) a(n)
c     jancae_clear2 ( a,n,m )
c       clear 2nd order tensor (matrix) a(n,m)
c     jancae_clear3 ( a,n,m,l )
c       clear 3rd order tensor a(n,m,l)
c     jancae_setunitm ( a,n )
c       set unit 2nd oder tensor [I]
c     jancae_print1 ( text,a,n )
c       print 1st order tensor with title (text)
c     jancae_print2 ( text,a,n,m )
c       print 2nd order tensor with title (text)
c     jancae_mv (v,a,u,nv,nu)
c       mutiply matrix to vector v(nv)=a(nv,nu)*u(nu)
c     jancae_mm (a,b,c,na1,na2,nbc)
c       mutiply matrix and matrix a(na1,na2)=b(na1,nbc)*c(nbc,na2)
c     jancae_vvs ( s,u,v,n )
c       calc. scalar (inner) product of v(n) & u(n)
c     jancae_minv ( b,a,n,d )
c       calc. inverse matrix b(n,n)=a(n,n)^-1 and det(a)
c       branch to following routines
c       jancae_ludcmp( a,n,indx,d ) : LU-decomposition
c       jancae_lubksb(a,n,indx,b)   : backward subsitution
c       jancae_minv2 ( b,a,deta )   : for 2*2 matrix
c       jancae_minv3 ( b,a,deta )   : for 3*3 matrix
c     jancae_eigen_sym3 ( es,ev,a )
c       calc. eigen value and normalized eigen vector (3*3sym)
c     jancae_printinfo.
c       print inc.info. and element info.
c     jancae_printinout
c       print input and output of user subroutine
c
c
c-----------------------------------------------------------------------
c     zero clear vector a(n)
c
      subroutine jancae_clear1 ( a,n )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n)
c
      do i = 1,n
        a(i) = 0.0d0
      end do
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     zero clear matrix a(n,m)
c
      subroutine jancae_clear2 ( a,n,m )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n,m)
c
      do i = 1,n
        do j = 1,m
          a(i,j) = 0.0d0
        end do
      end do
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     zero clear matrix a(n,m,l)
c
      subroutine jancae_clear3 ( a,n,m,l )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n,m,l)
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
      end
c
c
c
c-----------------------------------------------------------------------
c     set unit 2nd oder tensor [I]
c
      subroutine jancae_setunitm ( a,n )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n,n)
c
      call jancae_clear2 ( a,n,n )
      do i = 1,n
        a(i,i) = 1.0d0
      end do
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     print vector a(n) with text
c
      subroutine jancae_print1 ( text,a,n )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n)
      character text*32
c
      write (6,*) text
      write (6,9000) (a(i),i=1,n)
 9000 format (6e16.8)
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     print matrix a(n,m) with text
c
      subroutine jancae_print2 ( text,a,n,m )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n,m)
      character text*32
c
      write (6,*) text
      do i = 1,n
        write (6,9000) (a(i,j),j=1,m)
      end do
 9000 format (6e16.8)
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     calculate multiplication of matrix and vector
c     {v}=[a]{u}
c
      subroutine jancae_mv (v,a,u,nv,nu)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension v(nv),a(nv,nu),u(nu)
c
      call jancae_clear1 ( v,nv )
      do i = 1,nv
        do j = 1,nu
          v(i) = v(i) + a(i,j)*u(j)
        end do
      end do
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     calculate multiplication of matrix and matrix
c     [a]=[b][c]
c
      subroutine jancae_mm (a,b,c,na1,na2,nbc)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(na1,na2),b(na1,nbc),c(nbc,na2)
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
      end
c
c
c
c-----------------------------------------------------------------------
c     calculate scaler product of vectors
c     s={u}T{v}
c
      subroutine jancae_vvs ( s,u,v,n )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension v(n),u(n)
c
      s = 0.0
      do i = 1,n
        s = s + u(i)*v(i)
      end do
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     calculate inverse matrix using LU decomposition
c     [b]=[a]-1
c
      subroutine jancae_minv ( b,a,n,d )
c
c     Ref. http://astr-www.kj.yamagata-u.ac.jp/~shibata/kbg/
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n,n),b(n,n)
      dimension indx(n),y(n),c(n,n),aorg(n,n)
      character text*32
      logical   check
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
          if ( anorm .lt. abs(a(i,j)) ) anorm = abs(a(i,j))
        end do
      end do
      do i = 1,n
        do j = 1,n
          a(i,j) = a(i,j) / anorm
        end do
      end do
c
      if ( n .eq. 2 ) then
        call jancae_minv2 ( b,a,d,eps )
        goto 100
      else if ( n .eq. 3 ) then
        call jancae_minv3 ( b,a,d,eps )
        goto 100
      end if
c
      call jancae_ludcmp ( a,n,indx,d,eps )
c                                                 ---- check determinant
      if ( abs(d) .le. eps ) then
         write (6,*) 'determinant det[a] error',d
         write (6,*) 'stop in minv'
         call jancae_exit(9000)
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
      end
c
c
c
c-----------------------------------------------------------------------
c     LU decomposition
c
      subroutine jancae_ludcmp( a,n,indx,d,eps )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n,n),indx(n)
      dimension vtemp(n)
      character text*32
c
      d = 1.0d0
      do i = 1,n
        aamax = 0.0d0
        do j = 1,n
          if ( abs(a(i,j)) .gt. aamax ) aamax = abs(a(i,j))
        end do
        if ( aamax .le. eps ) then
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
          if ( dum .ge. aamax ) then
            imax = i
            aamax = dum
          end if
        end do
        if ( j .ne. imax ) then
          do k = 1,n
            dum = a(imax,k)
            a(imax,k) = a(j,k)
            a(j,k) = dum
          end do
          d = -d
          vtemp(imax) = vtemp(j)
        end if
        indx(j) = imax
c       if ( abs(a(i,j)) .le. eps ) a(i,j) = eps     !2010.07.02 c.out
        if ( j .ne. n ) then
          ajj = a(j,j)                               !2010.07.02 add
          if ( abs(ajj) .le. eps ) ajj = eps         !2010.07.02 add
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
      end
c
c
c
c-----------------------------------------------------------------------
c     LU backward substitution
c
      subroutine jancae_lubksb( a,n,indx,b,eps )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension a(n,n),b(n),indx(n)
c
      ii = 0
      do i = 1,n
        ll = indx(i)
        sum = b(ll)
        b(ll) = b(i)
        if ( ii .ne. 0 ) then
          do j = ii,i-1
            sum = sum - a(i,j)*b(j)
          end do
        else if ( abs(sum) .ge. eps ) then
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
      end
c
c
c
c-----------------------------------------------------------------------
c     calculate inverse matrix[2,2]
c
      subroutine jancae_minv2 ( b,a,deta,eps )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension b(2,2),a(2,2)
c
c
      deta = a(1,1)*a(2,2) - a(1,2)*a(2,1)
      if ( abs(deta) .le. eps ) then
         write (6,*) 'determinant det[a] error',deta
         write (6,*) 'stop in minv2'
         call jancae_exit(9000)
      end if
c
      detai = 1.0d0 / deta
      b(1,1) = a(2,2) * detai
      b(1,2) = -a(1,2) * detai
      b(2,1) = -a(2,1) * detai
      b(2,2) = a(1,1) * detai
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     calculate inverse matrix[3,3]
c
      subroutine jancae_minv3 ( b,a,deta,eps )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
      dimension b(3,3),a(3,3)
c
      deta = a(1,1) * (a(2,2)*a(3,3) - a(2,3)*a(3,2)) +
     &       a(1,2) * (a(2,3)*a(3,1) - a(2,1)*a(3,3)) +
     &       a(1,3) * (a(2,1)*a(3,2) - a(2,2)*a(3,1))
      if ( abs(deta) .le. eps ) then
         write (6,*) 'determinant det[a] error',deta
         write (6,*) 'stop in minv3'
         call jancae_exit(9000)
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
      end
c
c
c
c-----------------------------------------------------------------------
c     print informations for debug (intro)
c
      subroutine jancae_printinfo ( inc,nnrm,nshr )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /jancae1/ne,ip,lay
c
      nttl = nnrm + nshr
c
      write (6,*) '----- JANCAE.UMMDp debug info. -----'
      write (6,*) 'increment=',inc
      write (6,*) 'elem,ip,lay=',ne,ip,lay
      write (6,*) 'nttl,nnrm,nshr=',nttl,nnrm,nshr
      nerr = 0
      if ( nnrm .eq. 3 ) then
        if ( nshr .eq. 3 ) then
          write (6,*) '3d solid element'
        else if ( nshr .eq. 1 ) then
          write (6,*) 'plane strain or axi-sym solid element'
        else
          nerr = nerr + 1
        end if
      else if ( nnrm .eq. 2 ) then
        if ( nshr .eq. 1 ) then
          write (6,*) 'plane stress or thin shell element'
        else if ( nshr .eq. 3 ) then
          write (6,*) 'thick shell element'
        else
          nerr = nerr + 1
        end if
      else
        nerr = nerr + 1
      end if
      if ( nerr .ne. 0 ) then
        write (6,*) 'no supported element type',nnrm,nshr
        call jancae_exit (9000)
      end if
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c     print informations for debug (input/output)
c
      subroutine jancae_printinout ( io,s,de,d,nttl,
     &                               stv,nstv )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension s(nttl),de(nttl),d(nttl,nttl),stv(nstv)
      character text*32
c
      if ( io .eq. 0 ) then
        text = 'initial stresses'
      else
        text = 'updated stresses'
      end if
      call jancae_print1 ( text,s,nttl )
c
      if ( io .eq. 0 ) then
        text = 'initial internal state var.'
      else
        text = 'updated internal state var.'
      end if
      call jancae_print1 ( text,stv,nstv )
c
      if ( io.eq.0 ) then
        text = 'driving strain increment'
        call jancae_print1 ( text,de,nttl )
      else
        text = 'tangent modulus matrix'
        call jancae_print2 ( text,d,nttl,nttl )
      end if
c
      return
      end
c
c
c
c-----------------------------------------------------------------------
c    calculate eigen value and eigen vector by jacobi method
c
c    ( this subroutine is copied from http://www.flagshyp.com/ )
c
c    input  :  a(3,3)  : symmetric matrix to be analyzed
c    output :  es(i)   : i-th eigen value
c              ev(i,3) : normalized eigen vector for i-th eigen value
c
      subroutine jancae_eigen_sym3 ( es,ev,a )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension   a(3,3),es(3),ev(3,3)
      dimension   w(3,3),prc(3,3)
c
      msweep = 100
      eps = 1.0d-8
c
c                                                       ---- preparation
      ax = 0.0d0
      er = 0.0d0
      do i = 1,3
        do j = 1,3
          if ( abs(a(i,j)) .gt. ax ) ax = abs(a(i,j))
          er = er + abs(a(i,j)-a(j,i))
        end do
      end do
      if ( er/ax .gt. eps ) then
        write (6,*) 'a is not symmetric'
        write (6,*) 'stop in jancae_eigen_sym3'
        call jancae_exit (9000)
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
c
c                                                   ---- starts sweeping
c
      do is = 1,msweep
c
        sum = 0.0d0
        do ip = 1,2
          do iq = ip+1,3
            sum = sum + abs( w(ip,iq) )
          end do
        end do
c       write (6,*) 'ite',is,sum,eps
c
c            ---- if the sum of off-diagonal terms is zero evaluates the
c                                                     esches and returns
c
        if ( abs(sum) .lt. eps ) then
          do i = 1,3
            do j = 1,3
              ev(i,j) = prc(j,i)
            end do
            es(i) = es(i)*ax
          end do
          return
        end if
c
c                             ---- performs the sweep in three rotations
c                                         ---- one per off diagonal term
c
        do ip = 1,2
          do iq = ip+1,3
            od = 100.0d0 * abs( w(ip,iq) )
            if ( abs(od) .gt. eps ) then
              hd = es(iq) - es(ip)
c
c                                      ---- evaluates the rotation angle
c
              theta = 0.5d0 * hd / w(ip,iq)
              t = 1.0d0/(abs(theta) + sqrt(1.0d0+theta**2))
              if ( theta .lt. 0.0d0 ) t = -t
c
c                                   ---- re-evaluates the diagonal terms
c
              c = 1.0d0 / sqrt(1.0d0+t**2)
              s = t * c
              tau = s / (1.0d0+c)
              h = t * w(ip,iq)
              es(ip) = es(ip) - h
              es(iq) = es(iq) + h
c
c                     ---- re-evaluates the remaining off-diagonal terms
c
              ir = 6 - ip - iq
              g = w( min(ir,ip),max(ir,ip) )
              h = w( min(ir,iq),max(ir,iq) )
              w( min(ir,ip),max(ir,ip) ) = g - s*(h+g*tau)
              w( min(ir,iq),max(ir,iq) ) = h + s*(g-h*tau)
c
c                                          ---- rotates the eigenvectors
c
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
c
c                              ---- if convergence is not achieved stops
c
      write (6,*) 'did not converge in eigen calculation.'
      write (6,*) 'msweep=',msweep
      write (6,*) 'eps=',eps
      write (6,*) 'sum=',sum
      write (6,*) 'stop in jancae_eigen_sym3'
      call jancae_exit(9000)
c
      end
c
c
c
c-----------------------------------------------------------------------
c     checking existence of file named "flname"
c
       logical function jancae_file_exist ( flname )
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character   flname*16
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
      end
c
c
c

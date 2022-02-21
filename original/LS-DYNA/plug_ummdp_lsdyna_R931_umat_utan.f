      subroutine umat41 (cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,failel,crv,cma,qmat,elsiz,idele)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
c     isotropic elastic material (sample user subroutine)
c
c     Variables
c
c     cm(1)=first material constant, here young's modulus
c     cm(2)=second material constant, here poisson's ratio
c        .
c        .
c        .
c     cm(n)=nth material constant
c
c     eps(1)=local x  strain increment
c     eps(2)=local y  strain increment
c     eps(3)=local z  strain increment
c     eps(4)=local xy strain increment
c     eps(5)=local yz strain increment
c     eps(6)=local zx strain increment
c
c     sig(1)=local x  stress
c     sig(2)=local y  stress
c     sig(3)=local z  stress
c     sig(4)=local xy stress
c     sig(5)=local yz stress
c     sig(6)=local zx stress
c
c     hsv(1)=1st history variable
c     hsv(2)=2nd history variable
c        .
c        .
c        .
c        .
c     hsv(n)=nth history variable
c
c     dt1=current time step size
c     capa=reduction factor for transverse shear
c     etype:
c        eq."solid" for solid elements
c        eq."sph" for smoothed particle hydrodynamics
c        eq."sld2d" for shell forms 14, and 15 (2D solids - axisymmeteric)
c        eq."sldax" for shell forms 13 (2D solids - plane strain)
c        eq."shl_t" for shell forms 25, 26, and 27 (shells with thickness stretch)
c        eq."shell" for all other shell elements plus thick shell forms 1 and 2
c        eq."tshel" for thick shell forms 3 and 5
c        eq."hbeam" for beam element forms 1 and 11
c        eq."tbeam" for beam element form 3 (truss)
c        eq."dbeam" for beam element form 6 (discrete)
c        eq."beam " for all other beam elements
c
c     tt=current problem time.
c
c     temper=current temperature
c
c     failel=flag for failure, set to .true. to fail an integration point,
c            if .true. on input the integration point has failed earlier
c
c     crv=array representation of curves in keyword deck
c
c     cma=additional memory for material data defined by LMCA at 
c       6th field of 2nd crad of *DATA_USER_DEFINED
c
c     elsiz=characteristic element size
c
c     idele=element id
c
c     All transformations into the element local system are
c     performed prior to entering this subroutine.  Transformations
c     back to the global system are performed after exiting this
c     routine.
c
c     All history variables are initialized to zero in the input
c     phase. Initialization of history variables to nonzero values
c     may be done during the first call to this subroutine for each
c     element.
c
c     Energy calculations for the dyna3d energy balance are done
c     outside this subroutine.
c
      include 'nlqparm'
      include 'nhisparm.inc'
c
      integer nnc,lczc,
     &    ns11, ns12, ns13, ns14, ns15, ns16,
     &    nh11, nh12, nh13, nh14, nh15, nh16,
     &    nt11, nt12, nt13, nt14, nt15, nt16,
     &    nb11, nb12, nb13, nb14, nb15, nb16,
     &    nu11, nu12, nu13, nu14, nu15, nu16,
     &    nd11, nd12, nd13, nd14, nd15, nd16,
     &    nd17, nd18, nd19, nd20, nd21, nd22,
     &    ncf26, ncf27, ncf28,
     &    n_em_26, n_em_27, n_em_28,
     &    nscf08, nscf09, nscf10, nscf13, nscf14,
     &    nhcf11, nhcf12, nhcf13, nhcf14, nhcf15,
     &    nh_em_11, nh_em_12, nh_em_13, nh_em_14, nh_em_15,
     &    lccf1h, lccf1s, lc_em_1h
      common/sorter/nnc, lczc,
     &    ns11, ns12, ns13, ns14, ns15, ns16,
     &    nh11, nh12, nh13, nh14, nh15, nh16,
     &    nt11, nt12, nt13, nt14, nt15, nt16,
     &    nb11, nb12, nb13, nb14, nb15, nb16,
     &    nu11, nu12, nu13, nu14, nu15, nu16,
     &    nd11, nd12, nd13, nd14, nd15, nd16,
     &    nd17, nd18, nd19, nd20, nd21, nd22,
     &    ncf26, ncf27, ncf28,
     &    n_em_26, n_em_27, n_em_28,
     &    nscf08, nscf09, nscf10, nscf13, nscf14,
     &    nhcf11, nhcf12, nhcf13, nhcf14, nhcf15,
     &    nh_em_11, nh_em_12, nh_em_13, nh_em_14, nh_em_15,
     &    lccf1h, lccf1s, lc_em_1h
c
      integer
     1 nh01,nh02,nh03,nh04,nh05,nh06,nh07,nh08,nh09,nh10,
     2 nb01,nb02,nb03,nb04,nb05,nb06,nb07,nb08,nb09,nb10,
     3 ns01,ns02,ns03,ns04,ns05,ns06,ns07,ns08,ns09,ns10,
     4 nt01,nt02,nt03,nt04,nt05,nt06,nt07,nt08,nt09,nt10,
     5 nu01,nu02,nu03,nu04,nu05,nu06,nu07,nu08,nu09,nu10,
     6 nd01,nd02,nd03,nd04,nd05,nd06,nd07,nd08,nd09,nd10,
     7 ntbqs
      common/bk05/
     1 nh01,nh02,nh03,nh04,nh05,nh06,nh07,nh08,nh09,nh10,
     2 nb01,nb02,nb03,nb04,nb05,nb06,nb07,nb08,nb09,nb10,
     3 ns01,ns02,ns03,ns04,ns05,ns06,ns07,ns08,ns09,ns10,
     4 nt01,nt02,nt03,nt04,nt05,nt06,nt07,nt08,nt09,nt10,
     5 nu01,nu02,nu03,nu04,nu05,nu06,nu07,nu08,nu09,nu10,
     6 nd01,nd02,nd03,nd04,nd05,nd06,nd07,nd08,nd09,nd10,
     7 ntbqs
c
      integer n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
      common/bk07/n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,n11,n12,n13,n14,n15,
     1 n16,n17,n18,n19,n20,n21,n22,n23,n24,n25,n26,n27,n28,n29,n30,n31,
     2 n32,n33,n34,n35,n36,n37,n38,n39,n40,n41,n42,n43,n44,n44a,n45,
     3 n46,n47,n48,n49,n50,n51,n52,n53,n54,n55,n56,n57,n58,n59,n60,n61,
     4 n62,n63,n64,n65,n66,n67,n68,n69,n70,n71,n72,n73,n74,n75,n76,n77,
     5 n78,n79,n80,n81,n82,locend,iname,lendf
c
      integer*8 dm_hmtnum,dm_hnfegp
      integer*8 dm_bmtnum,dm_bnfegp
      integer*8 dm_smtnum,dm_snfegp
      integer*8 dm_tmtnum,dm_tnfegp
c
      common/dmbk05/
     1 dm_hmtnum,dm_hnfegp,
     2 dm_bmtnum,dm_bnfegp,
     3 dm_smtnum,dm_snfegp,
     4 dm_tmtnum,dm_tnfegp
c
      include 'bk06.inc'
      include 'iounits.inc'
      include 'memaia.inc'
c
      integer i_mem
      common/dynmem/i_mem(1)
      real r_mem(1)
      equivalence (i_mem,r_mem)
      integer*4 dm_x, dm_v, dm_xms, dm_me1
      integer*4 dm_x0,dm_v0,dm_xms0,dm_me10
      common /dynmem1/ dm_x, dm_v, dm_xms, dm_me1,
     &                 dm_x0,dm_v0,dm_xms0,dm_me10
      common/bk13/lc0,lc1h,lc1b,lc1s,lc1t,lc2,lc3,lc4,lc5,lc6,lc7,lc9,
     1 lc10,lc11,lc12,lc13,lc14,lc15,lc16,lc17,lc18,lb0,lb1,lb2,
     2 lc7a,lc7b,lc7c,lc7d,lc7e,lc7f,lc7g,lc7h,lc7i,lc7j,lc7k,lc7l
c
      common/bk00/numnp,numpc,numlp,neq,ndof,nlcur,numcl,numvc,
     + ndtpts,nelmd,nmmat,numelh,numelb,numels,numelt,numdp,
     + grvity,idirgv,nodspc,nspcor,numelh10,numels8
c
c     common /jancae2/ddsdde(6,6)
      common /jancae3/nyf,nih,nkh  !ti160112
c
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*),qmat(3,3)
c     dimension sigjan(6),s2(6),dpe(6),epsde(6)
c     dimension x1(100,6),x2(100,6),pe(6) !ht150908
c     dimension ddsdde(6,6)
      real, allocatable :: sigjan(:),s2(:),dpe(:),epsde(:)
      real, allocatable :: x1(:,:),x2(:,:),pe(:) !ht150908
      real, allocatable :: ddsdde(:,:)
      dimension ddsdde_dyna(6,6)
c     parameter (mxprop=40) !ht150908
      parameter (mxprop=100) !ht150908
      dimension prop(mxprop)
c
      logical failel
      character*5 etype
c
      mxpbs=100
      nprop=mxprop
c     if (ncycle.eq.1) then
c       if (cm(16).ne.1234567) then
c       call usermsg('mat41')
c       end if
c     endif
c
c     if (etype.eq.'solid') then
c       i_element_type=2
c     elseif (etype.eq.'sld2d') then
c       i_element_type=4
c     elseif (etype.eq.'shell') then
c       i_element_type=4
c     else
c       write(*,*)'error etype=',etype
c       stop
c     endif
c
c     i_elem_id=lqfint(idele,i_element_type,ierror)
c     call load_ip(a(lc1s),i_elem_id,ip)
c     call load_iop(ia(n1+nmmat+ip-1),iop)
c
c     write(*,*)'n1=',n1
c     write(*,*)'ns13=',ns13
c     write(*,*)'nmmat=',nmmat
c     write(*,*)'i_elem_id',i_elem_id
c     write(*,*)'idele',idele
c     write(*,*)'ip=',ip
c     write(*,*)'iop=',iop
c     do icheck=601,800
c       write(*,*)'a(',icheck,')=',a(icheck)
c     enddo
c     do icheck=ns13-100,ns13+100
c       write(*,*)'a(',icheck,')=',a(icheck)
c     enddo
c     stop
c
      if (etype.eq.'solid') then
        nnrm = 3
        nshr = 3
      elseif (etype.eq.'sld2d') then
        nnrm = 3
        nshr = 1
      elseif (etype.eq.'shell') then
c       if (iop .eq. 12 ) then
          nnrm = 2
          nshr = 1
c       else
c         nnrm = 2
c         nshr = 3
c       endif
      else
        write(*,*)'error etype=',etype
      endif
c
c     write(*,*)'nnrm=',nnrm
c     write(*,*)'nshr=',nshr
c
      nttl = nnrm + nshr
      ntensor = nttl
c
      allocate( sigjan(nttl),s2(nttl),dpe(nttl),epsde(nttl) )
      allocate( x1(100,nttl),x2(100,nttl),pe(nttl) )!ht150908
      allocate( ddsdde(nttl,nttl) )
c
c     do i=1,6
c     write(*,*)'eps=',eps(i)
c     enddo
c     --- for JANCAE coded by Ida ---
cti160112
c     nyf=cm(9)
c     nih=cm(10)
c     nkh=cm(11)
c     write(*,*)'in umat41 nyf=',nyf
c     write(*,*)'in umat41 nih=',nih
c     write(*,*)'in umat41 nkh=',nkh
cti160112
c
      call jancae_debugmode ( nvbs )
c
c                                               --- output detailed information
c     idebug = cm(40)
c     if ( idebug .ne. 0 ) nvbs = idebug
c
      if ( nvbs.ge.4 ) then
        kinc = 1
        call jancae_printinfo  ( kinc,nnrm,nshr )
        call jancae_printinout ( 0,sig,eps,ddsdde,ntensor,
     &                             hsv,100 )
      endif
c
      do i=1,nprop
        prop(i)=0.0
      enddo
      do i=1,nprop
        prop(i)=cm(i)
      enddo
c       --- ht150908
c                                                  ---- set material properties
c     call jancae_set_prop ( prop,nprop )
c
c     npbs=nint(prop(1))                    ! number of partial back stresses
c     if ( npbs.gt.mxpbs ) then
c       write (6,*) 'npbs > mxpbs error in UMAT41'
c       write (6,*) 'npbs =',npbs
c       write (6,*) 'mxpbs=',mxpbs
c       call jancae_exit ( 9000 )
c     endif
c
      call jancae_prop_dim ( prop,mxprop,
     &                       ndela,ndyld,ndihd,ndkin,
     &                       npbs )
c  
c                                                           -- check nhisvar
      call jancae_check_nisv ( nhisvar,ntensor,npbs ) !ht150908
c                                    ---- copy current internal state variables
      call jancae_isvprof ( isvrsvd,isvsclr )
c
      call jancae_isv2pex ( isvrsvd,isvsclr,
     &                        hsv,100,
     &                        epsp,pe,x1,nttl,mxpbs,npbs ) !ht150908
c
c     call loadshlelf (a(ns06),mx,iop)
c
c     write(*,*)'ns06=',ns06
c     write(*,*)'mx=',mx
c     write(*,*)'iop=',iop
c     do i=1,10000
c       write(*,*)'a',i,'=',a(i)
c     enddo
c
c     stop
c     compute shear modulus, g
c
c     g2 =abs(cm(1))/(1.+cm(2))
c     g  =.5*g2
      b  = cm(2)
      g  = cm(3)
      bg23=b-2.*g/3.
      bg43=b+4.*g/3.
c
c-----ht180207: this g comes from input deck
c               i use this g for out of plane shear modulus
c
c      write(6,888)cm(1)
c 888  format("young's modulus =",e16.8)
c      write(6,889)cm(2)
c 889  format("poisson's ratio =",e16.8)
c      write(6,890)g
c 890  format("shear modulus g =",e16.8)
c
c
c
c
 10   format(/
     1 ' *** Error element type ',a,' can not be',
     2 '           run with the current material model.')
c
c     mjac=-1
      mjac=1
c
      if (etype .eq. 'solid') then
        sigjan(1)=sig(1)
        sigjan(2)=sig(2)
        sigjan(3)=sig(3)
        sigjan(4)=sig(4)
        sigjan(5)=sig(5)
        sigjan(6)=sig(6)
        epsde(1)=eps(1)
        epsde(2)=eps(2)
        epsde(3)=eps(3)
        epsde(4)=eps(4)
        epsde(5)=eps(5)
        epsde(6)=eps(6)
      elseif (etype .eq. 'sld2d') then
        sigjan(1)=sig(1)
        sigjan(2)=sig(2)
        sigjan(3)=sig(3)
        sigjan(4)=sig(4)
        epsde(1)=eps(1)
        epsde(2)=eps(2)
        epsde(3)=eps(3)
        epsde(4)=eps(4)
      elseif (etype .eq. 'shell') then
c       if (iop .eq. 12) then
        sigjan(1)=sig(1)
        sigjan(2)=sig(2)
        sigjan(3)=sig(4)
        epsde(1)=eps(1)
        epsde(2)=eps(2)
        epsde(3)=eps(4)
c       write(*,*)'sigjan1=',sigjan(1)
c       write(*,*)'sigjan2=',sigjan(2)
c       write(*,*)'sigjan3=',sigjan(3)
c       else
c       sigjan(1)=sig(1)
c       sigjan(2)=sig(2)
c       sigjan(3)=sig(4)
c       sigjan(4)=sig(5)
c       sigjan(5)=sig(6)
c       epsde(1)=eps(1)
c       epsde(2)=eps(2)
c       epsde(3)=eps(4)
c       epsde(4)=eps(5)
c       epsde(5)=eps(6)
c       endif
      endif
c
c
c
c
c-----ht180207--------------------------------DYNA's shell is thin shell but nshr=3
c
c     if ( (nnrm.eq.2).and.(nshr.eq.3) ) then
c                    * thick to thin shell
c       nnrm=2
c       nshr=1
c       nttl=nnrm+nshr
c
c       deallocate( sigjan,s2,dpe,epsde )
c       deallocate( x1,x2,pe )!ht150908
c       deallocate( ddsdde )
c       allocate( sigjan(nttl),s2(nttl),dpe(nttl),epsde(nttl) )
c       allocate( x1(100,nttl),x2(100,nttl),pe(nttl) )!ht150908
c       allocate( ddsdde(nttl,nttl) )
c     endif
c
      call jancae_plasticity  ( sigjan,s2,epsde,
     &                          epsp,dp,dpe,de33,
     &                          x1,x2,mxpbs,
     &                          ddsdde,
     &                          nnrm,nshr,nttl,
     &                          nvbs,mjac,
     &                          prop,nprop )  !ht150908
c
c     do i=1,nttl
c       do j=1,nttl
c         write(*,*)'ddsdde',i,j,'=',ddsdde(i,j)
c       enddo
c     enddo
c
c     if (etype .eq. 'shell') then
c     if (iop .ne. 12 ) then
c
c       do i=4,5
c         dpe(i)=0.0
c         do n=1,mxpbs
c           x2(n,i)=0.0
c         enddo
c       enddo
c       ddsdde(4,4)=g
c       ddsdde(5,5)=g
c     endif
c     endif
c
      if ( etype .eq. 'shell') then
        sig(1)=s2(1)
        sig(2)=s2(2)
        sig(3)=0.0
        sig(4)=s2(3)
c       if ( iop .ne. 12 ) then
        sig(5)=sig(5)+g*eps(5)
        sig(6)=sig(6)+g*eps(6)
c       endif
        eps(3) = de33
      elseif ( etype .eq. 'solid' .or. etype .eq. 'sld2d') then
        sig(1)=s2(1)
        sig(2)=s2(2)
        sig(3)=s2(3)
        sig(4)=s2(4)
        if ( etype .eq. 'solid') then
        sig(5)=s2(5)
        sig(6)=s2(6)
        endif
      endif
c     do itensor = 1,nttl
c       sig(itensor)=s2(itensor)
c     enddo
c
      epsp=epsp+dp
c
      call jancae_isvprof ( isvrsvd,isvsclr )
c
c                                 ---- update plast.strain comp.
      hsv(isvrsvd+1)=epsp
c                                     -- plastic strain comp.
      do i=1,nttl
        hsv(isvrsvd+isvsclr+i)=hsv(isvrsvd+isvsclr+i)+dpe(i)
      enddo
c
c                                -- partial back stress comp.
      if ( npbs.ne.0 ) then   !ht150908
        do nb=1,npbs          !ht150908
          do i=1,nttl
            it=isvrsvd+isvsclr+nttl*nb+i
            hsv(it)=x2(nb,i)
          enddo
        enddo
      endif
c
c
c     changing LS-DYNA MATRIX(6,6)
      do i=1,6
        do j=1,6
          ddsdde_dyna(i,j)=0.0d0
        enddo
      enddo
      if ( etype .eq. 'solid') then
      do i=1,6
        do j=1,6
          ddsdde_dyna(i,j)=ddsdde(i,j)
        enddo
      enddo
      elseif ( etype .eq. 'sld2d') then
      do i=1,4
        do j=1,4
          ddsdde_dyna(i,j)=ddsdde(i,j)
        enddo
      enddo
      ddsdde_dyna(5,5)=g
      ddsdde_dyna(6,6)=g
      elseif ( etype .eq. 'shell') then
      do i=1,3
        do j=1,3
          i2=i
          j2=j
          IF (i2.eq.3) i2=4
          if (j2.eq.3) j2=4
          ddsdde_dyna(i2,j2) = ddsdde(i,j)
        enddo
      enddo
      ddsdde_dyna(3,3)=9.0*b*g/(3*b+g) ! dummy : young's modulus 
      ddsdde_dyna(5,5)=g
      ddsdde_dyna(6,6)=g
c     if ( iop .eq. 12) then
c     do i=1,2
c       do j=1,2
c         ddsdde_dyna(i,j) = ddsdde(i,j)
c       enddo
c     enddo
c     else
c     ddsdde_dyna(4,4)=ddsdde(3,3)
c     ddsdde_dyna(1,4)=ddsdde(1,3)
c     ddsdde_dyna(4,1)=ddsdde(3,1)
c     ddsdde_dyna(2,4)=ddsdde(2,3)
c     ddsdde_dyna(4,2)=ddsdde(3,2)
c     ddsdde_dyna(1,3)=bg23
c     ddsdde_dyna(3,1)=bg23
c     ddsdde_dyna(3,2)=bg23
c     ddsdde_dyna(3,3)=bg43
c     ddsdde_dyna(5,5)=g
c     ddsdde_dyna(6,6)=g
c     endif
      endif
c
c     write(*,*)'no_hsvs=',no_hsvs
c     write(*,*)'ii_start=',ii_start
      ii_start=30
      do i = 1,6
        do j = 1,6
          ii_start = ii_start + 1
          hsv(ii_start)=ddsdde_dyna(i,j)
        enddo
      enddo
c
 999  format(6e16.8)
      return
      end
c
      subroutine utan41(cm,eps,sig,epsp,hsv,dt1,capa,etype,tt,
     1 temper,es,crv,failel,cma,qmat)
c
c******************************************************************
c|  Livermore Software Technology Corporation  (LSTC)             |
c|  ------------------------------------------------------------  |
c|  Copyright 1987-2008 Livermore Software Tech. Corp             |
c|  All rights reserved                                           |
c******************************************************************
c
      include 'nlqparm'
      dimension cm(*),eps(*),sig(*),hsv(*),crv(lq1,2,*),cma(*)
      dimension es(6,*),qmat(3,3)
      logical failel
      character*5 etype
      character*32 text
c
c
      factor=1.
      if (failel) factor=1.e-8
c
      ii_start=30
      do i = 1,6
        do j = 1,6
          ii_start = ii_start + 1
          es(i,j)=hsv(ii_start)
        enddo
      enddo
c
      call jancae_debugmode ( nvbs )
c
c     idebug = cm(40)
c     if ( idebug .ne. 0 ) nvbs = idebug
c
      if ( nvbs.ge.5 ) then
        text='Restore ES matrix'
        call jancae_print2 ( text,es,6,6 )
      endif
c
c     do i=1,6
c       write(*,*)'eps',i,'=',eps(i)
c     enddo
c
      return
      end
c-------------------------------------------------------------
c     set internal state variables profile
c
      subroutine jancae_isvprof ( isvrsvd,isvsclr )
c-------------------------------------------------------------
c
c               no reserved variables
      isvrsvd=0
c               
      isvsclr=1
c
      return
      end
c
c
c-------------------------------------------------------------
c     exit program by error
c
      subroutine jancae_exit (nexit)
c-------------------------------------------------------------
      common /jancae1/ne,ip,lay
c                                  nexit : exit code
      write (6,*) 'error code :',nexit
      write (6,*) 'element no.           :',ne
      write (6,*) 'integration point no. :',ip
      write (6,*) 'layer no.             :',lay
c
      call adios(nexit)
c
      return
      end
c
c     subroutine load_ip(ixs,i,ip)
c
c     dimension ixs(5,*)
c
c     ip=ixs(1,i)
c
c     return
c     end
c
c     subroutine load_iop(ia_iop,iop)
c
c     iop=ia_iop
c
c     return
c     end

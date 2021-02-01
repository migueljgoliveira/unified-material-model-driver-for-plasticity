c
c      JANCAE.UMMDp/MSC.Marc
c
c
c
c
c
c
c
c *************   user subroutine for defining material behavior  **************
c
c
CAUTION : Due to calculation of the Deformation gradients, Stretch Tensors and
c         Rotation tensors at previous and current states, the analysis can be
c         computationally expensive. Please use the user subroutine ->  hypela
c         if these kinematic quantities are not needed in the constitutive model
c
c
c IMPORTANT NOTES :
c
c (1) F,R,U are only available for continuum and membrane elements (not for
c     shells and beams).
c
c (2) For total Lagrangian formulation use the -> 'Elasticity,1' card(=
c     total Lagrange with large disp) in the parameter section of input deck.
c     For updated Lagrangian formulation use the -> 'Plasticity,3' card(=
c     update+finite+large disp+constant d) in the parameter section of
c     input deck.
c
c
c     d            stress strain law to be formed
c     g            change in stress due to temperature effects
c     e            total elastic strain
c     de           increment of strain
c     s            stress - should be updated by user
c     t            state variables (comes in at t=n, must be updated
c                                   to have state variables at t=n+1)
c     dt           increment of state variables
c     ngens        size of stress - strain law
c     n            element number
c     nn           integration point number
c     kcus(1)      layer number
c     kcus(2)      internal layer number
c     matus(1)     user material identification number
c     matus(2)     internal material identification number
c     ndi          number of direct components
c     nshear       number of shear components
c     disp         incremental displacements
c     dispt        displacements at t=n   (at assembly,        lovl=4) and
c                  displacements at t=n+1 (at stress recovery, lovl=6)
c     coord        coordinates
c     ncrd         number of coordinates
c     ndeg         number of degrees of freedom
c     itel         dimension of F and R, either 2 or 3
c     nnode        number of nodes per element
c     jtype        element type
c     lclass       element class
c     ifr          set to 1 if R has been calculated
c     ifu          set to 1 if strech has been calculated
c
c --- for md_hypela2
c     nstats       number of state variables
c     matnamec     material name
c     isunit       system units
c                  =0 not entered
c                  =1 - SI-m unit  (N,m,S,C)
c                  =2 - SI-mm unit (N,mm,S,C)
c                  =3 - US(British) unit (lbf,inch,S,F)
c     nid          number of auxiallary integer data
c     nrd          number of auxiallary real data
c     ncd          number of auxiallary character data
c     idata        integer data
c     rdata        real data
c     cdata        character data
c
c     at t=n   :
c
c     ffn          deformation gradient
c     frotn        rotation tensor
c     strechn      square of principal stretch ratios, lambda(i)
c     eigvn(i,j)   i principal direction components for j eigenvalues
c
c     at t=n+1 :
c
c     ffn1         deformation gradient
c     frotn1       rotation tensor
c     strechn1     square of principal stretch ratios, lambda(i)
c     eigvn1(i,j)  i principal direction components for j eigenvalues
c
c     The following operation obtains U (stretch tensor) at t=n+1 :
c
c     call scla(un1,0.d0,itel,itel,1)
c     do 3 k=1,3
c      do 2 i=1,3
c       do 1 j=1,3
c        un1(i,j)=un1(i,j)+dsqrt(strechn1(k))*eigvn1(i,k)*eigvn1(j,k)
c1      continue
c2     continue
c3    continue
c
c
c------------------------------------------------------start of program
c
      subroutine md_hypela2 (
     &                   d,g,e,de,s,t,dt,ngens,n,nn,kcus,matus,ndi,
     &                   nshear,
     &                   disp,dispt,coord,ffn,frotn,strechn,eigvn,ffn1,
     &                   frotn1,strechn1,eigvn1,ncrd,itel,ndeg,ndm,
     &                   nnode,jtype,lclass,ifr,ifu,
     &                   nstats,matnamec,isunit,
     &                   nid,nrd,ncd,idata,rdata,cdata )
c
c
      include 'implicit'
c                                  for real*8 (a-h,o-z)
      include 'concom'
c                                  for inc,ncycle,lovl
c     include 'hards'
c                                  for nstats
c      ( nstats is added to arguments of md_hypela2. )
c
      common /jancae1/ne,ip,lay
      common /jancae3/prop
c
      dimension e(*),de(*),t(*),dt(*),g(*),d(ngens,*),s(*)
      dimension n(2),coord(ncrd,*),disp(ndeg,*),matus(2),
     *          dispt(ndeg,*),ffn(itel,*),frotn(itel,*),
     *          strechn(itel),eigvn(itel,*),ffn1(itel,*),
     *          frotn1(itel,*),strechn1(itel),eigvn1(itel,*),
     *          kcus(2),lclass(2)
c
      character*16 cdata
      character*24 matnamec
      dimension idata(*),rdata(*),cdata(*)     
c
      parameter (mxpbs=10)
      dimension s2(ngens),dpe(ngens),x1(mxpbs,ngens),x2(mxpbs,ngens),
     &          pe(ngens),t2(nstats)
c
      parameter (mxprop=100)
      dimension prop(mxprop)
c
      character text*32
c
c     write (6,*) '----------------',n(1),nn,kcus
c     write (6,*) 'frotn'
c     do i=1,itel
c       write (6,*) (frotn(i,j),j=1,itel)
c     enddo
c     write (6,*) 'frotn1'
c     do i=1,itel
c       write (6,*) (frotn1(i,j),j=1,itel)
c     enddo
c     write (6,*) 'ffn'
c     do i=1,itel
c       write (6,*) (ffn(i,j),j=1,itel)
c     enddo
c     write (6,*) 'ffn1'
c     do i=1,itel
c       write (6,*) (ffn1(i,j),j=1,itel)
c     enddo
c
c     text='frotn'
c     call jancae_print2 ( text,frotn,itel,itel )
c     text='frotn1'
c     call jancae_print2 ( text,frotn1,itel,itel )
c
c     text='ffn'
c     call jancae_print2 ( text,ffn,itel,itel )
c     text='ffn1'
c     call jancae_print2 ( text,ffn1,itel,itel )
c
c
      nprop=mxprop
c                                                      --- variable name change
c                        ne  : element no.
c                        ip  : integration point no.
c                        lay : layer no. of shell
      ne =n(1)
      ip =nn
      lay=kcus(1)
c                                                --- set debug and verbose mode
      call jancae_debugmode ( nvbs )
c                                               --- output detailed information
      if ( nvbs.ge.4 ) then
        if ( lovl.eq.4 ) write (6,*) 'assembly phase ---------'
        if ( lovl.eq.6 ) write (6,*) 'recovery phase ---------'
        call jancae_printinfo  ( inc,ndi,nshear )
        call jancae_printinout ( 0,s,de,d,ngens,
     &                           t,nstats )
      endif
c                                                  ---- set material properties
      if ( .true. ) then
        do i=1,nprop
          prop(i)=0.0
        enddo
        do i=1,nrd
          prop(i)=rdata(i)
        enddo
      else
        call jancae_set_prop ( prop,nprop )
      endif
c
      call jancae_prop_dim ( prop,mxprop,
     &                       ndela,ndyld,ndihd,ndkin,
     &                       npbs )
c                                                               -- check nstats
      call jancae_check_nisv ( nstats,ngens,npbs )
c                                    ---- copy current internal state variables
      call jancae_isvprof ( isvrsvd,isvsclr )
      call jancae_isv2pex ( isvrsvd,isvsclr,
     &                      t,nstats,
     &                      p,pe,x1,ngens,mxpbs,npbs )
c                                                             *** update stress 
c                                                       and set tangent modulus
      mjac= 1
c
      if ( .true. ) then
      call jancae_plasticity  ( s,s2,de,
     &                          p,dp,dpe,de33,
     &                          x1,x2,mxpbs,
     &                          d,
     &                          ndi,nshear,ngens,
     &                          nvbs,mjac,
     &                          prop,nprop )
      else
c                                                      * use hyplas example
      call jancae_hyplasvm    ( s,s2,de,
     &                          p,dp,dpe,de33,
     &                          x1,x2,mxpbs,
     &                          d,
     &                          ndi,nshear,ngens,
     &                          nvbs,mjac,
     &                          prop,nprop )
      endif
c                                                            ---- update stress
      do i=1,ngens
        s(i)=s2(i)
      enddo
c                                                          ---- state variables
c                                                     * inc. of eq.plast,strain
      dt(isvrsvd+1)=dp
c                                                  * inc. of plast.strain comp.
      do i=1,ngens
        it=isvrsvd+isvsclr+i
        dt(it)=dpe(i)
      enddo
c                                                   * inc. of back stress comp.
      if ( npbs.ne.0 ) then
        do nb=1,npbs
          do i=1,ngens
            it=isvrsvd+isvsclr+ngens*nb+i
            dt(it)=x2(nb,i)-x1(nb,i)
          enddo
        enddo
      endif
c                                  ----  if debug mode, output return arguments
      if ( nvbs.ge.4 ) then
        do i=1,nstats
          t2(i)=t(i)+dt(i)
        enddo
        call jancae_printinout ( 1,s,de,d,ngens,
     &                           t2,nstats )
      endif
c
      return
      end
c
c
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c
c     define a variable for contour plotting (user subroutine).
c
c     v            variable to be put onto the post file
c     s (idss)     stress array
c     sp           stresses in preferred direction
c     etot         total strain (generalized)
c     eplas        total plastic strain
c     ecreep       total creep strain
c     t            array of state variable (temperature first)
c     m(1)         user element number
c     m(2)         internal element number
c     m(3)         material id
c     m(4)         internal material id
c     nn           integration point number
c     kcus(1)      layer number
c     kcus(2)      internal layer number
c     ndi          number of direct stress components
c     nshear       number of shear stress components
c     jpltcd       the absolute value of the user's entered post code
c
c
c------------------------------------------------------start of program
c
      subroutine plotv(v,s,sp,etot,eplas,ecreep,t,m,nn,kcus,ndi,
     *                 nshear,jpltcd)
c
      include 'implicit'
      include 'hards'
c                                  for nstats
      common /jancae1/ne,ip,lay
      common /jancae3/prop
c
      dimension s(*),etot(*),eplas(*),ecreep(*),sp(*),m(4),kcus(2),
     *          t(*)
c
      dimension dseds(ndi+nshear),d2seds2(ndi+nshear,ndi+nshear)

      parameter (mxpbs=10)
      dimension x(mxpbs,ndi+nshear),xsum(ndi+nshear),
     &          eta(ndi+nshear),pe(ndi+nshear)
c
      parameter (mxprop=100)
      dimension prop(mxprop)
      real*8,allocatable,dimension(:) :: prela,pryld,prihd,prkin
c
      ne=m(1) 
      ip=nn
      lay=kcus(1)
      njp=jpltcd
      ngens=ndi+nshear
c
c      variables list :
c        plot codes and location of state variables arrays
c        cal : be calculated in this subroutine
c        nd  : ndi
c        ng  : ngens
c
c      jpltcd svcode          name of variable
c
c           1 2               equivalent plastic strain
c           2 cal             equivalent stress (X-direction equivalent)
c           3 cal             flow stress (function of plastic strain)
c
c          11 2+1             plastic strain component pe11 
c          12 2+2                 ..                   pe22
c          13 2+3                 ..                   pe33
c          14 2+nd+1              ..                   pg12=pe12+pe21
c          15 2+nd+2              ..                   pg23=pe23+pe32
c          16 2+nd+3              ..                   pg31=pe31+pe13
c
c          21 cal             back stress component x11
c          22 cal                 ..                x22
c          23 cal                 ..                x33
c          24 cal                 ..                x12
c          25 cal                 ..                x23
c          26 cal                 ..                x31
c
c  (2+i)*10+1 2+ng*(i+1)+1    partial back stress component xi-11
c  (2+i)*10+2 2+ng*(i+1)+2        ..                        xi-22
c  (2+i)*10+3 2+ng*(i+1)+3        ..                        xi-33
c  (2+i)*10+4 2+ng*(i+1)+nd+1     ..                        xi-12
c  (2+i)*10+5 2+ng*(i+1)+nd+2     ..                        xi-23
c  (2+i)*10+6 2+ng*(i+1)+nd+3     ..                        xi-31
c
c
c
c                                                  ---- set material properties
c     call jancae_set_prop ( prop,nprop )
c 
      call jancae_prop_dim ( prop,mxprop,
     &                       ndela,ndyld,ndihd,ndkin,
     &                       npbs )
c
      allocate( prela(ndela) )
      allocate( pryld(ndyld) )
      allocate( prihd(ndihd) )
      allocate( prkin(ndkin) )
      k=0
      do i=1,ndela
        k=k+1
        prela(i)=prop(k)
      enddo
      do i=1,ndyld
        k=k+1
        pryld(i)=prop(k)
      enddo
      do i=1,ndihd
        k=k+1
        prihd(i)=prop(k)
      enddo
      do i=1,ndkin
        k=k+1
        prkin(i)=prop(k)
      enddo
c
      call jancae_isvprof ( isvrsvd,isvsclr )
      call jancae_isv2pex ( isvrsvd,isvsclr,
     &                      t,nstats,
     &                      p,pe,x,ngens,mxpbs,npbs )
      do i=1,ngens
        xsum(i)=0.0
      enddo
      if ( npbs.ne.0 ) then
        do i=1,ngens
          do nb=1,npbs
            xsum(i)=xsum(i)+x(nb,i)
          enddo
        enddo
      endif
c
c                                                         --- eq. plastic strain
      if ( njp.eq.1 ) then
        v=p
c                                                                 --- eq. stress
      else if ( njp.eq.2 ) then
        do i=1,ngens
          eta(i)=sp(i)-xsum(i)
        enddo
        call jancae_yfunc ( se,dseds,d2seds2,0,
     &                      eta,ngens,ndi,nshear,
     &                      pryld,ndyld )


        v=se
c                                                                --- flow stress
      else if ( njp.eq.3 ) then
        call jancae_hardencurve ( sy,dsydp,d2sydp2,
     &                            0,p,prihd,ndihd )
        v=sy
c                                                       --- plastic strain comp.
      else if ( (njp.ge.11).and.(njp.le.16) ) then
        v=0.0
        n0=10
        nc=njp-n0
        if ( nc.le.3 ) then
          if ( nc.le.ndi ) then
            v=pe(nc)
          else
            v=-pe(1)-pe(2)
          endif
        else
          if ( nc-3.le.nshear ) then
            v=pe(ndi+nc-3)
          endif
        endif
      endif
c
c                                                      * for kinematic hardening
      if ( npbs.ne.0 ) then
c                                                          --- back stress comp.
        if ( (njp.ge.21).and.(njp.le.26) ) then
          v=0.0
          n0=20
          nc=njp-n0
          if ( nc.le.3 ) then
            if ( nc.le.ndi ) then
              v=xsum(nc)
            else
              v=-xsum(1)-xsum(2)
            endif
          else
            if ( nc-3.le.nshear ) then
              v=xsum(ndi+nc-3)
            endif
          endif
        endif
c                                                 ---- partial back stress comp.
        do nb=1,npbs
          k=(2+nb)*10
          if ( (njp.ge.k+1).and.(njp.le.k+6) ) then
            v=0.0
            nc=njp-k
            if ( nc.le.3 ) then
              if ( nc.le.ndi ) then
                v=x(nb,nc)
              else
                v=-x(nb,1)-x(nb,2)
              endif
            else
              if ( nc-3.le.nshear ) then
                v=x(nb,ndi+nc-3)
              endif
            endif
          endif
        enddo
c
      endif 
c
      return
      end
c
c
c
c-------------------------------------------------------------
c     set internal state variables profile
c
      subroutine jancae_isvprof ( isvrsvd,isvsclr )
c-------------------------------------------------------------
      include 'implicit'
c                                  for real*8 (a-h,o-z)
c               t(1) is reserved for temperature
      isvrsvd=1
c               t(2) is for eq.plast.strain
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
      include 'implicit'
      common /jancae1/ne,ip,lay
c                                  for real*8 (a-h,o-z)
c                                  nexit : exit code
c
      write (6,*) 'element no.           :',ne
      write (6,*) 'integration point no. :',ip
      write (6,*) 'layer no.             :',lay
      call quit(nexit)
c
      return
      end
c
c
c

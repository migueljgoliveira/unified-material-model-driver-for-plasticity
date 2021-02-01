      SUBROUTINE UCMAT2 [DLLEXPORT] (NG,NEL,IPT,IT2D,IDEATH,ITE,IUPDTL,
     1                   STRESS,EPS,STRAIN,DEPS,DEPST,THSTR1,THSTR2,
     2                   KTR,INTER,SCP,ARRAY,LGTH1,IARRAY,LGTH2,LGTH3,
     3                   LGTH4,D,ALFA,CTD,ALFAA,CTDD,CTI,TMP1,TMP2,
     4                   TIME,ETIMV,ETIMV2,DT,PHIST,PRST,RN,ANGLE,
     5                   DPSP,TGRAD,INTEG,ISUBM,INDNL,LCTI,NSYMST,IERR,
     6                   DP,NELP,DPJE1D,DPJE2D,AKAPPA,PBAR,
     7                   NNODE,NDM2,NODNUM,YZ,DCA,OMEGAJ,WVTEMP,
     8                   IOUT,KEY)                                  !ty150908
c
c################################################################
c 
c About this SUBROUTINE:
c 
c  # NAME
c     ADINA SUBROUTINE for JANCAE-MMSM-MPWG
c
c  
c  # VERSION NUMBER
c     Ver. 2.8
c
c
c  # UPDATE DATE
c     15/09/08
c
c
c  # EDIT BY
c     T.Yamanashi, Newtonworks
c
c
c  # OBJECT
c     jancae_plasticity
c
c    
c  # ENVELOPMENT
c     ADINA Version: 9.0.7
c     OS: Windows 7 64bits
c     Compiler: Intel Visual Fortran Compiler XE 12.1
c     * recommended compiler is Intel Fortran 11.0
c
c
c  # SPECIFICATION
c     1. SET-AXES-MAT command must set.
c     2. IT2D parameter must set.
c        axissymetric : 0
c        plane strain : 1
c        plane stress : 2
c
c################################################################
c
c Parameters:
c
c  # DEPS(6): increment of strain comp.
c       - input from ADINA
c         components defined with global coordinate
c
c  # STRESS(6): stress comp.
c       - output to ADINA
c         components defined with global coordinate
c
c  # ARRAY(1): equivalent plastic strain     (mod ht180310)
c       - store for next step
c
c  # ARRAY(2)...(1+ngens): plastic strain comp.  (mod ht180310)
c       - store for next step
c         components defined with local coordinate 
c
c  # ARRAY(1+ngens)...(1+2*ngens): stress comp. (mod ht180310)
c       - store for next step , input stress
c         components defined with local coordinate 
c
c  # D(6,6): consistent matrix
c       - output to ADINA
c         components defined with global coordinate
c
c  # dlocal(6,6): consistent matrix
c         components defined with local coordinate        
c
c  # DCA(3,2): initial element material axes(a,b)
c       - input from ADINA
c
c  # dca3(3,3): initial element material axes(a,b,c) 
c       - a,b is same to DCA(3,2).
c         c is vector product of a and b.
c
c  # qt(6,6): transfom matrix
c             strain to local coordinate from global
c
c  # qtt(6,6): transpose of qt
c
c  # LGTH3: length of working real array 
c           to write to porthole file
c
c  # LGTH1: length of working real array 
c           for storing history dependent variables
c           (in this subroutine, equal to LGTH3+36)
c
c  # NEL: element number (== ne)
c              
c  # IPT: integration point number (== ip)
c
c  # lay: shell layer number (== 1)
c
c  # DPSP(6):
c       - for this subroutine, this parameter must set to 0.0.
c
c  # ndi: no. of normal components
c         3 for UCMAT3
c
c  # nshear: no. of shear  components
c         3 for UCMAT3
c
c  # ngens: total number of components =nnrm+nshr
c         6 for UCMAT3
c            
c
c################################################################
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      common /jancae1/ne,ip,lay
c
      parameter (mxpbs=10)   !ht150908
      parameter (mxprop=100)  !ht150908
c
      DIMENSION STRESS(4),DEPS(4),D(4,4),DCA(2),DPSP(4)
      DIMENSION ARRAY(*),IARRAY(*)
      DIMENSION CTI(99)
c
      DIMENSION s1(4),s2(4),de(4),dpe(4),dca3(3,3)
      DIMENSION x1(mxpbs,4),x2(mxpbs,4)    !ht150908
      DIMENSION qt(4,4),qtt(4,4),dlocal(4,4)
      DIMENSION qt2(3,3),qtt2(3,3),dlocal2(3,3)
      DIMENSION prop(mxprop)     !ht150908
      DIMENSION tarray(16)       !ty180308
c
      ne =NEL
      ip =IPT
      lay=1
c
c
      nporp=mxprop
c
c
c     IT2D:0   axi-symetric
c     IT2D:1   plane strain
c     IT2D:2   plane stress
c
      if (IT2D.le.1) then
        ndi=3
      else
        ndi=2
      endif
      nshear=1
      ngens=ndi+nshear
c
      call jancae_isvprof ( isvrsvd,isvsclr )
c
c ---------------------------- K E Y = 0 -------------------------------
c
c  This KEY is only for setting of working arrays. 
c
      maxisv=isvrsvd+isvsclr+ngens*2
      if ( LGTH1.lt.maxisv ) then
        write (6,*) 'LGTH1=',LGTH1
        write (6,*) 'isvrsvd=',isvrsvd
        write (6,*) 'isvsclr=',isvsclr
        write (6,*) 'ngens  =',ngens
        write (6,*) 'maxisv =',maxisv
        write (6,*) 'LGT1 is too small'
        call jancae_exit (1000)
      endif
c
c      IF (KEY.EQ.0) THEN
c          LGTH3=13
c          LGTH1=(LGTH3+4)+ngens*ngens
c          LGTH2=0
c          LGTH4=0
c          RETURN
c      ENDIF
c
c ----------------------------------------------------------------------
c
      GO TO (1,2,3,4), KEY
c
c ---------------------------- K E Y = 1 -------------------------------
c
    1 CONTINUE
c
      do i=1,LGTH1
         ARRAY(I)=0.D0
      enddo
      do i=1,LGTH2
         IARRAY(I)=0
      enddo
c
      if (IT2D.le.1) then
        qt(1,1)=DCA(1)*DCA(1)
        qt(1,2)=DCA(2)*DCA(2)
        qt(1,3)=0.0d0
        qt(1,4)=DCA(1)*DCA(2)
        qt(2,1)=DCA(2)*DCA(2)
        qt(2,2)=DCA(1)*DCA(1)
        qt(2,3)=0.0d0
        qt(2,4)=-DCA(2)*DCA(1)
        qt(3,1)=0.0d0
        qt(3,2)=0.0d0
        qt(3,3)=1.0d0
        qt(3,4)=0.0d0
        qt(4,1)=-2.0d0*DCA(1)*DCA(2)
        qt(4,2)=2.0d0*DCA(1)*DCA(2)
        qt(4,3)=0.0d0
        qt(4,4)=DCA(1)*DCA(1)-DCA(2)*DCA(2)
        do i=1,4
        do j=1,4
           qtt(j,i)=qt(i,j)
        enddo
        enddo
      else
        qt2(1,1)=DCA(1)*DCA(1)
        qt2(1,2)=DCA(2)*DCA(2)
        qt2(1,3)=DCA(1)*DCA(2)
        qt2(2,1)=DCA(2)*DCA(2)
        qt2(2,2)=DCA(1)*DCA(1)
        qt2(2,3)=-DCA(2)*DCA(1)
        qt2(3,1)=-2.0d0*DCA(1)*DCA(2)
        qt2(3,2)=2.0d0*DCA(1)*DCA(2)
        qt2(3,3)=DCA(1)*DCA(1)-DCA(2)*DCA(2)
        do i=1,3
        do j=1,3
           qtt2(j,i)=qt2(i,j)
        enddo
        enddo
      endif
c

c
      RETURN
c
c ---------------------------- K E Y = 2 -------------------------------
c
    2 CONTINUE   
c
      if (KTR.eq.0) goto 100
c
      ioffset=isvrsvd+isvsclr+ngens
      do i=1,ngens
         s1(i)=ARRAY(i+ioffset)
      enddo
c
c---- ht171230
c
      do i=1,LCTI
        prop(i)=CTI(i)
      enddo
      call jancae_prop_dim ( prop,nprop,
     &                       ndela,ndyld,ndihd,ndkin,
     &                       npbs )
c
c------ ty171230
c
      mjac=1
c
      if (IT2D.le.1) then
c
        tarray(1)=DEPS(3)
        tarray(2)=DEPS(4)
        DEPS(3)=tarray(2)
        DEPS(4)=tarray(1)
c
        call jancae_mv (de,qt,DEPS,ngens,ngens)
c
        call jancae_plasticity ( s1,s2,de,
     &                           ARRAY(1),dpp,dpe,de33,
     &                           x1,x2,mxpbs,
     &                           dlocal,
     &                           ndi,nshear,ngens,
     &                           nvbs,mjac,
     &                           prop,nprop ) ! ht180310
c
      else
c
        call jancae_mv (de,qt2,DEPS,ngens,ngens)
c
        call jancae_plasticity ( s1,s2,de,
     &                           ARRAY(1),dpp,dpe,de33,
     &                           x1,x2,mxpbs,
     &                           dlocal2,
     &                           ndi,nshear,ngens,
     &                           nvbs,mjac,
     &                           prop,nprop )  ! ht180310
c
      endif
c
c
c ### store updated eq. plastic strain in ARRAY ###
c
      ARRAY(1)=ARRAY(1)+dpp    ! ht180310
c
c ### stire updated plastic strain in ARRAY ###
c
      do i=1,ngens
         ARRAY(1+i)=ARRAY(1+i)+dpe(i)   ! ht180310
      enddo
c
c ### store updated stress(local) in ARRAY ###
c
      do i=1,ngens
         ARRAY(i+1+ngens)=s2(i)  ! ht180310
      enddo
c
c ### transform updated stress to global from local ###
c
      if (IT2D.le.1) then
c
        call jancae_mv (STRESS,qtt,s2,ngens,ngens)
c 
        tarray(3)=STRESS(3)
        tarray(4)=STRESS(4)
        STRESS(3)=tarray(4)
        STRESS(4)=tarray(3)
      else
c
        call jancae_mv (STRESS,qtt2,s2,ngens,ngens)
c
      endif

c
c ### transform consitent matrix to global from local ###
c
      do i=1,ngens
      do j=1,ngens
         D(i,j)=0.0d0
      enddo
      enddo
c
      if (IT2D.le.1) then
c
        do i=1,ngens
        do j=1,ngens
        do k=1,ngens
        do l=1,ngens
           D(i,j)=D(i,j)+qtt(i,k)*dlocal(k,l)*qt(l,j)
        enddo
        enddo
        enddo
        enddo
c
        tarray(5)=D(1,3)
        tarray(6)=D(1,4)
        D(1,3)=tarray(6)
        D(1,4)=tarray(5)
        tarray(7)=D(2,3)
        tarray(8)=D(2,4)
        D(2,3)=tarray(8)
        D(2,4)=tarray(7)
        tarray(9)=D(3,1)
        tarray(10)=D(4,1)
        D(3,1)=tarray(10)
        D(4,1)=tarray(9)
        tarray(11)=D(3,2)
        tarray(12)=D(4,2)
        D(3,2)=tarray(12)
        D(4,2)=tarray(11)
        tarray(13)=D(3,3)
        tarray(14)=D(4,4)
        D(3,3)=tarray(14)
        D(4,4)=tarray(13)
        tarray(15)=D(3,4)
        tarray(16)=D(4,3)
        D(3,4)=tarray(16)
        D(4,3)=tarray(15)
c
      else
c
        do i=1,ngens
        do j=1,ngens
        do k=1,ngens
        do l=1,ngens
           D(i,j)=D(i,j)+qtt2(i,k)*dlocal2(k,l)*qt2(l,j)
        enddo
        enddo
        enddo
        enddo
c
      endif
c
c ### store consitent matrix in ARRAY ###
c
c      kdn=(LGTH3+4)+1
c      do i=1,ngens
c      do j=1,ngens
c         ARRAY(kdn)=D(i,j)
c         kdn=kdn+1
c      enddo
c      enddo
c
c ### PARAMETER:DPSP don't need this subroutine.
c
      do i=1,ngens
         DPSP(i)=0.0D0
      enddo
c
c
  100 continue
c
      RETURN
c
c ---------------------------- K E Y = 3 -------------------------------
c
    3 CONTINUE
c
c
c ### substitute ARRAY for D(i,j) ###
c
c      kdn=(LGTH3+4)+1
c      do i=1,ngens
c      do j=1,ngens
c         D(i,j)=ARRAY(kdn)
c         kdn=kdn+1
c      enddo
c      enddo
c
      RETURN
C
C ---------------------------- K E Y = 4 -------------------------------
C
    4 CONTINUE
c
      RETURN                                           
C
      END
c
c################################################################
c
c###############################
c########### from 
c###############################
c
c-------------------------------------------------------------
c     set internal state variables profile
c
      subroutine jancae_isvprof ( isvrsvd,isvsclr )
c-------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c
c               no reserved variables
       isvrsvd=0
c               statev(1) is for eq.plast.strain
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
      implicit real*8 (a-h,o-z)
      common /jancae1/ne,ip,lay
c
      write (6,*) 'ne= ',ne
      write (6,*) 'ip= ',ip
      write (6,*) 'lay=',lay
      write (6,*) 'exit code :',nexit
c     stop
c      call end  ! ht180310
c
      return
      end
c
c
c
c

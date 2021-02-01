      SUBROUTINE UCMAT3 [DLLEXPORT] (NG,NEL,IPT,IDEATH,ITE,IUPDTL,
     1                   STRESS,EPS,STRAIN,DEPS,DEPST,THSTR1,THSTR2,
     2                   KTR,INTER,SCP,ARRAY,LGTH1,IARRAY,LGTH2,LGTH3,
     3                   LGTH4,D,ALFA,CTD,ALFAA,CTDD,CTI,TMP1,TMP2,
     4                   TIME,ETIMV,ETIMV2,DT,PHIST,PRST,RN,PHIST1,
     5                   DPSP,TGRAD,INTEG,ISUBM,INDNL,LCTI,NSYMST,IERR,
     6                   DP,NELP,DPJE1D,DPJE2D,AKAPPA,PBAR,
     7                   NNODE,NODNUM,XYZ,DCA,OMEGAJ,WVTEMP,IOUT,KEY)  !ty150911
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
c     15/09/11
c
c
c  # EDIT BY
c     T.Yamanashi
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
c  # ARRAY(1): equivalent plastic strain (mod ht180310)
c       - store for next step
c
c  # ARRAY(2)...(7): plastic strain comp.  (mod ht180310)
c       - store for next step
c         components defined with local coordinate
c
c  # ARRAY(8)...(13): stress comp.
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
      parameter (mxpbs=10)     ! ht150908
      parameter (mxprop=100)   ! ht150908
c
      DIMENSION STRESS(6),DEPS(6),D(6,6),DCA(3,2),DPSP(6)
      DIMENSION ARRAY(*),IARRAY(*)
      DIMENSION CTI(99)
c
      DIMENSION s1(6),s2(6),de(6),dpe(6),dca3(3,3)
      DIMENSION x1(mxpbs,6),x2(mxpbs,6)   ! ht150908
      DIMENSION qt(6,6),qtt(6,6),dlocal(6,6)
      DIMENSION prop(mxprop)    ! ht150908
c
      ne =NEL
      ip =IPT
      lay=1
c
      nprop=mxprop
c
      ndi=3
      nshear=3
      ngens=6
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
c      IF (KEY.EQ.0) THEN
c          LGTH3=13
c          LGTH1=LGTH3+36
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
      do i=1,2
      do j=1,3
         dca3(j,i)=DCA(j,i)
      enddo
      enddo
      dca3(1,3)=DCA(2,1)*DCA(3,2)-DCA(3,1)*DCA(2,2)
      dca3(2,3)=DCA(3,1)*DCA(1,2)-DCA(1,1)*DCA(3,2)
      dca3(3,3)=DCA(1,1)*DCA(2,2)-DCA(2,1)*DCA(1,2)
c
      qt(1,1)=dca3(1,1)*dca3(1,1)
      qt(1,2)=dca3(2,1)*dca3(2,1)
      qt(1,3)=dca3(3,1)*dca3(3,1)
      qt(1,4)=dca3(1,1)*dca3(2,1)
      qt(1,5)=dca3(1,1)*dca3(3,1)
      qt(1,6)=dca3(2,1)*dca3(3,1)
      qt(2,1)=dca3(1,2)*dca3(1,2)
      qt(2,2)=dca3(2,2)*dca3(2,2)
      qt(2,3)=dca3(3,2)*dca3(3,2)
      qt(2,4)=dca3(1,2)*dca3(2,2)
      qt(2,5)=dca3(1,2)*dca3(3,2)
      qt(2,6)=dca3(2,2)*dca3(3,2)
      qt(3,1)=dca3(1,3)*dca3(1,3)
      qt(3,2)=dca3(2,3)*dca3(2,3)
      qt(3,3)=dca3(3,3)*dca3(3,3)
      qt(3,4)=dca3(1,3)*dca3(2,3)
      qt(3,5)=dca3(1,3)*dca3(3,3)
      qt(3,6)=dca3(2,3)*dca3(3,3)
      qt(4,1)=dca3(1,1)*dca3(1,2)+dca3(1,1)*dca3(1,2)
      qt(4,2)=dca3(2,1)*dca3(2,2)+dca3(2,1)*dca3(2,2)
      qt(4,3)=dca3(3,1)*dca3(3,2)+dca3(3,1)*dca3(3,2)
      qt(4,4)=dca3(1,1)*dca3(2,2)+dca3(1,2)*dca3(2,1)
      qt(4,5)=dca3(1,1)*dca3(3,2)+dca3(1,2)*dca3(3,1)
      qt(4,6)=dca3(2,1)*dca3(3,2)+dca3(2,2)*dca3(3,1)
      qt(5,1)=dca3(1,3)*dca3(1,1)+dca3(1,3)*dca3(1,1)
      qt(5,2)=dca3(2,3)*dca3(2,1)+dca3(2,3)*dca3(2,1)
      qt(5,3)=dca3(3,3)*dca3(3,1)+dca3(3,3)*dca3(3,1)
      qt(5,4)=dca3(1,3)*dca3(2,1)+dca3(1,1)*dca3(2,3)
      qt(5,5)=dca3(1,3)*dca3(3,1)+dca3(1,1)*dca3(3,3)
      qt(5,6)=dca3(2,3)*dca3(3,1)+dca3(2,1)*dca3(3,3)
      qt(6,1)=dca3(1,2)*dca3(1,3)+dca3(1,2)*dca3(1,3)
      qt(6,2)=dca3(2,2)*dca3(2,3)+dca3(2,2)*dca3(2,3)
      qt(6,3)=dca3(3,2)*dca3(3,3)+dca3(3,2)*dca3(3,3)
      qt(6,4)=dca3(1,2)*dca3(2,3)+dca3(1,3)*dca3(2,2)
      qt(6,5)=dca3(1,2)*dca3(3,3)+dca3(1,3)*dca3(3,2)
      qt(6,6)=dca3(2,2)*dca3(3,3)+dca3(2,3)*dca3(3,2)
c
      do i=1,6
      do j=1,6
         qtt(j,i)=qt(i,j)
      enddo
      enddo
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
      do i=1,6
         s1(i)=ARRAY(i+ioffset)
      enddo
c
c------ ht171230
c
      do i=1,LCTI
        prop(i)=CTI(i)
      enddo
      call jancae_prop_dim ( prop,nprop,
     &                       ndela,ndyld,ndihd,ndkin,
     &                       npbs )
c
c------ ht171230
c
      call jancae_mv (de,QT,DEPS,6,6)
c
      mjac=1
c
      call jancae_plasticity ( s1,s2,de,
     &                         ARRAY(1),dpp,dpe,de33,
     &                         x1,x2,mxpbs,
     &                         dlocal,
     &                         ndi,nshear,ngens,
     &                         nvbs,mjac,
     &                         prop,nprop )    ! ht180310
c
c
c ### store updated eq. plastic strain in ARRAY ###
c
      ARRAY(1)=ARRAY(1)+dpp  ! ht180310
c
c ### store updated plastic strain in ARRAY ### 
c
      do i=1,6
         ARRAY(1+I)=ARRAY(1+I)+dpe(I)  ! ht180310
      enddo
c
c ### store updated stress(local) in ARRAY ###
c
      do i=1,6
         ARRAY(i+7)=s2(i)
      enddo
c
c ### transform updated stress to global from local ###
c
      call jancae_mv (STRESS,qtt,s2,6,6)
c
c ### transform consitent matrix to global from local ###
c
      do i=1,6
      do j=1,6
         D(i,j)=0.0d0
      enddo
      enddo
c
      do i=1,6
      do j=1,6
      do k=1,6
      do l=1,6
         D(i,j)=D(i,j)+qtt(i,k)*dlocal(k,l)*qt(l,j)
      enddo
      enddo
      enddo
      enddo
c
c      enddo
c
c ### PARAMETER:DPSP don't need this subroutine.
c
      do i=1,6
         DPSP(I)=0
      enddo
c
c
  100 continue
c
      RETURN
C
C ---------------------------- K E Y = 3 -------------------------------
C
    3 CONTINUE
      RETURN
C
C ---------------------------- K E Y = 4 -------------------------------
C
    4 CONTINUE
      RETURN                                           
C
      END
c
c################################################################
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

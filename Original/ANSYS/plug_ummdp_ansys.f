c******** dummy *********
      subroutine usermat3d ()
        write (6,*) 'stop in usermat3d'
        call jancae_exit ( 9000 )
        return
      end
c******** dummy *********
      subroutine usermatps ()
        write (6,*) 'stop in usermatps'
        call jancae_exit ( 9000 )
        return
      end
c******** dummy *********
      subroutine usermatbm ()
        write (6,*) 'stop in usermatbm'
        call jancae_exit ( 9000 )
        return
      end
c******** dummy *********
      subroutine usermat1d ()
        write (6,*) 'stop in usermat1d'
        call jancae_exit ( 9000 )
        return
      end
c
c
c
c
c
*deck,usermat      USERDISTRIB  parallel                                gal
      subroutine usermat(
     &                   matId, elemId,kDomIntPt, kLayer, kSectPt,
     &                   ldstep,isubst,keycut,
     &                   nDirect,nShear,ncomp,nStatev,nProp,
     &                   Time,dTime,Temp,dTemp,
     &                   stress,ustatev,dsdePl,sedEl,sedPl,epseq,
     &                   Strain,dStrain, epsPl, prop, coords, 
     &                   var0, defGrad_t, defGrad,
     &                   tsstif, epsZZ,cutFactor,
     &                   var1, var2, var3, var4, var5,
     &                   var6, var7 )
c*************************************************************************
c     *** primary function ***
c
c           user defined material constitutive model
c
c      Attention:
c           User must define material constitutive law properly
c           according to the stress state such as 3D, plane strain
c           and axisymmetry, plane stress and 3D/1D beam.
c
c           A 3D material constitutive model can be used for
c           plane strain and axisymmetry cases.
c
c           When using shell elements, a plane stress algorithm
c           must be used.
c
c                                             gal July, 1999
c
c       The following demonstrates a USERMAT subroutine for 
c       a plasticity model, which is the same as TB, BISO,
c       for different stress states. 
c       See "ANSYS user material subroutine USERMAT" for detailed
c       description of how to write a USERMAT routine.
c
c       This routine calls four routines,
c       usermat3d.F, usermatps.F usermatbm.F and usermat1d.F, w.r.t.
c       the corresponding stress states.
c       Each routine can be also a usermat routine for the specific 
c       element.
c
c*************************************************************************
c Copyright ANSYS.  All Rights Reserved.
c
c     input arguments
c     ===============
c      matId     (int,sc,i)               material #
c      elemId    (int,sc,i)               element #
c      kDomIntPt (int,sc,i)               "k"th domain integration point
c      kLayer    (int,sc,i)               "k"th layer
c      kSectPt   (int,sc,i)               "k"th Section point
c      ldstep    (int,sc,i)               load step number
c      isubst    (int,sc,i)               substep number
c      nDirect   (int,sc,in)              # of direct components
c      nShear    (int,sc,in)              # of shear components
c      ncomp     (int,sc,in)              nDirect + nShear
c      nstatev   (int,sc,i)               Number of state variables
c      nProp     (int,sc,i)               Number of material constants
c
c      Temp      (dp,sc,in)               temperature at beginning of
c                                         time increment
c      dTemp     (dp,sc,in)               temperature increment 
c      Time      (dp,sc,in)               time at beginning of increment (t)
c      dTime     (dp,sc,in)               current time increment (dt)
c
c      Strain   (dp,ar(ncomp),i)          Strain at beginning of time increment
c      dStrain  (dp,ar(ncomp),i)          Strain increment
c      prop     (dp,ar(nprop),i)          Material constants defined by TB,USER
c      coords   (dp,ar(3),i)              current coordinates
c      defGrad_t(dp,ar(3,3),i)            Deformation gradient at time t
c      defGrad  (dp,ar(3,3),i)            Deformation gradient at time t+dt
c
c     input output arguments              
c     ======================             
c      stress   (dp,ar(ncomp),io)         stress
c      ustatev   (dp,ar(nstatev),io)      user state variables
c      sedEl    (dp,sc,io)                elastic work
c      sedPl    (dp,sc,io)                plastic work
c      epseq    (dp,sc,io)                equivalent plastic strain
c      epsPl   (dp,ar(ncomp),io)          plastic strain
c      var?     (dp,sc,io)                not used, they are reserved arguments 
c                                         for further development
c
c     output arguments
c     ================
c      keycut   (int,sc,o)                loading bisect/cut control
c                                         0 - no bisect/cut
c                                         1 - bisect/cut 
c                                         (factor will be determined by solution control)
c      dsdePl   (dp,ar(ncomp,ncomp),o)    material jacobian matrix
c      tsstif   (dp,ar(2),o)              transverse shear stiffness
c                                         tsstif(1) - Gxz
c                                         tsstif(2) - Gyz
c                                         tsstif(1) is also used to calculate hourglass
c                                         stiffness, this value must be defined when low
c                                         order element, such as 181, 182, 185 with uniform 
c                                         integration is used.
c      epsZZ    (dp,sc,o)                 strain epsZZ for plane stress,
c                                         define it when accounting for thickness change
c                                         in shell and plane stress states
c      cutFactor(dp,sc,o)                 time step size cut-back factor

c                                         define it if a smaller step size is wished

c                                         recommended value is 0~1c
c
c*************************************************************************
c
c      ncomp   6   for 3D  (nshear=3)
c      ncomp   4   for plane strain or axisymmetric (nShear = 1)
c      ncomp   3   for plane stress (nShear = 1)
c      ncomp   3   for 3d beam      (nShear = 2)
c      ncomp   1   for 1D (nShear = 0)
c
c      stresses and strains, plastic strain vectors
c          11, 22, 33, 12, 23, 13    for 3D
c          11, 22, 33, 12            for plane strain or axisymmetry
c          11, 22, 12                for plane stress
c          11, 13, 12                for 3d beam
c          11                        for 1D
c
c      material jacobian matrix
c        3D
c           dsdePl    |  1111   1122   1133   1112   1123   1113 |
c           dsdePl    |  2211   2222   2233   2212   2223   2213 |
c           dsdePl    |  3311   3322   3333   3312   3323   3313 |
c           dsdePl    |  1211   1222   1233   1212   1223   1213 |
c           dsdePl    |  2311   2322   2333   2312   2323   2313 |
c           dsdePl    |  1311   1322   1333   1312   1323   1313 |
c        plane strain or axisymmetric (11, 22, 33, 12)
c           dsdePl    |  1111   1122   1133   1112 |
c           dsdePl    |  2211   2222   2233   2212 |
c           dsdePl    |  3311   3322   3333   3312 |
c           dsdePl    |  1211   1222   1233   1212 |
c        plane stress (11, 22, 12)
c           dsdePl    |  1111   1122   1112 |
c           dsdePl    |  2211   2222   2212 |
c           dsdePl    |  1211   1222   1212 |
c        3d beam (11, 13, 12)
c           dsdePl    |  1111   1113   1112 |
c           dsdePl    |  1311   1313   1312 |
c           dsdePl    |  1211   1213   1212 |
c        1d
c           dsdePl    |  1111 |
c
c*************************************************************************
#include "impcom.inc"
c
      INTEGER          
     &                 matId, elemId,
     &                 kDomIntPt, kLayer, kSectPt,
     &                 ldstep,isubst,keycut,
     &                 nDirect,nShear,ncomp,nStatev,nProp
      DOUBLE PRECISION 
     &                 Time,    dTime,   Temp,    dTemp,
     &                 sedEl,   sedPl,   epseq,   epsZZ, cutFactor
      DOUBLE PRECISION 
     &                 stress  (ncomp  ), ustatev (nStatev),
     &                 dsdePl  (ncomp,ncomp),
     &                 Strain  (ncomp  ), dStrain (ncomp  ), 
     &                 epsPl   (ncomp  ), prop    (nProp  ), 
     &                 coords  (3),       
     &                 defGrad (3,3),     defGrad_t(3,3),
     &                 tsstif  (2)
c
c     EXTERNAL         usermat3d, usermatps, usermatbm, usermat1d
c
      DOUBLE PRECISION var0, var1, var2, var3, var4, var5,
     &                 var6, var7
c
c
c***************** User defined part *************************************
c
      INTEGER ne,ip,lay,nvbs,isvrsvd,isvsclr,
     &        mjac,i,is,n,
     &        mxpbs,npbs,mxprop,jprop,
     &        ndela,ndyld,ndihd,ndkin
c
      DOUBLE PRECISION p,pe,x1,s2,dp,dpe,x2,depsZZ,rprop,
     &                 young,poas
c
      common /jancae1/ne,ip,lay
c
      parameter (mxpbs=10)
      parameter (mxprop=100)
      dimension s2(ncomp),dpe(ncomp),
     &          x1(mxpbs,ncomp),x2(mxpbs,ncomp),
     &          pe(ncomp)
      dimension rprop(mxprop)
c
c
c                                                --- check element types 
      if ( ((ncomp.eq.3).and.(nDirect.ne.2)).or.(ncomp.lt.3) ) then
        write (6,*) 'UMMDp does NOT support beam elements',
     &               ncomp,nDirect,nShear
        write (6,*) 'stop in USERMAT'
        call jancae_exit ( 9000 )
      endif
c
c                                               --- variable name change
c                     ne  : element no.
c                     ip  : integration point no. 
c                     lay : layer no. of shell
c
c     in ANSYS ip :=  (ip_in_plane*-1)*maxip_thickness+ip_thickness
c
      ne =elemId
      ip =kDomIntPt
c     lay=kLayer
      lay=kSectPt
c
      jprop=mxprop
c                                         --- set debug and verbose mode
      call jancae_debugmode ( nvbs )
c
c                                        --- output detailed information
      if ( nvbs.ge.4 ) then
         call jancae_printinfo ( isubst,nDirect,nShear )
         call jancae_printinout ( 0,stress,dStrain,dsdePl,ncomp,
     &                          ustatev,nstatev )
      endif
c
c                                ---- set material properties from tbdata input
      do i=1,nProp
        rprop(i)=prop(i)
      enddo
      call jancae_prop_dim ( rprop,jprop,
     &                       ndela,ndyld,ndihd,ndkin,
     &                       npbs )
c     call jancae_set_prop ( rprop,jprop )
c     npbs=nint(rprop(1))                      ! number of partial back stresses
c     if ( npbs.gt.mxpbs ) then
c       write (6,*) 'npbs > mxpbs error in umat'
c       write (6,*) 'npbs =',npbs
c       write (6,*) 'mxpbs=',mxpbs
c       call jancae_exit ( 9000 )
c     endif
c
c                                                        -- check nstatv
      call jancae_check_nisv ( nstatev,ncomp,npbs )
c                                  copy current internal state variables
      call jancae_isvprof ( isvrsvd,isvsclr )
      call jancae_isv2pex ( isvrsvd,isvsclr,
     &                      ustatev,nstatev,
     &                      p,pe,x1,ncomp,mxpbs,npbs )
c                                                      *** update stress
c                                                and set tangent modulus
      mjac=1
c     mjac=-1
      call jancae_plasticity ( stress,s2,dStrain,
     &                         p,dp,dpe,depsZZ,
     &                         x1,x2,mxpbs,
     &                         dsdePl,
     &                         nDirect,nShear,ncomp,
     &                         nvbs,mjac,
     &                         rprop,jprop   )
c
c     *** for ANSYS special start ***
      keycut=0
      cutFactor=0.0
c                                  ----- update elastic and plastic work
      do i=1,ncomp
        sedEl=sedEl+0.5d0*(stress(i)+s2(i))*(dStrain(i)-dpe(i))
        sedPl=sedPl+0.5d0*(stress(i)+s2(i))*            dpe(i)
      enddo
c                                       ----- transverse shear stiffness
      if ( (nDirect.eq.2).and.(ncomp.eq.3) ) then
        if ( nint(rprop(1)).eq.0 ) then
          young=rprop(1+1)
          poas =rprop(1+2)
          tsstif(1)=0.5d0*young/(1.0d0+poas)
          tsstif(2)=tsstif(1)
        else
          write (6,*) 'USERMAT unknown elatic model :',nint(rprop(1))
          write (6,*) 'tsstif(trans. shear stiff.) must be defined.'
          call jancae_exit ( 9000 )
        endif
      endif
c     *** for ANSYS special end ***
c
c                                                     ---- update stress
      do i=1,ncomp
         stress(i)=s2(i)
      enddo
c                                            ---- update eq.plast,strain
      ustatev(isvrsvd+1)=p+dp
      epseq=p+dp
c                                         ---- update plast.strain comp.
      do i=1,ncomp
         is=isvrsvd+isvsclr+i
         ustatev(is)=ustatev(is)+dpe(i)
      enddo
c
      do i=1,ncomp
        is=isvrsvd+isvsclr+i
        epsPl(i)=ustatev(is)
      enddo
c                                                       --- update epsZZ
      ustatev(isvrsvd+2) = ustatev(isvrsvd+2) + depsZZ
      epsZZ = ustatev(isvrsvd+2)
c                                        --- update of back stress comp.
      if ( npbs.ne.0 ) then
         do n=1,npbs
            do i=1,ncomp
               is=isvrsvd+isvsclr+ncomp*n+i
               ustatev(is)=x2(n,i)
            enddo
         enddo
      endif
c                            ---- if debug mode, output return arguments
      if ( nvbs.ge.4 ) then
         call jancae_printinout (1,stress,dStrain,dsdePl,ncomp,
     &                           ustatev,nstatev )
      endif
c
      return
      end
c
c
c-------------------------------------------------------------
c     set internal state variables profile
c
      subroutine jancae_isvprof ( isvrsvd,isvsclr )
c-------------------------------------------------------------
c
c               no reserved variables
      isvrsvd=0
c               ustatev(1) is for eq.plast.strain
c               ustatev(2) is for strain in z-direction
      isvsclr=2
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
c
      common /jancae1/ne,ip,lay
c                                  nexit : exit code
      write (6,*) 'error code :',nexit
      write (6,*) 'element no.           :',ne
      write (6,*) 'integration point no. :',ip
      write (6,*) 'layer no.             :',lay
c
      call systop (0)
c
      return
      end
c
c
c
c

       subroutine umat(stress, statev, ddsdde, sse, spd, scd, rpl,
     1  ddsddt, drplde, drpldt, stran, dstran, time, dtime, temp,
     2  dtemp, predef, dpred, cmname, ndi, nshr, ntens, nstatv, props,
     3  nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, noel,
     4  npt, layer, kspt, kstep, kinc)
c
      include 'aba_param.inc'
c
      character*8 cmname
c
      dimension stress(ntens), statev(nstatv), ddsdde(ntens, ntens),
     1 ddsddt(ntens), drplde(ntens), stran(ntens), dstran(ntens),
     2 predef(1), dpred(1), props(nprops), coords(3), drot(3, 3),
     3 dfgrd0(3, 3), dfgrd1(3, 3)

c    local arrays
c ----------------------------------------------------------------
c    eelas  - elastic strains
c    eplas  - plastic strains
c    flow   - direction of plastic flow
c ----------------------------------------------------------------
c
      dimension eelas(6),eplas(6),flow(6), hard(3)
c
      parameter(zero=0.d0, one=1.d0, two=2.d0, three=3.d0, six=6.d0,
     1          enumax=.4999d0, newton=10, toler=1.0d-6)
c
c ----------------------------------------------------------------
c    umat for isotropic elasticity and isotropic mises plasticity
c    cannot be used for plane stress
c ----------------------------------------------------------------
c    props(1) - e
c    props(2) - nu
c    props(3..) - syield an hardening data
c    calls hardsub for curve of yield stress vs. plastic strain
c ----------------------------------------------------------------
c
c    elastic properties
c
      emod=props(1)
      enu=min(props(2), enumax)
      ebulk3=emod/(one-two*enu)
      eg2=emod/(one+enu)
      eg=eg2/two
      eg3=three*eg
      elam=(ebulk3-eg2)/three
c
c    elastic stiffness
c

      do k1=1, ndi
        do k2=1, ndi
          ddsdde(k2, k1)=elam
        end do
        ddsdde(k1, k1)=eg2+elam
      end do
      do k1=ndi+1, ntens
        ddsdde(k1, k1)=eg
      end do

c    recover elastic and plastic strains and rotate forward
c    also recover equivalent plastic strain
c
      call rotsig(statev(      1), drot, eelas, 2, ndi, nshr)
      call rotsig(statev(ntens+1), drot, eplas, 2, ndi, nshr)
      eqplas=statev(1+2*ntens)
c
c    calculate predictor stress and elastic strain
c
      do k1=1, ntens
        do k2=1, ntens
          stress(k2)=stress(k2)+ddsdde(k2, k1)*dstran(k1)
        end do
        eelas(k1)=eelas(k1)+dstran(k1)
      end do
c
c    calculate equivalent von mises stress
c
      smises=(stress(1)-stress(2))**2+(stress(2)-stress(3))**2
     1                               +(stress(3)-stress(1))**2
      do k1=ndi+1,ntens
        smises=smises+six*stress(k1)**2
      end do
      smises=sqrt(smises/two)
c
c    get yield stress from the specified hardening curve
c
      nvalue=nprops/2-1
      call uhard(syiel0, hard, eqplas, eqplasrt,time,dtime,temp,
     1     dtemp,noel,npt,layer,kspt,kstep,kinc,
     2     cmname,nstatv,statev,numfieldv,
     3     predef,dpred,nvalue,props(3))
c
c    determine if actively yielding
c
      if (smises.gt.(one+toler)*syiel0) then
c
c      actively yielding
c      separate the hydrostatic from the deviatoric stress
c      calculate the flow direction
c
        shydro=(stress(1)+stress(2)+stress(3))/three
        do k1=1,ndi
          flow(k1)=(stress(k1)-shydro)/smises
        end do
        do k1=ndi+1, ntens
          flow(k1)=stress(k1)/smises
        end do

c      solve for equivalent von mises stress
c      and equivalent plastic strain increment using newton iteration
c
        syield=syiel0
        deqpl=zero
        do kewton=1, newton
          rhs=smises-eg3*deqpl-syield
          deqpl=deqpl+rhs/(eg3+hard(1))
          call uhard(syield,hard,eqplas+deqpl,eqplasrt,time,dtime,temp,
     1     dtemp,noel,npt,layer,kspt,kstep,kinc,
     2     cmname,nstatv,statev,numfieldv,
     3     predef,dpred,nvalue,props(3))
          if(abs(rhs).lt.toler*syiel0) goto 10
        end do
c
c      write warning message to the .msg file
c
        write(7,2) newton
    2     format(//,30x,'***warning - plasticity algorithm did not ',
     1                  'converge after ',i3,' iterations')
   10   continue

c      update stress, elastic and plastic strains and
c      equivalent plastic strain
c
        do k1=1,ndi
          stress(k1)=flow(k1)*syield+shydro
          eplas(k1)=eplas(k1)+three/two*flow(k1)*deqpl
          eelas(k1)=eelas(k1)-three/two*flow(k1)*deqpl
        end do
        do k1=ndi+1,ntens
          stress(k1)=flow(k1)*syield
          eplas(k1)=eplas(k1)+three*flow(k1)*deqpl
          eelas(k1)=eelas(k1)-three*flow(k1)*deqpl
        end do
        eqplas=eqplas+deqpl
c
c      calculate plastic dissipation
c
        spd=deqpl*(syiel0+syield)/two

c      formulate the jacobian (material tangent)
c      first calculate effective moduli
c
        effg=eg*syield/smises
        effg2=two*effg
        effg3=three/two*effg2
        efflam=(ebulk3-effg2)/three
        effhrd=eg3*hard(1)/(eg3+hard(1))-effg3
c...
      IF (PROPS(7).LT..001) GO TO 99
c...
        do k1=1, ndi
          do k2=1, ndi
            ddsdde(k2, k1)=efflam
          end do
          ddsdde(k1, k1)=effg2+efflam
        end do
        do k1=ndi+1, ntens
          ddsdde(k1, k1)=effg
        end do
        do k1=1, ntens
          do k2=1, ntens
            ddsdde(k2, k1)=ddsdde(k2, k1)+effhrd*flow(k2)*flow(k1)
          end do
        end do
c...
   99 CONTINUE
c...
      endif
c
c    store elastic and (equivalent) plastic strains
c    in state variable array
c
      do k1=1, ntens
        statev(k1)=eelas(k1)
        statev(k1+ntens)=eplas(k1)
      end do
      statev(1+2*ntens)=eqplas
c
      return
      end


      subroutine uhard(syield,hard,eqplas,eqplasrt,time,dtime,temp,
     1     dtemp,noel,npt,layer,kspt,kstep,kinc,
     2     cmname,nstatv,statev,numfieldv,
     3     predef,dpred,nvalue,table)

      include 'aba_param.inc'

      character*80 cmname
      dimension hard(3),statev(nstatv),time(*),
     1          predef(numfieldv),dpred(*)
c
      dimension table(2, nvalue)
c
      parameter(zero=0.d0)
c
c    set yield stress to last value of table, hardening to zero
c
      syield=table(1, nvalue)
      hard(1)=zero
c    if more than one entry, search table
c
      if(nvalue.gt.1) then
        do k1=1, nvalue-1
          eqpl1=table(2,k1+1)
          if(eqplas.lt.eqpl1) then
            eqpl0=table(2, k1)
            if(eqpl1.le.eqpl0) then
              write(7, 1)
    1         format(//, 30x, '***error - plastic strain must be `,
     1                        `entered in ascending order')
              call xit
            endif
c
c          current yield stress and hardening
c
            deqpl=eqpl1-eqpl0
            syiel0=table(1, k1)
            syiel1=table(1, k1+1)
            dsyiel=syiel1-syiel0
            hard(1)=dsyiel/deqpl
            syield=syiel0+(eqplas-eqpl0)*hard(1)
            goto 10
          endif
        end do
   10   continue
      endif
      return
      end

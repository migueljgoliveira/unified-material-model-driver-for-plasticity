c-----------------------------------------------(yld2004-18p)
c     Barlat YLD2004-18p yield function and its dfferentials
c     ( IJP v.21(2005) p1009-1039. )
c
      subroutine jancae_yld2004_18p ( s,se,dseds,d2seds2,nreq,
     &                                pryld,ndyld )
c------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension s(6),dseds(6),d2seds2(6,6),pryld(ndyld)
c
      dimension sp1(6),sp2(6),cp1(6,6),cp2(6,6),cl(6,6),ctp1(6,6)
      dimension ctp2(6,6),psp1(3),psp2(3),hp1(3),hp2(3)
c     dimension aap1(3),aap2(3),bbp1(3),bbp2(3)
      dimension dfadpsp1(3),dfadpsp2(3),dpsdhp1(3,3),dpsdhp2(3,3)
      dimension dfadhp1(3),dfadhp2(3)
      dimension dhdsp1(3,6),dhdsp2(3,6),dsdsp1(6,6), dsdsp2(6,6)
      dimension dfads(6)
      dimension d2fadpsp11(3,3),d2fadpsp22(3,3)
      dimension d2fadpsp12(3,3),d2fadpsp21(3,3)
      dimension d2fadhp11(3,3),d2fadhp22(3,3)
      dimension d2fadhp12(3,3),d2fadhp21(3,3)
      dimension d2psdhp11(3,3,3),d2psdhp22(3,3,3)
      dimension d2hdsp11(3,6,6),d2hdsp22(3,6,6)
      dimension d2fads2(6,6)
      dimension delta(3,3)
      dimension xx1(3,6),xx2(3,6)
c
c   variables 
c
c     sp1(6),sp2(6) :linear-transformed stress
c     cp1(6,6),cp2(6,6) : matrix for anisotropic parameters
c     cl(6,6)  : matrix for transforming Cauchy stress to deviator
c     ctp1(6,6): matrix for transforming Cauchy stress to sp1
c     ctp2(6,6): matrix for transforming Cauchy stress to sp2
c     dc    : coefficient of equivalent stress(se=(fai/dc)**(1/am))
c     aap1(3),aap2(3)  :for dc
c     ppp1,ppp2: for dc
c     qqp1,qqp2: for dc
c     ttp1,ttp2: for dc
c     bbp1(3),bbp2(3) : for dc
c     fai      : yield fuction
c     psp1(3)  : principal values of sp1
c     psp2(3)  : principal values of sp2
c     hp1(3)   : invariants of sp1
c     hp2(3)   : invariants of sp2
c     cep1,cep2 : coefficient of charactaristic equation p
c     ceq1,ceq2 : coefficient of charactaristic equation q
c     cet1,cet2 : arccos ( q/p^(3/2) )
c
c     dfadpsp1(3),dfadpsp2(3) : d(fai)/d(psp)
c     dfadhp1(3),dfadhp2(3) : d(fai)/d(hp)
c     dpsdhp1(3,3),dpsdhp2(3,3) : d(psp)/d(hp)
c     dhdsp1(3,6),dhdsp2(3,6) : d(hp)/d(sp)
c     dsdsp1(6,6), dsdsp2(6,6) : d(sp)/d(s)
c     dfads(6) : d(fai)/d(s)
c
c     d2fadpsp11(3,3),d2fadpsp22(3,3) : d2(fai)/d(psp)2
c     d2fadpsp12(3,3),d2fadpsp21(3,3) : d2(fai)/d(psp)2
c     d2psdhp11(3,3,3),d2psdhp22(3,3,3) :d2(psp)/d(hp)2
c     d2hdsp11(3,6,6),d2hdsp22(3,6,6) : d2(hp)/d(sp)2
c     d2fads2(6,6) : d2(fai)/d(s)2
c
c     eps2,eps3 : values for calculation of d(psp)/d(hp),d2(psp)/d(hp)2
c
      pi=acos(-1.0d0)
      eps2 = 1.0d-15
      eps3 = 1.0d-8
      del = 1.0d-4
c                                                 Kronecker Delta
      call jancae_clear2( delta,3,3 )
      do i=1,3
        delta(i,i)=1.0d0
      end do
c                                           set anisotropic parameters
      call jancae_clear2( cp1,6,6 )
      call jancae_clear2( cp2,6,6 )
      cp1(1,2) = -pryld(1+1)
      cp1(1,3) = -pryld(1+2)
      cp1(2,1) = -pryld(1+3)
      cp1(2,3) = -pryld(1+4)
      cp1(3,1) = -pryld(1+5)
      cp1(3,2) = -pryld(1+6)
      cp1(4,4) =  pryld(1+7)  ! for tau_xy (c'66 in original paper)
      cp1(5,5) =  pryld(1+8)  ! for tau_yz (c'44 in original paper)
      cp1(6,6) =  pryld(1+9)  ! for tau_zx (c'55 in original paper)
      cp2(1,2) = -pryld(1+10)
      cp2(1,3) = -pryld(1+11)
      cp2(2,1) = -pryld(1+12)
      cp2(2,3) = -pryld(1+13)
      cp2(3,1) = -pryld(1+14)
      cp2(3,2) = -pryld(1+15)
      cp2(4,4) =  pryld(1+16) ! for tau_xy (c"66 in original paper)
      cp2(5,5) =  pryld(1+17) ! for tau_yz (c"44 in original paper)
      cp2(6,6) =  pryld(1+18) ! for tau_zx (c"55 in original paper)
      am       =  pryld(1+19)
      ami=1.0d0/am
c
c                  set matrix for transforming Cauchy stress to deviator
      call jancae_clear2( cl,6,6 )
      do i=1,3
        do j=1,3
          if(i==j) then
            cl(i,j)=2.0d0
          else
            cl(i,j)=-1.0d0
          end if
        end do
      end do
      do i=4,6      
        cl(i,i) = 3.0d0
      end do
      do i=1,6
        do j=1,6
          cl(i,j)=cl(i,j)/3.0d0
        end do
      end do
c
c                       matrix for transforming Cauchy stress to sp1,sp2
      call jancae_mm (ctp1,cp1,cl,6,6,6)
      call jancae_mm (ctp2,cp2,cl,6,6,6)
c                                    coefficient of equivalent stress dc
      call jancae_yld2004_18p_coef (cp1,cp2,pi,am,dc)
c                                       calculation of equivalent stress
      call jancae_yld2004_18p_yf (ctp1,ctp2,s,am,ami,dc,pi,
     1       sp1,sp2,psp1,psp2,hp1,hp2,cetpq1,cetpq2,fai,se)
c
c                                                  1st order differential
      if ( nreq.ge.1 ) then
c                                                        -d(fai)/d(psp)
        do i = 1,3
          dfadpsp1(i) = 0.0d0
          dfadpsp2(i) = 0.0d0
          do j= 1,3
            dfadpsp1(i)=dfadpsp1(i)+(psp1(i)-psp2(j))
     1                *abs(psp1(i)-psp2(j))**(am-2.0d0)
            dfadpsp2(i)=dfadpsp2(i)+(psp1(j)-psp2(i))
     1                *abs(psp1(j)-psp2(i))**(am-2.0d0)
          end do
          dfadpsp1(i)=dfadpsp1(i)*am
          dfadpsp2(i)=dfadpsp2(i)*(-am)
        end do
c                                            - d(psp)/d(hp)&d(fai)/d(hp)
        call jancae_clear2( dpsdhp1,3,3 )
        call jancae_clear2( dpsdhp2,3,3 )
        call jancae_clear1( dfadhp1,3 )
        call jancae_clear1( dfadhp2,3 )
c                                                  -- theta'<>0 & <>pai
        if(abs(cetpq1-1.0d0)>=eps2 .and. abs(cetpq1+1.0d0)>=eps2) then
          do i=1,3
            call jancae_yld2004_18p_dpsdhp(i,psp1,hp1,dpsdhp1)
          end do
c                                                          -- theta'=0
        else if(abs(cetpq1-1.0d0)<eps2) then
          i=1
          call jancae_yld2004_18p_dpsdhp(i,psp1,hp1,dpsdhp1)
          do i=2,3
            do j=1,3
              dpsdhp1(i,j)=-0.5d0*(dpsdhp1(1,j)-3.0d0*delta(1,j))
            end do
           end do
c                                                         -- theta'=pai
        else
          i=3
          call jancae_yld2004_18p_dpsdhp(i,psp1,hp1,dpsdhp1)
          do i=1,2
            do j=1,3
              dpsdhp1(i,j)=-0.5d0*(dpsdhp1(3,j)-3.0d0*delta(1,j))
            end do
           end do
        end if
c                                                  -- theta''<>0 & <>pai
        if(abs(cetpq2-1.0d0)>=eps2 .and. abs(cetpq2+1.0d0)>=eps2) then
          do i=1,3
            call jancae_yld2004_18p_dpsdhp(i,psp2,hp2,dpsdhp2)
          end do
c                                                          -- theta''=0
        else if(abs(cetpq2-1.0d0)<eps2)then
          i=1
          call jancae_yld2004_18p_dpsdhp(i,psp2,hp2,dpsdhp2)
          do i=2,3
            do j=1,3
              dpsdhp2(i,j)=-0.5d0*(dpsdhp2(1,j)-3.0d0*delta(1,j))
            end do
           end do
c                                                          -- theta''=pai
        else
          i=3
          call jancae_yld2004_18p_dpsdhp(i,psp2,hp2,dpsdhp2)
          do i=1,2
            do j=1,3
              dpsdhp2(i,j)=-0.5d0*(dpsdhp2(3,j)-3.0d0*delta(1,j))
            end do
           end do
        end if 
c
          do i=1,3
            do j=1,3
              dfadhp1(i)=dfadhp1(i)+dfadpsp1(j)*dpsdhp1(j,i)
              dfadhp2(i)=dfadhp2(i)+dfadpsp2(j)*dpsdhp2(j,i)
            end do
          end do
c                                                          - d(hp)/d(sp)
        call jancae_clear2( dhdsp1,3,6 )
        call jancae_clear2( dhdsp2,3,6 )
        do i = 1,3
          j=mod(i,3)+1
          k=mod(i+1,3)+1
          l=mod(i,3)+4
          dhdsp1(1,i)=1.0d0/3.0d0
          dhdsp2(1,i)=1.0d0/3.0d0
          dhdsp1(2,i)=-1.0d0/3.0d0*(sp1(j)+sp1(k))
          dhdsp2(2,i)=-1.0d0/3.0d0*(sp2(j)+sp2(k))
          dhdsp1(3,i)=1.0d0/2.0d0*(sp1(j)*sp1(k)-sp1(l)**2)
          dhdsp2(3,i)=1.0d0/2.0d0*(sp2(j)*sp2(k)-sp2(l)**2)
        end do
        do i=4,6
          k=mod(i+1,3)+1
          l=mod(i,3)+4
          m=mod(i+1,3)+4
          dhdsp1(2,i)=2.0d0/3.0d0*sp1(i)
          dhdsp2(2,i)=2.0d0/3.0d0*sp2(i)
          dhdsp1(3,i)=sp1(l)*sp1(m)-sp1(k)*sp1(i)
          dhdsp2(3,i)=sp2(l)*sp2(m)-sp2(k)*sp2(i)
        end do 
c                                                            -d(sp)/d(s)
        do i=1,6
          do j=1,6
            dsdsp1(i,j)=ctp1(i,j)
            dsdsp2(i,j)=ctp2(i,j)
          end do
        end do
c                                                           -d(fai)/d(s)
        call jancae_clear1( dfads,6 )
        call jancae_clear2( xx1,3,6 )
        call jancae_clear2( xx2,3,6 )
        do l=1,6
          do j=1,3
            do k=1,6
              xx1(j,l)=xx1(j,l)+dhdsp1(j,k)*dsdsp1(k,l)
              xx2(j,l)=xx2(j,l)+dhdsp2(j,k)*dsdsp2(k,l)
            end do
            dfads(l)=dfads(l)
     1               +dfadhp1(j)*xx1(j,l)
     2               +dfadhp2(j)*xx2(j,l)            
          end do
        end do
c                                                           -d(se)/d(s)
        dsedfa=fai**(ami-1.0d0)/am/dc**ami
        do i=1,6
          dseds(i)=dsedfa*dfads(i)
        end do
c
      endif
c                                                 2nd order differential
      if ( nreq.ge.2 ) then
c                                                      -d2(fai)/d(psp)2
        call jancae_clear2( d2fadpsp11,3,3 )
        call jancae_clear2( d2fadpsp22,3,3 )
        call jancae_clear2( d2fadpsp12,3,3 )
        call jancae_clear2( d2fadpsp21,3,3 )
        do i=1,3
          d2fadpsp11(i,i)=am*(am-1.0d0)
     1                   *(abs(psp1(i)-psp2(1))**(am-2.0d0)
     2                   + abs(psp1(i)-psp2(2))**(am-2.0d0)
     3                   + abs(psp1(i)-psp2(3))**(am-2.0d0))
          d2fadpsp22(i,i)=am*(am-1.0d0)
     1                   *(abs(psp1(1)-psp2(i))**(am-2.0d0)
     2                   + abs(psp1(2)-psp2(i))**(am-2.0d0)
     3                   + abs(psp1(3)-psp2(i))**(am-2.0d0))
          do j=1,3
            d2fadpsp12(i,j)=-am*(am-1.0d0)
     1                     * abs(psp1(i)-psp2(j))**(am-2.0d0)
            d2fadpsp21(i,j)=-am*(am-1.0d0)
     1                     * abs(psp1(j)-psp2(i))**(am-2.0d0)
          end do
        end do
c                                                    - d2(psp)/d(hp)2 
        call jancae_clear3( d2psdhp11,3,3,3 )
        call jancae_clear3( d2psdhp22,3,3,3 )
c
        if(abs(cetpq1-1.0d0)>=eps3 .and. abs(cetpq2-1.0d0)>=eps3 .and.
     1     abs(cetpq1+1.0d0)>=eps3 .and. abs(cetpq2+1.0d0)>=eps3) then
c
          do i=1,3
            call jancae_yld2004_18p_d2psdhp(i,psp1,hp1,d2psdhp11)
            call jancae_yld2004_18p_d2psdhp(i,psp2,hp2,d2psdhp22)
          end do
c                                                     - d2(fai)/d(hp)2
          call jancae_clear2(d2fadhp11,3,3)
          call jancae_clear2(d2fadhp12,3,3)
          call jancae_clear2(d2fadhp21,3,3)
          call jancae_clear2(d2fadhp22,3,3)
c                                                 -- d2(fai)/d(hd)d(hd)
          do iq=1,3
            do m=1,3
              do ip=1,3
                do l=1,3
                  d2fadhp11(iq,m)=d2fadhp11(iq,m)
     1                   +d2fadpsp11(ip,l)*dpsdhp1(l,m)*dpsdhp1(ip,iq)
                end do
                d2fadhp11(iq,m)=d2fadhp11(iq,m)
     1                   +dfadpsp1(ip)*d2psdhp11(ip,iq,m)
              end do
            end do
          end do
c                                               -- d2(fai)/d(hdd)d(hdd)
          do iq=1,3
            do m=1,3
              do ip=1,3
                do l=1,3
                  d2fadhp22(iq,m)=d2fadhp22(iq,m)
     1                   +d2fadpsp22(ip,l)*dpsdhp2(l,m)*dpsdhp2(ip,iq)
                end do
                d2fadhp22(iq,m)=d2fadhp22(iq,m)
     1                   +dfadpsp2(ip)*d2psdhp22(ip,iq,m)
              end do
            end do
          end do
c                             -- d2(fai)/d(hdd)d(hd) & d2(fai)/d(hd)d(hdd)
          do iq=1,3
            do m=1,3
              do ip=1,3
                do l=1,3
                  d2fadhp12(iq,m)=d2fadhp12(iq,m)
     1                   +d2fadpsp12(ip,l)*dpsdhp2(l,m)*dpsdhp1(ip,iq)
                  d2fadhp21(iq,m)=d2fadhp21(iq,m)
     1                   +d2fadpsp21(ip,l)*dpsdhp1(l,m)*dpsdhp2(ip,iq)
                end do
               end do
            end do
          end do
c                                                          - d2(hp)/d(sp)2
          call jancae_clear3( d2hdsp11,3,6,6)
          call jancae_clear3( d2hdsp22,3,6,6)
          do i=1,3
            j=mod(i,3)+1
            k=mod(i+1,3)+1
            dummy=-1.0d0/3.0d0
            d2hdsp11(2,i,j)=dummy
            d2hdsp11(2,j,i)=dummy
            d2hdsp22(2,i,j)=dummy
            d2hdsp22(2,j,i)=dummy
            d2hdsp11(3,i,j)=sp1(k)/2.0d0
            d2hdsp11(3,j,i)=sp1(k)/2.0d0
            d2hdsp22(3,i,j)=sp2(k)/2.0d0
            d2hdsp22(3,j,i)=sp2(k)/2.0d0
          end do
          do i=4,6
            j=mod(i,3)+4
            k=mod(i+1,3)+1
            m=mod(i+1,3)+4
            dummy=2.0d0/3.0d0
            d2hdsp11(2,i,i)=dummy
            d2hdsp22(2,i,i)=dummy
            d2hdsp11(3,i,i)=-sp1(k)
            d2hdsp22(3,i,i)=-sp2(k)
            d2hdsp11(3,i,j)=sp1(m)
            d2hdsp11(3,j,i)=sp1(m)
            d2hdsp22(3,i,j)=sp2(m)
            d2hdsp22(3,j,i)=sp2(m)
          end do
          do i=1,3
            j=mod(i,3)+4
            d2hdsp11(3,i,j)=-sp1(j)
            d2hdsp11(3,j,i)=-sp1(j)
            d2hdsp22(3,i,j)=-sp2(j)
            d2hdsp22(3,j,i)=-sp2(j)
          end do
c                                                      -d2(fai)/d(s)2
          call jancae_clear2( d2fads2,6,6)
          call jancae_clear2( xx1,3,6 )
          call jancae_clear2( xx2,3,6 )
          do i=1,3
            do j=1,6
              do ip=1,6
                xx1(i,j)=xx1(i,j)+dhdsp1(i,ip)*dsdsp1(ip,j)
                xx2(i,j)=xx2(i,j)+dhdsp2(i,ip)*dsdsp2(ip,j)
              end do
            end do
          end do
          do i=1,6
            do j=1,6
              do iq=1,3
                do m=1,3
                  d2fads2(i,j)=d2fads2(i,j)
     1                       +d2fadhp11(iq,m)*xx1(iq,i)*xx1(m,j)
     2                       +d2fadhp12(iq,m)*xx1(iq,i)*xx2(m,j)
     3                       +d2fadhp21(iq,m)*xx2(iq,i)*xx1(m,j)
     4                       +d2fadhp22(iq,m)*xx2(iq,i)*xx2(m,j)
                end do
                do n=1,6
                  do ir=1,6
                    d2fads2(i,j)=d2fads2(i,j)
     1                       +d2hdsp11(iq,ir,n)*dfadhp1(iq)
     2                                  *dsdsp1(ir,i)*dsdsp1(n,j)
     3                       +d2hdsp22(iq,ir,n)*dfadhp2(iq)
     4                                  *dsdsp2(ir,i)*dsdsp2(n,j)
                  end do
                end do
              end do
            end do
          end do
c                                                          -d2(se)/d(s)2
          d2sedfa2=ami*(ami-1.0d0)*fai**(ami-2.0d0)/dc**ami
          do i=1,6
            do j=1,6
              d2seds2(i,j)=d2sedfa2*dfads(i)*dfads(j)
     1                  +dsedfa*d2fads2(i,j)
            end do
          end do
        else
c                                                 -numerical differential
          call jancae_yld2004_18p_nu2 ( ctp1,ctp2,s,se,am,ami,
     1                                  dc,pi,del,d2seds2)
        end if
      endif
c
      return
      end
c
c
c
c-----------------------------------------------(yld2004-18p)
c     calculate  coefficient of equivalent stress dc
c
      subroutine jancae_yld2004_18p_coef (cp1,cp2,pi,am,dc)
c------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension cp1(6,6),cp2(6,6),bbp1(3),bbp2(3)
c
      call jancae_yld2004_18p_coef_sub(cp1,pi,bbp1)
      call jancae_yld2004_18p_coef_sub(cp2,pi,bbp2)
      dc=0.0d0
      do i=1,3
        do j=1,3
          dc=dc+(abs(bbp1(i)-bbp2(j)))**am
        end do
      end do  
c
      return
      end
c
c-----------------------------------------------(yld2004-18p)
c     calculate  coefficient of equivalent stress dc 2
c
      subroutine jancae_yld2004_18p_coef_sub (cp,pi,bbp)
c------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension cp(6,6),aap(3),bbp(3)
c
c                                               - coefficients aap
      aap(1)=(cp(1,2)+cp(1,3)-2.0d0*cp(2,1)
     1        +cp(2,3)-2.0d0*cp(3,1)+cp(3,2))/9.0d0
      aap(2)=((2.0d0*cp(2,1)-cp(2,3))*(cp(3,2)-2.0d0*cp(3,1))
     2       +(2.0d0*cp(3,1)-cp(3,2))*(cp(1,2)+cp(1,3))
     3        +(cp(1,2)+cp(1,3))*(2.0d0*cp(2,1)-cp(2,3)))/2.7d1
      aap(3)=(cp(1,2)+cp(1,3))*(cp(2,3)-2.0d0*cp(2,1))
     1        *(cp(3,2)-2.0d0*cp(3,1))/5.4d1
c                                         - coefficients ppp,qqp,ttp
      ppp=aap(1)**2+aap(2)
      qqp=(2.0d0*aap(1)**3+3.0d0*aap(1)*aap(2)+2.0d0*aap(3))/2.0d0
      ttp=acos(qqp/ppp**(3.0d0/2.0d0))
c                                                  - coefficients bbp
      bbp(1)=2.0d0*sqrt(ppp)*cos(ttp/3.0d0)+aap(1)
      bbp(2)=2.0d0*sqrt(ppp)*cos((ttp+4.0d0*pi)/3.0d0)+aap(1)
      bbp(3)=2.0d0*sqrt(ppp)*cos((ttp+2.0d0*pi)/3.0d0)+aap(1)
c
      return
      end
c
c-----------------------------------------------(yld2004-18p)
c     calculate yield function
c
      subroutine jancae_yld2004_18p_yf (ctp1,ctp2,s,am,ami,dc,pi,
     1         sp1,sp2,psp1,psp2,hp1,hp2,cetpq1,cetpq2,fai,se)
c------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension ctp1(6,6),ctp2(6,6),s(6)
      dimension sp1(6),sp2(6),psp1(3),psp2(3),hp1(3),hp2(3)
c
      call jancae_yld2004_18p_yfsub(ctp1,s,pi,
     1                              sp1,psp1,hp1,cetpq1)
      call jancae_yld2004_18p_yfsub(ctp2,s,pi,
     1                              sp2,psp2,hp2,cetpq2)
c                                                         yield function
      fai=0.0d0
      do i=1,3
        do j=1,3
          fai=fai+(abs(psp1(i)-psp2(j)))**am
        end do
      end do
c                                                      equivalent stress
      se=(fai/dc)**ami
c
      return
      end
c
c
c-----------------------------------------------(yld2004-18p)
c     calculate yield function2
c
      subroutine jancae_yld2004_18p_yfsub(ctp,s,pi,
     1                                    sp,psp,hp,cetpq)
c------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension ctp(6,6),s(6)
      dimension sp(6),psp(3),hp(3)
c
c                                              linear-transformed stress
      call jancae_mv (sp,ctp,s,6,6)
c                                                       invariants of sp
      hp(1)=(sp(1)+sp(2)+sp(3))/3.0d0
      hp(2)=(sp(5)**2+sp(6)**2+sp(4)**2
     1       -sp(2)*sp(3)-sp(3)*sp(1)-sp(1)*sp(2))/3.0d0
      hp(3)=(2.0d0*sp(5)*sp(6)*sp(4)+sp(1)*sp(2)*sp(3)
     1       -sp(1)*sp(5)**2-sp(2)*sp(6)**2-sp(3)*sp(4)**2)/2.0d0
c                                coefficients of charactaristic equation
      hpq=sqrt(hp(1)**2+hp(2)**2+hp(3)**2)
      if ( hpq.gt.1.0e-16 ) then
        cep=hp(1)**2+hp(2)
        ceq=(2.0d0*hp(1)**3+3.0d0*hp(1)*hp(2)+2.0d0*hp(3))/2.0d0
        cetpq=ceq/cep**(3.0d0/2.0d0)
        if ( cetpq >  1.0d0 ) cetpq= 1.0d0
        if ( cetpq < -1.0d0 ) cetpq=-1.0d0
        cet=acos(cetpq)
c                                                principal values of sp1
        psp(1)=2.0d0*sqrt(cep)*cos(cet/3.0d0)+hp(1)
        psp(2)=2.0d0*sqrt(cep)*cos((cet+4.0d0*pi)/3.0d0)+hp(1)
        psp(3)=2.0d0*sqrt(cep)*cos((cet+2.0d0*pi)/3.0d0)+hp(1)
      else
        cetpq=0.0
        do i=1,3
           psp(i)=0.0
        enddo
      endif
c
      return
      end
c
c
c-----------------------------------------------(yld2004-18p)
c     numerical differential for 2nd order differentials
c
      subroutine jancae_yld2004_18p_nu2 (ctp1,ctp2,s,se,am,ami,
     1                                   dc,pi,del,d2seds2)
c------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension s(6),ctp1(6,6),ctp2(6,6),d2seds2(6,6)
      dimension sp1(6),sp2(6),psp1(3),psp2(3),hp1(3),hp2(3),s0(6)
c
      s0(:) = s(:)
      do i=1, 6
        do j=1, 6
          if(i == j) then
            s0(i)=s(i)-del
            call jancae_yld2004_18p_yf (ctp1,ctp2,s0,am,ami,dc,pi,
     1         sp1,sp2,psp1,psp2,hp1,hp2,cetpq1,cetpq2,fai,sea)
            s0(i)=s(i)+del
            call jancae_yld2004_18p_yf (ctp1,ctp2,s0,am,ami,dc,pi,
     1         sp1,sp2,psp1,psp2,hp1,hp2,cetpq1,cetpq2,fai,seb)
            s0(i)=s(i)
            abc1=(se-sea)/del
            abc2=(seb-se)/del
            d2seds2(i,j) = (abc2-abc1)/del
          else
            s0(i)=s(i)-del
            s0(j)=s(j)-del
            call jancae_yld2004_18p_yf (ctp1,ctp2,s0,am,ami,dc,pi,
     1         sp1,sp2,psp1,psp2,hp1,hp2,cetpq1,cetpq2,fai,seaa)
            s0(i)=s(i)+del
            s0(j)=s(j)-del
            call jancae_yld2004_18p_yf (ctp1,ctp2,s0,am,ami,dc,pi,
     1         sp1,sp2,psp1,psp2,hp1,hp2,cetpq1,cetpq2,fai,seba)
            s0(i)=s(i)-del
            s0(j)=s(j)+del
            call jancae_yld2004_18p_yf (ctp1,ctp2,s0,am,ami,dc,pi,
     1         sp1,sp2,psp1,psp2,hp1,hp2,cetpq1,cetpq2,fai,seab)
            s0(i)=s(i)+del
            s0(j)=s(j)+del
            call jancae_yld2004_18p_yf (ctp1,ctp2,s0,am,ami,dc,pi,
     1         sp1,sp2,psp1,psp2,hp1,hp2,cetpq1,cetpq2,fai,sebb)
            s0(i)=s(i)
            s0(j)=s(j)
            abc1 = (seba - seaa)/(2.0d0*del)
            abc2 = (sebb - seab)/(2.0d0*del)
            d2seds2(i,j) = (abc2 - abc1)/(2.0d0*del)
          end if
        end do
      end do
c
      return
      end
c
c
c-----------------------------------------------(yld2004-18p)
c     calculate d(psp)/d(hp)
c
      subroutine jancae_yld2004_18p_dpsdhp(i,psp,hp,dpsdhp)
c------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension psp(3),hp(3),dpsdhp(3,3)
c
      dummy=psp(i)**2-2.0d0*hp(1)*psp(i)-hp(2)
      dpsdhp(i,1)=psp(i)**2/dummy
      dpsdhp(i,2)=psp(i)/dummy
      dpsdhp(i,3)=2.0d0/3.0d0/dummy
c
      return
      end
c
c
c-----------------------------------------------(yld2004-18p)
c     calculate d2(psp)/d(hp)2
c
      subroutine jancae_yld2004_18p_d2psdhp(i,psp,hp,d2psdhp)
c------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension psp(3),hp(3),d2psdhp(3,3,3)
c
      dummy=(psp(i)**2-2.0d0*hp(1)*psp(i)-hp(2))**3
      d2psdhp(i,1,1)=2.0d0*psp(i)**3
     1        *(psp(i)**2-3.0d0*hp(1)*psp(i)-2.0d0*hp(2))/dummy
      d2psdhp(i,2,2)=-2.0d0*psp(i)*(hp(1)*psp(i)+hp(2))/dummy
      d2psdhp(i,3,3)=-8.0d0/9.0d0*(psp(i)-hp(1))/dummy
      d2psdhp(i,1,2)=psp(i)**2
     1        *(psp(i)**2-4.0d0*hp(1)*psp(i)-3.0d0*hp(2))/dummy
      d2psdhp(i,2,1)=d2psdhp(i,1,2)
      d2psdhp(i,2,3)=-2.0d0*(psp(i)**2+hp(2))/3.0d0/dummy
      d2psdhp(i,3,2)=d2psdhp(i,2,3)
      d2psdhp(i,3,1)=-4.0d0*psp(i)*(hp(1)*psp(i)+hp(2))/3.0d0/dummy
      d2psdhp(i,1,3)=d2psdhp(i,3,1)
c
      return
      end
c
c

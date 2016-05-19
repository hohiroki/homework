
!====================================================================================
      subroutine time_step
      !===================
      
      use main  
      
!
!     compute spectral radii
!     ----------------------
!
      
      do j=jb,jm-1
      do i=ib,im-1
          
           aco2=abs(ga*p(i,j)/ro(i,j))
           aco=sqrt(aco2) 
          
!       i direction
          arx=0.5*(aix(i,j)+aix(i+1,j))
          ary=0.5*(aiy(i,j)+aiy(i+1,j))
          vibc=0.5*(vib(i,j)+vib(i+1,j))
          sii=arx*arx+ary*ary
          armx=sqrt(sii)
          velpro=vx(i,j)*arx+vy(i,j)*ary+vibc
          acoa=armx*aco
          eox1=abs(velpro)+acoa

!       j direction
          arx=0.5*(ajx(i,j)+ajx(i,j+1))
          ary=0.5*(ajy(i,j)+ajy(i,j+1))
          vjbc=0.5*(vjb(i,j)+vjb(i,j+1))
          sjj=arx*arx+ary*ary
          army=sqrt(sjj)
          velpro=vx(i,j)*arx+vy(i,j)*ary+vjbc
          acoa=army*aco
          eoy1=abs(velpro)+acoa
          
          ste=cfl*vol(i,j)/(eox1+eoy1)
          step(i,j)=ste
!
      end do
      end do
!
      if(icgl.eq.1) then
!
      ttg=1000.
      
     
      do j=jb,jm-1
       do i=ib,im-1
      if(ttg.gt.step(i,j)) ttg=step(i,j)
!
      end do
      end do
!
      
      do j=jb,jm-1
      do i=ib,im-1
!
      step(i,j)=ttg
!
      end do
      end do
!
      end if
      
      return
      end
!======================================================================================
      
     subroutine interpolation_to_interface
     !===================================
      use main
      
      ! interpolate the primitive variables to cell interfaces
      
      ! iface
     
      do j=jb,jm-1
      do i=ib,im
          alsinv=1./(alagmx (i,j)+alagmx (i-1,j))
          uil (i,j)=(alagmx (i,j)*vx (i-1,j)+alagmx (i-1,j)*vx (i,j))&
     &                  *alsinv
          vil (i,j)=(alagmx (i,j)*vy (i-1,j)+alagmx (i-1,j)*vy (i,j))&
     &                  *alsinv
          
          pil (i,j)=(alagmx (i,j)*p  (i-1,j)+alagmx (i-1,j)*p  (i,j))&
     &                  *alsinv
          ril (i,j)=(alagmx (i,j)*ro (i-1,j)+alagmx (i-1,j)*ro (i,j))&
     &                  *alsinv
     
      end do
      end do
!
!     j face
!
     
      do j=jb,jm
      do i=ib,im-1
          alsinv=1./(alagmy (i,j)+alagmy (i,j-1))
          ujl (i,j)=(alagmy (i,j)*vx (i,j-1)+alagmy (i,j-1)*vx (i,j))&
     &                  *alsinv
          vjl (i,j)=(alagmy (i,j)*vy (i,j-1)+alagmy (i,j-1)*vy (i,j))&
     &                  *alsinv
          
          pjl (i,j)=(alagmy (i,j)*p  (i,j-1)+alagmy (i,j-1)*p  (i,j))&
     &                  *alsinv
          rjl (i,j)=(alagmy (i,j)*ro (i,j-1)+alagmy (i,j-1)*ro (i,j))&
     &                  *alsinv
                 
      end do
      end do
      
      
      return
      end
      
!===========================================================================================      
      subroutine muscl_interpolation 
     !===============================
      use main  
      
      ! MUSCL interpolation at cell interface except boundaries
      
      alimiter(var1,var2)=(var1*(var2*var2+2.*epsm*epsm)+&
     &                     var2*(2.*var1*var1+epsm*epsm))/&
     &            (2.*var1*var1-var1*var2+2.*var2*var2+3.*epsm*epsm) &
     &            *0.5*abs(sign(one,var1)+sign(one,var2)) 
     
   !   alimiter(var1,var2)=2./3.*var1+1./3.*var2
      epsm=small1
      
    ! store cell interface value to facilitate the limiting process   
      
      
      do j=jb,jm-1
      do i=ib,im
      
         uir(i,j)=uil(i,j)
         vir(i,j)=vil(i,j)
         pir(i,j)=pil(i,j)
         rir(i,j)=ril(i,j)
      
      end do
      end do
      
     
      do j=jb,jm
       do i=ib,im-1
      
         ujr(i,j)=ujl(i,j)
         vjr(i,j)=vjl(i,j)
         pjr(i,j)=pjl(i,j)
         rjr(i,j)=rjl(i,j)
         
      end do
      end do
      
      
!
!
!     compute limited cell-interface value
!
      
      do j=jb,jm-1
      do i=ib,im-1
      
      !limiting u
      
      dr=(uil(i+1,j)-  vx(i,j))*2.
      dl=(vx(i,j)    -uir(i,j))*2.
      
      uir(i,j)  =vx(i,j)-0.5*alimiter(dl,dr)
      uil(i+1,j)=vx(i,j)+0.5*alimiter(dr,dl)
      
      dr=(ujl(i,j+1)-  vx(i,j))*2.
      dl=(vx(i,j)    -ujr(i,j))*2.
      
      ujr(i,j)  =vx(i,j)-0.5*alimiter(dl,dr)
      ujl(i,j+1)=vx(i,j)+0.5*alimiter(dr,dl)
      
      
      !limiting v
      
      dr=(vil(i+1,j)-  vy(i,j))*2.
      dl=(vy(i,j)    -vir(i,j))*2.
      
      vir(i,j)  =vy(i,j)-0.5*alimiter(dl,dr)
      vil(i+1,j)=vy(i,j)+0.5*alimiter(dr,dl)
      
      dr=(vjl(i,j+1)-  vy(i,j))*2.
      dl=(vy(i,j)    -vjr(i,j))*2.
      
      vjr(i,j)  =vy(i,j)-0.5*alimiter(dl,dr)
      vjl(i,j+1)=vy(i,j)+0.5*alimiter(dr,dl)
      
      !limiting p
          
      dr=(pil(i+1,j)-  p(i,j))*2.
      dl=(p(i,j)    -pir(i,j))*2.
      
      pir(i,j)  =p(i,j)-0.5*alimiter(dl,dr)
      pil(i+1,j)=p(i,j)+0.5*alimiter(dr,dl)
      
      dr=(pjl(i,j+1)-  p(i,j))*2.
      dl=(p(i,j)    -pjr(i,j))*2.
      
      pjr(i,j)  =p(i,j)-0.5*alimiter(dl,dr)
      pjl(i,j+1)=p(i,j)+0.5*alimiter(dr,dl)
      
      !limiting rou
      
      dr=(ril(i+1,j)-  ro(i,j))*2.
      dl=(ro(i,j)    -rir(i,j))*2.
      
      rir(i,j)  =ro(i,j)-0.5*alimiter(dl,dr)
      ril(i+1,j)=ro(i,j)+0.5*alimiter(dr,dl)
      
      dr=(rjl(i,j+1)-  ro(i,j))*2.
      dl=(ro(i,j)    -rjr(i,j))*2.
      
      rjr(i,j)  =ro(i,j)-0.5*alimiter(dl,dr)
      rjl(i,j+1)=ro(i,j)+0.5*alimiter(dr,dl)
      
      end do
      end do
     
      return
      end
      
!===========================================================================================      
       
      subroutine inviscid_fluxes_roe
!     =========================
      
      use main
      
          
      epsm=0.1
      
      do i=ib-1,im+1
!
      do j=jb-1,jm+1
      
      hro(i,j)=0.
      hre(i,j)=0.
      hrx(i,j)=0.
      hry(i,j)=0.
      
      
      
      end do
      end do
      
          
! evaluate the inviscid fluxes   
!
! -   i direction
!
!
      do j=jb,jm-1
!
      do i=ib,im
      
      
      
      
!    
!     left&right states
!
      
      
      ul=uil(i,j)
      ur=uir(i,j) 
     
      vl=vil(i,j)
      vr=vir(i,j)
     
     
      pl=max(pil(i,j),small)
      pr=max(pir(i,j),small)
     
      rl=max(ril(i,j),small)
      rr=max(rir(i,j),small)
           
      vaml=ul*ul+vl*vl
      vamr=ur*ur+vr*vr
      hl=pl*ga/(rl*ga1)+0.5*vaml
      hr=pr*ga/(rr*ga1)+0.5*vamr 
      acl=sqrt(ga*pl/rl)
      acr=sqrt(ga*pr/rr)
!
!     Roe average
!
      
      rrorl=sqrt(rr/rl)
      rrorlp1=1.+rrorl
      
      rm=sqrt(rr*rl)
      um=(ul+ur*rrorl)/rrorlp1
      vm=(vl+vr*rrorl)/rrorlp1
     
      hm=(hl+hr*rrorl)/rrorlp1
      vm2=um*um+vm*vm
      am2=ga1*abs(hm-0.5*vm2)
      am=sqrt(am2)
      
                  
! 
!     
!     surface area vectors
!
      sav1=aix(i,j)
      sav2=aiy(i,j)
      
      sav=aim(i,j)+small
      sav1n=sav1/sav
      sav2n=sav2/sav
      
!
!     engenvalues
!      
         
           du=(ur-ul)
           dv=(vr-vl)
           
           dunormal=(sav1n*du+sav2n*dv)
           drou=rr-rl
           dp=pr-pl
           dru=rm*du+um*drou
           drv=rm*dv+vm*drou
          
           dre=dp/ga1+0.5*vm2*drou+rm*um*du+rm*vm*dv
           
           sos=am
           sosi=1./sos
      
           qn0=sav1n*um+sav2n*vm
           qn=qn0+vib(i,j)/sav
           aqn=abs(qn)
           amn=qn*sosi
           am0=min(abs(amn),1.)*sign(1.,amn)
           if(abs(amn-1.)<epsm) am0=(-(1.-epsm)**2+2.*(1.+epsm)*amn-amn*amn)/(4.*epsm)
           if(abs(amn+1.)<epsm) am0=((1.-epsm)**2+2.*(1.+epsm)*amn+amn*amn)/(4.*epsm)
           aam0=abs(am0)
           if(aam0<epsm) aam0=(epsm*epsm+aam0*aam0)/(2.*epsm)
           if(aqn<epsm*sos) aqn=(epsm*epsm*sos*sos+aqn*aqn)/(2.*epsm*sos)
           omam0=1.-aam0
           
           sosip=sosi
           sosp=sos 
           
           
!
!     inviscid flux
!      
      coef=aim(i,j)
      dpc1=sosip*omam0*dp
      dpc2=am0*dp
      dqnc1=rm*am0*dunormal
      dqnc2=rm*sosp*omam0*dunormal
      
      
      adr =coef*(              dpc1              +dqnc1+aqn*drou)
      adru=coef*(sav1n*dpc2+um*dpc1+sav1n*dqnc2+um*dqnc1+aqn*dru)
      adrv=coef*(sav2n*dpc2+vm*dpc1+sav2n*dqnc2+vm*dqnc1+aqn*drv)  
      adre=coef*(  qn0*dpc2+hm*dpc1+  qn0*dqnc2+hm*dqnc1+aqn*dre)  
       
      
      ulnormal=(sav1*ul+sav2*vl+vib(i,j))
      urnormal=(sav1*ur+sav2*vr+vib(i,j))
      rulnormal=rl*ulnormal
      rurnormal=rr*urnormal
      
!
!     flux terms
!
      droi=-0.5*(rulnormal+rurnormal-adr)
      dxi=-0.5*(rulnormal*ul+rurnormal*ur&
     &            +sav1*(pl+pr)-adru) 
      dyi=-0.5*(rulnormal*vl+rurnormal*vr&
     &            +sav2*(pl+pr)-adrv)  
      drei=-0.5*(rulnormal*hl+rurnormal*hr&
     &            -vib(i,j)*(pl+pr)-adre)    
    
      
      if((i==83).and.(j==102))then
	  
	  write(2,*)sav1n,sav2n,dpc1,dpc2,dqnc1,dqnc2,adr,qn0,hm,aqn,adre,adru,adrv,droi,drei,dxi,dyi
	  endif
      hro(i,j)=hro(i,j)-droi
      hre(i,j)=hre(i,j)-drei
      hrx(i,j)=hrx(i,j)-dxi
      hry(i,j)=hry(i,j)-dyi
     
      hro(i-1,j)=hro(i-1,j)+droi
      hre(i-1,j)=hre(i-1,j)+drei
      hrx(i-1,j)=hrx(i-1,j)+dxi
      hry(i-1,j)=hry(i-1,j)+dyi
      
      
      
      end do
      end do
      
!
! -   j direction
!
      
!
      do j=jb,jm
      do i=ib,im-1
      
      
!     left&right states
! 
      ul=ujl(i,j)
      ur=ujr(i,j) 
     
      vl=vjl(i,j)
      vr=vjr(i,j)
     
     
      pl=max(pjl(i,j),small)
      pr=max(pjr(i,j),small)
     
      rl=max(rjl(i,j),small)
      rr=max(rjr(i,j),small)
          
      vaml=ul*ul+vl*vl
      vamr=ur*ur+vr*vr
      hl=pl*ga/(rl*ga1)+0.5*vaml
      hr=pr*ga/(rr*ga1)+0.5*vamr 
      acl=sqrt(ga*pl/rl)
      acr=sqrt(ga*pr/rr)
!
!     Roe average
!
      
      rrorl=sqrt(rr/rl)
      rrorlp1=1.+rrorl
      
      rm=sqrt(rr*rl)
      um=(ul+ur*rrorl)/rrorlp1
      vm=(vl+vr*rrorl)/rrorlp1
      hm=(hl+hr*rrorl)/rrorlp1
      vm2=um*um+vm*vm
      am2=ga1*abs(hm-0.5*vm2)
      am=sqrt(am2)
      
                  
! 
!     
!     surface area vectors
!
      sav1=ajx(i,j)
      sav2=ajy(i,j)
      sav=ajm(i,j)+small
      sav1n=sav1/sav
      sav2n=sav2/sav
          
!
!     engenvalues
!      
           du=(ur-ul)
           dv=(vr-vl)
           dunormal=(sav1n*du+sav2n*dv)
           drou=rr-rl
           dp=pr-pl
           dru=rm*du+um*drou
           drv=rm*dv+vm*drou
           dre=dp/ga1+0.5*vm2*drou+rm*um*du+rm*vm*dv 
           
           sos=am
           sosi=1./sos
       
           qn0=sav1n*um+sav2n*vm
           qn=qn0+vjb(i,j)/sav
           aqn=abs(qn)
           amn=qn*sosi
           am0=min(abs(amn),1.)*sign(1.,amn)
           if(abs(amn-1.)<epsm) am0=(-(1.-epsm)**2+2.*(1.+epsm)*amn-amn*amn)/(4.*epsm)
           if(abs(amn+1.)<epsm) am0=((1.-epsm)**2+2.*(1.+epsm)*amn+amn*amn)/(4.*epsm)
           aam0=abs(am0)
           if(aam0<epsm) aam0=(epsm*epsm+aam0*aam0)/(2.*epsm)
           if(aqn<epsm*sos) aqn=(epsm*epsm*sos*sos+aqn*aqn)/(2.*epsm*sos)
           omam0=1.-aam0
           
           sosip=sosi
           sosp=sos 
           
          
           
           
           
           
!
!     inviscid flux
!      
      coef=ajm(i,j)
      dpc1=sosip*omam0*dp
      dpc2=am0*dp
      dqnc1=rm*am0*dunormal
      dqnc2=rm*sosp*omam0*dunormal
      
      adr =coef*(              dpc1              +dqnc1+aqn*drou)
      adru=coef*(sav1n*dpc2+um*dpc1+sav1n*dqnc2+um*dqnc1+aqn*dru)
      adrv=coef*(sav2n*dpc2+vm*dpc1+sav2n*dqnc2+vm*dqnc1+aqn*drv)  
      adre=coef*(  qn0*dpc2+hm*dpc1+  qn0*dqnc2+hm*dqnc1+aqn*dre)    
      
      ulnormal=(sav1*ul+sav2*vl+vjb(i,j))
      urnormal=(sav1*ur+sav2*vr+vjb(i,j))
      rulnormal=rl*ulnormal
      rurnormal=rr*urnormal
      
!
!     flux terms
!
      droj=-0.5*(rulnormal+rurnormal-adr)
      dxj=-0.5*(rulnormal*ul+rurnormal*ur&
     &            +sav1*(pl+pr)-adru) 
      dyj=-0.5*(rulnormal*vl+rurnormal*vr&
     &            +sav2*(pl+pr)-adrv) 
      drej=-0.5*(rulnormal*hl+rurnormal*hr&
     &            -vjb(i,j)*(pl+pr)-adre)
           
      
      hro(i,j)=hro(i,j)-droj
      hre(i,j)=hre(i,j)-drej
      hrx(i,j)=hrx(i,j)-dxj
      hry(i,j)=hry(i,j)-dyj
      
      hro(i,j-1)=hro(i,j-1)+droj
      hre(i,j-1)=hre(i,j-1)+drej
      hrx(i,j-1)=hrx(i,j-1)+dxj
      hry(i,j-1)=hry(i,j-1)+dyj
      
      
      
      end do
      end do     
!  
      return
      end 
     

!========================================================================================
       subroutine inviscid_fluxes_lax
!     =========================
      
      use main
      
          
      epsm=0.15
      
     
!
      do j=jb-1,jm+1
       do i=ib-1,im+1
      
      hro(i,j)=0.
      hre(i,j)=0.
      hrx(i,j)=0.
      hry(i,j)=0.
      
      
      
      end do
      end do
      
          
! evaluate the inviscid fluxes   
!
! -   i direction
!
!
      do j=jb,jm-1
!
      do i=ib,im
      
      
      
      
!    
!     left&right states
!
      
      
      ul=uil(i,j)
      ur=uir(i,j) 
     
      vl=vil(i,j)
      vr=vir(i,j)
     
     
      pl=max(pil(i,j),small)
      pr=max(pir(i,j),small)
     
      rl=max(ril(i,j),small)
      rr=max(rir(i,j),small)
           
      vaml=ul*ul+vl*vl
      vamr=ur*ur+vr*vr
      hl=pl*ga/(rl*ga1)+0.5*vaml
      hr=pr*ga/(rr*ga1)+0.5*vamr 
      el=pl/(rl*ga1)+0.5*vaml
      er=pr/(rr*ga1)+0.5*vamr 
      acl=sqrt(ga*pl/rl)
      acr=sqrt(ga*pr/rr)

                  
! 
!     
!     surface area vectors
!
      sav1=aix(i,j)
      sav2=aiy(i,j)
      
      sav=aim(i,j)+small
      sav1n=sav1/sav
      sav2n=sav2/sav
      

      ulnormal=(sav1*ul+sav2*vl+vib(i,j))
      urnormal=(sav1*ur+sav2*vr+vib(i,j))
      
      

	  coex=max(abs(ulnormal),abs(urnormal))+max(acl,acr)*sav

      rulnormal=rl*ulnormal
      rurnormal=rr*urnormal


      
!
!     flux terms
!
      droi=-0.5*(rulnormal+rurnormal-coex*(rr-rl))
       dxi=-0.5*(rulnormal*ul+rurnormal*ur&
     &            +sav1*(pl+pr)-coex*(rr*ur-rl*ul)) 
       dyi=-0.5*(rulnormal*vl+rurnormal*vr&
     &            +sav2*(pl+pr)-coex*(rr*vr-rl*vl))  
      drei=-0.5*(rulnormal*hl+rurnormal*hr&
     &            -vib(i,j)*(pl+pr)-coex*(rr*er-rl*el))    
    
      
     
      hro(i,j)=hro(i,j)-droi
      hre(i,j)=hre(i,j)-drei
      hrx(i,j)=hrx(i,j)-dxi
      hry(i,j)=hry(i,j)-dyi
     
      hro(i-1,j)=hro(i-1,j)+droi
      hre(i-1,j)=hre(i-1,j)+drei
      hrx(i-1,j)=hrx(i-1,j)+dxi
      hry(i-1,j)=hry(i-1,j)+dyi
      
      
      
      end do
      end do
      
!
! -   j direction
!
     
!
      do j=jb,jm
       do i=ib,im-1
      
      
!     left&right states
! 
      ul=ujl(i,j)
      ur=ujr(i,j) 
     
      vl=vjl(i,j)
      vr=vjr(i,j)
     
     
      pl=max(pjl(i,j),small)
      pr=max(pjr(i,j),small)
     
      rl=max(rjl(i,j),small)
      rr=max(rjr(i,j),small) 
          
      vaml=ul*ul+vl*vl
      vamr=ur*ur+vr*vr
      hl=pl*ga/(rl*ga1)+0.5*vaml
      hr=pr*ga/(rr*ga1)+0.5*vamr 
	    el=pl/(rl*ga1)+0.5*vaml
      er=pr/(rr*ga1)+0.5*vamr 
      acl=sqrt(ga*pl/rl)
      acr=sqrt(ga*pr/rr)

                  
! 
!     
!     surface area vectors
!
      sav1=ajx(i,j)
      sav2=ajy(i,j)
      sav=ajm(i,j)+small
      sav1n=sav1/sav
      sav2n=sav2/sav
          

!
!     inviscid flux
!      
       
      
      ulnormal=(sav1*ul+sav2*vl+vjb(i,j))
      urnormal=(sav1*ur+sav2*vr+vjb(i,j))
	    coey=max(abs(ulnormal),abs(urnormal))+max(acl,acr)*sav
      rulnormal=rl*ulnormal
      rurnormal=rr*urnormal
      
      
!
!     flux terms
!
      droj=-0.5*(rulnormal+rurnormal-coey*(rr-rl))
       dxj=-0.5*(rulnormal*ul+rurnormal*ur&
     &            +sav1*(pl+pr)-coey*(rr*ur-rl*ul)) 
       dyj=-0.5*(rulnormal*vl+rurnormal*vr&
     &            +sav2*(pl+pr)-coey*(rr*vr-rl*vl))  
      drej=-0.5*(rulnormal*hl+rurnormal*hr&
     &            -vjb(i,j)*(pl+pr)-coey*(rr*er-rl*el)) 
           
      
      hro(i,j)=hro(i,j)-droj
      hre(i,j)=hre(i,j)-drej
      hrx(i,j)=hrx(i,j)-dxj
      hry(i,j)=hry(i,j)-dyj
      
      hro(i,j-1)=hro(i,j-1)+droj
      hre(i,j-1)=hre(i,j-1)+drej
      hrx(i,j-1)=hrx(i,j-1)+dxj
      hry(i,j-1)=hry(i,j-1)+dyj
      
      
      
      end do
      end do     
!  
      return
      end

!==============================================================================
      subroutine Runge_Kutta(nrk)
      !=========================
      
      ! compute the right hand side of the Navier-Stokes equations
      
      use main
      
      real:: rkpa(2)
      
      rkpa(1)=1.
      rkpa(2)=0.5
!      
!     update variables
!
      
      
      do  j=jb,jm-1
      do  i=ib,im-1
      
      ro(i,j)=(1.-rkpa(nrk))*rom1(i,j)*vol(i,j)+rkpa(nrk)*ro(i,j)*vol(i,j)+&
     &              rkpa(nrk)*step(i,j)*hro(i,j)
     
      roe(i,j)=(1.-rkpa(nrk))*roem1(i,j)*vol(i,j)+rkpa(nrk)*roe(i,j)*vol(i,j)+&
     &              rkpa(nrk)*step(i,j)*hre(i,j)
     
      rovx(i,j)=(1.-rkpa(nrk))*rovxm1(i,j)*vol(i,j)+rkpa(nrk)*rovx(i,j)*vol(i,j)+&
     &              rkpa(nrk)*step(i,j)*hrx(i,j)
     
       rovy(i,j)=(1.-rkpa(nrk))*rovym1(i,j)*vol(i,j)+rkpa(nrk)*rovy(i,j)*vol(i,j)+&
     &              rkpa(nrk)*step(i,j)*hry(i,j)
     
      end do
      end do

      
      return
      end
      
!===============================================================================================
      subroutine update_variables
!     ===========================

      use main                
      

!
!
!     update variables
!
!
      
      do  j=jb,jm-1
      do  i=ib,im-1
      
!
      volinv=1./vol(i,j)
      ro(i,j)=ro(i,j)*volinv
      roe(i,j)=roe(i,j)*volinv
      rovx(i,j)=rovx(i,j)*volinv
      rovy(i,j)=rovy(i,j)*volinv
      
      roinv=1./ro(i,j)
      
      vx(i,j)=rovx(i,j)*roinv
      vy(i,j)=rovy(i,j)*roinv
      
      
      vsrh=.5*(vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j))
      p(i,j)=ga1*(roe(i,j)-vsrh*ro(i,j))
      ho(i,j)=ga*(roe(i,j)/ro(i,j)-vsrh)+vsrh
!
      end do
      end do
     
!
      return
      end

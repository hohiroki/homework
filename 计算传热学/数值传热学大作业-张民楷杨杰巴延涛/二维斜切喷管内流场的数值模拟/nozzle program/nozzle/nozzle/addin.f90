      subroutine inviscid_fluxes_rtroe
      
      use main
      
      rteps=6.0   
      epsm=0.1
      rtflk=0

      do i=ib-1,im+1
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
         
!     surface area vectors
!
      sav1=aix(i,j)
      sav2=aiy(i,j)
      
      sav=aim(i,j)+small
      sav1n=sav1/sav
      sav2n=sav2/sav
      
!
!     nondissipation flux terms
!
      ulnormal=(sav1*ul+sav2*vl+vib(i,j))
      urnormal=(sav1*ur+sav2*vr+vib(i,j))
      rulnormal=rl*ulnormal
      rurnormal=rr*urnormal

      droi=-0.5*(rulnormal+rurnormal)
      dxi=-0.5*(rulnormal*ul+rurnormal*ur&
     &            +sav1*(pl+pr)) 
      dyi=-0.5*(rulnormal*vl+rurnormal*vr&
     &            +sav2*(pl+pr))  
      drei=-0.5*(rulnormal*hl+rurnormal*hr&
     &            -vib(i,j)*(pl+pr))    

      
         
           du=(ur-ul)
           dv=(vr-vl)
           
           dunormal=(sav1n*du+sav2n*dv)
           drou=rr-rl
           dp=pr-pl
           dru=rm*du+um*drou
           drv=rm*dv+vm*drou
          
           dre=dp/ga1+0.5*vm2*drou+rm*um*du+rm*vm*dv
	   
!     Rotating direction

	   du2=du*du+dv*dv+small
	   adu=sqrt(du2)
	   u2ref=rteps*(vaml+vamr+6*small)
	   if(du2.ge.u2ref)then
!      
       rtflk=rtflk+1
       rtfli=i
	   rtflj=j
	   rtm=2;
	   rtc1=du/adu
	   rts1=dv/adu
	   rt1n=sav1n*rtc1+sav2n*rts1
	   flag=sign(1.,rt1n)
	   rtnc(1,1)=rtc1*flag
	   rtnc(2,1)=rts1*flag
	   rtnc(3,1)=rt1n*flag

	   rtc2=-rts1
           rts2=rtc1     
	   rt2n=sav1n*rtc2+sav2n*rts2
	   flag=sign(1.,rt2n)
	   rtnc(1,2)=rtc2*flag
	   rtnc(2,2)=rts2*flag
	   rtnc(3,2)=rt2n*flag

	   else
	   rtm=1;
	   rtnc(1,1)=sav1n
	   rtnc(2,1)=sav2n   
	   rtnc(3,1)=1.0d0

	   endif
	   
           do rti=1,rtm

	   rtc=rtnc(1,rti)
	   rts=rtnc(2,rti)
	   rtn=rtnc(3,rti)
!                 
!     engenvalues
!
           sos=am
           sosi=1./sos
      
           qn0=rtc*um+rts*vm
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
      
      
      adr =coef*(              dpc1           +dqnc1+aqn*drou)
      adru=coef*(rtc*dpc2+um*dpc1+rtc*dqnc2+um*dqnc1+aqn*dru)
      adrv=coef*(rts*dpc2+vm*dpc1+rts*dqnc2+vm*dqnc1+aqn*drv)  
      adre=coef*(qn0*dpc2+hm*dpc1+qn0*dqnc2+hm*dqnc1+aqn*dre)  
       
      if((i==83).and.(j==102))then
	 
	  write(2,*)rtc,rts,dpc1,dpc2,dqnc1,dqnc2,adr,qn0,hm,aqn,adre,adru,adrv,droi,drei,dxi,dyi
	  
	  endif
 
!    flux term
      droi=droi+0.5*rtn*adr
      drei=drei+0.5*rtn*adre
      dxi=dxi+0.5*rtn*adru    
      dyi=dyi+0.5*rtn*adrv
           
           enddo  

     
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

!     surface area vectors
!
      sav1=ajx(i,j)
      sav2=ajy(i,j)
      sav=ajm(i,j)+small
      sav1n=sav1/sav
      sav2n=sav2/sav
!
!     nondissipation flux terms
!
     
      ulnormal=(sav1*ul+sav2*vl+vjb(i,j))
      urnormal=(sav1*ur+sav2*vr+vjb(i,j))
      rulnormal=rl*ulnormal
      rurnormal=rr*urnormal
      

      droj=-0.5*(rulnormal+rurnormal)
      dxj=-0.5*(rulnormal*ul+rurnormal*ur&
     &            +sav1*(pl+pr)) 
      dyj=-0.5*(rulnormal*vl+rurnormal*vr&
     &            +sav2*(pl+pr)) 
      drej=-0.5*(rulnormal*hl+rurnormal*hr&
     &            -vjb(i,j)*(pl+pr))
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
           
	   du2=du*du+dv*dv+small
	   adu=sqrt(du2)
       
	   u2ref=rteps*(vaml+vamr+6*small)
	   if(du2.ge.u2ref)then
	   rtflk=rtflk+1
	   rtfli=i
	   rtflj=j

	   rtm=2;
	   rtc1=du/adu
	   rts1=dv/adu
	   rt1n=sav1n*rtc1+sav2n*rts1
	   flag=sign(1.,rt1n)
	   rtnc(1,1)=rtc1*flag
	   rtnc(2,1)=rts1*flag
	   rtnc(3,1)=rt1n*flag

	   rtc2=-rts1
           rts2=rtc1     
	   rt2n=sav1n*rtc2+sav2n*rts2
	   flag=sign(1.,rt2n)
	   rtnc(1,2)=rt2n*flag
	   rtnc(2,2)=rtc2*flag
	   rtnc(3,2)=rts2*flag

	   else
	   rtm=1;
	   rtnc(1,1)=sav1n
	   rtnc(2,1)=sav2n   
	   rtnc(3,1)=1.0

	   endif

           do rti=1,rtm

	   rtc=rtnc(1,rti)
	   rts=rtnc(2,rti)
	   rtn=rtnc(3,rti)

           sos=am
           sosi=1./sos
       
           qn0=rtc*um+rts*vm
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
      
      adr =coef*(            dpc1             +dqnc1+aqn*drou)
      adru=coef*(rtc*dpc2+um*dpc1+rtc*dqnc2+um*dqnc1+aqn*dru)
      adrv=coef*(rts*dpc2+vm*dpc1+rts*dqnc2+vm*dqnc1+aqn*drv)  
      adre=coef*(qn0*dpc2+hm*dpc1+qn0*dqnc2+hm*dqnc1+aqn*dre)    
      
 !    flux term
      droj=droj+0.5*rtn*adr
      drej=drej+0.5*rtn*adre
      dxj=dxj+0.5*rtn*adru    
      dyj=dyj+0.5*rtn*adrv
           
           end do         
      
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
      write(*,*)
      return
      end 


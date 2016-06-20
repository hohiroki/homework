     subroutine bc_ghostcell_value

      use main
	  real::machleft,pinleft,Tinleft, rouleft,aleft,uleft,w2left,w3left,w1,aleft2
	  real::ptotal,Ttotal
      
      ! set ghost-cell at boundary for various boundary conditions
      p0=pl0
	  rho0=rol0
      tp0=p0/rho0/rcpcv
	  
      !I-boundary
      
      !inlet
      
      ! 方法：在i=ib-1 这个虚拟网格上给定进口条件 
      

      !暂时用滞止条件  ma=0.4357 u=462.7784/
      i=ib-1
	 machleft=0.4357

     
	  pinleft=p0/(1.0+ga1/2*machleft**2)**rfga
	  Tinleft=tp0/(1.0+ga1/2*machleft**2)
	  aleft=sqrt(ga*rcpcv*Tinleft)
	  uleft=aleft*machleft
	 
	 w2left=uleft-2.0/ga1*aleft

      do j=jb,jm-1
  
  
           w1=vx(i+1,j)+2.0/ga1*sqrt(ga*p(i+1,j)/ro(i+1,j))       
		   vx(i,j)=(w1+w2left)/2.0
           vy(i,j)=0.0
		   Tinleft=tp0-vx(i,j)**2/2.0/cp
		   p(i,j)=p0*(Tinleft/tp0)**rfga
		   ro(i,j)=p(i,j)/rcpcv/Tinleft
 
 
          vmu(i,j)=2.*vmu(i+1,j)-vmu(i+2,j)
         vmul(i,j)=2.*vmul(i+1,j)-vmul(i+2,j)
       alagmx(i,j)=alagmx(i+1,j)
     
      end do

     
  
     !outlet，外推
      i=im
     
      do j=jb,jm-1
     
      p (i,j)=max(2.*p(i-1,j)-p(i-2,j),1.0e-10) 
      vx(i,j)=2.*vx(i-1,j)-vx(i-2,j)
      vy(i,j)=2.*vy(i-1,j)-vy(i-2,j)
      ro(i,j)=max(2.*ro(i-1,j)-ro(i-2,j),1.0e-10)
      vmu(i,j)=2.*vmu(i-1,j)-vmu(i-2,j)
      vmul(i,j)=2.*vmul(i-1,j)-vmul(i-2,j)
      alagmx(i,j)=alagmx(i-1,j)
  
      end do
      
     
      
      
      !j-boundary
      
      ! solid wall boundary
      
      j=jb
      
	  
      do i=ib,im-1
      
      xhalf=0.5*(x(i,j)+x(i+1,j))
      if(xhalf.lt.lx) then  !x>lx,均匀来流条件，否则固壁条件

	  ctmp=(vx(i,j)*ajy(i,j)-vy(i,j)*ajx(i,j))/ajm(i,j)/ajm(i,j)
	  vx(i,j-1)= 2*ajy(i,j)*ctmp-vx(i,j)
	  vy(i,j-1)=-2*ajx(i,j)*ctmp-vy(i,j)


      p (i,j-1)=p (i,j)
      ro(i,j-1)=ro(i,j)
      vmu(i,j-1)=vmu(i,j)
      vmul(i,j-1)=vmul(i,j)
      alagmy(i,j-1)=alagmy(i,j)
      else
      vx(i,j-1)=2*vx(i,j)-vx(i,j+1)
      vy(i,j-1)=2*vy(i,j)-vy(i,j+1)
	  p (i,j-1)=pr0
	  ro(i,j-1)=ror0
      vmu(i,j-1)=2*vmu(i,j)-vmu(i,j+1)
      vmul(i,j-1)=2*vmul(i,j)-vmul(i,j+1)
      alagmy(i,j-1)=alagmy(i,j)
      end if
      end do
      
      j=jm
      
      do i=ib,im-1
      xhalf=0.5*(x(i,j)+x(i+1,j))
      if(xhalf.lt.lm) then  !x>lx,均匀来流条件，否则固壁条件

	  ctmp=(vx(i,j-1)*ajy(i,j)-vy(i,j-1)*ajx(i,j))/ajm(i,j)/ajm(i,j)
	  vx(i,j)= 2*ajy(i,j)*ctmp-vx(i,j-1)
	  vy(i,j)=-2*ajx(i,j)*ctmp-vy(i,j-1)


      p (i,j)=p (i,j-1)
      ro(i,j)=ro(i,j-1)
      vmu(i,j)=vmu(i,j-1)
      vmul(i,j)=vmul(i,j-1)
      alagmy(i,j)=alagmy(i,j-1)
      else
      vx(i,j)=2*vx(i,j-1)-vx(i,j-2)
      vy(i,j)=2*vy(i,j-1)-vy(i,j-2)

	  p (i,j)=pr0
	  ro(i,j)=ror0
      vmu(i,j)=2*vmu(i,j-1)-vmu(i,j-2)
      vmul(i,j)=2*vmul(i,j-1)-vmul(i,j-2)
      alagmy(i,j)=alagmy(i,j-1)
      endif
      end do
       
      return
      end

!================================================================================== 

      subroutine bc_muscl_interpolation
      !================================
      use main
      
      ! set muscl_interpolated value at the cell interfaces on the boundary
      
      
      !I-boundary
      
      !inlet
      
      ! 
      i=ib
      
      do j=jb,jm-1
     
          uil(i,j)=vx(i-1,j)
          vil(i,j)=vy(i-1,j)
          pil(i,j)=p (i-1,j)
          ril(i,j)=ro(i-1,j)
      
      end do
      
      !outlet 
      i=im
      
      do j=jb,jm-1
		  uil(i,j)=vx(i-1,j)
          vil(i,j)=vy(i-1,j)
          pil(i,j)=p (i-1,j)
          ril(i,j)=ro(i-1,j)		

		  uir(i,j)=vx(i-1,j)
          vir(i,j)=vy(i-1,j)
          pir(i,j)=p (i-1,j)
          rir(i,j)=ro(i-1,j)       
      end do
     
      !j-boundary
      
      ! solid wall boundary
      
      j=jb
      
      do i=ib,im-1
      	
      xhalf=0.5*(x(i,j)+x(i+1,j))
      if(xhalf.lt.lx) then
     
      ujr(i,j)=ujl(i,j)
      vjr(i,j)=vjl(i,j)
      pjr(i,j)=pjl(i,j)
      rjr(i,j)=rjl(i,j)
      else

      ujl(i,j)=vx(i,j-1)
      vjl(i,j)=vy(i,j-1)
      pjl(i,j)=p (i,j-1)
      rjl(i,j)=ro(i,j-1)
      end if
 
      end do
      
      j=jm
      
! upper boundary 时间相关的边界条件（先判断激波位置，再决定用激波前或者激波后的参数）
      
      do i=ib,im-1
      xhalf=0.5*(x(i,j)+x(i+1,j))
      if(xhalf.lt.lm) then
      ujr(i,j-1)=ujl(i,j)
      vjr(i,j-1)=vjl(i,j)
      pjr(i,j-1)=pjl(i,j)
      rjr(i,j-1)=rjl(i,j)
      else

      ujr(i,j-1)=vx(i,j)
      vjr(i,j-1)=vy(i,j)
      pjr(i,j-1)=p (i,j)
      rjr(i,j-1)=ro(i,j)
      end if
 
      
      end do

       
      return
      end
     
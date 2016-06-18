
!=========================================================================================
     
      subroutine bc_ghostcell_value
     !===================================
      use main
      
      ! set ghost-cell at boundary for various boundary conditions
      
      
      !I-boundary
      
      !inlet
      
      ! 方法：在i=ib-1 这个虚拟网格上给定进口条件 
      
      !暂时用均匀进口条件
      
      i=ib-1
     
      do j=jb,jm-1
  
           ro(i,j)=max(2.*ro(i+1,j)-ro(i+2,j),1.0e-10)
           vx(i,j)=2.*vx(i+1,j)-vx(i+2,j)
           vy(i,j)=2.*vy(i+1,j)-vy(i+2,j)
            p(i,j)=max(2.*p(i+1,j)-p(i+2,j),1.0e-10)
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
      if(xhalf>=1./6.) then  !x<1/6,均匀来流条件，否则固壁条件
      vx(i,j-1)=vx(i,j)
      vy(i,j-1)=-vy(i,j)
      p (i,j-1)=p (i,j)
      ro(i,j-1)=ro(i,j)
      vmu(i,j-1)=vmu(i,j)
      vmul(i,j-1)=vmul(i,j)
      alagmy(i,j-1)=alagmy(i,j)
      else
      vx(i,j-1)=2*vx(i,j)-vx(i,j+1)
      vy(i,j-1)=2*vy(i,j)-vy(i,j+1)
      p (i,j-1)=2*p (i,j)-p (i,j+1)
      ro(i,j-1)=2*ro(i,j)-ro(i,j+1)
      vmu(i,j-1)=2*vmu(i,j)-vmu(i,j+1)
      vmul(i,j-1)=2*vmul(i,j)-vmul(i,j+1)
      alagmy(i,j-1)=alagmy(i,j)
      end if
      end do
      
      j=jm
      
      do i=ib,im-1
      
      vx(i,j)=2*vx(i,j-1)-vx(i,j-2)
      vy(i,j)=2*vy(i,j-1)-vy(i,j-2)
      p (i,j)=2*p (i,j-1)-p (i,j-2)
      ro(i,j)=2*ro(i,j-1)-ro(i,j-2)
      vmu(i,j)=2*vmu(i,j-1)-vmu(i,j-2)
      vmul(i,j)=2*vmul(i,j-1)-vmul(i,j-2)
      alagmy(i,j)=alagmy(i,j-1)
      
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
     
          a=8.0
          b=7.1447
          c=-4.125
          d=116.5
          uil(i,j)=b
          vil(i,j)=c
          pil(i,j)=d
          ril(i,j)=a
      
      end do
      
      !outlet 
      i=im
      
    !  do j=jb,jm-1
       
    !  end do
     
      !j-boundary
      
      ! solid wall boundary
      
      j=jb
      
      do i=ib,im-1
      	
      xhalf=0.5*(x(i,j)+x(i+1,j))
      if(xhalf>=1./6.) then
     
      ujr(i,j)=ujl(i,j)
      vjr(i,j)=vjl(i,j)
      pjr(i,j)=pjl(i,j)
      rjr(i,j)=rjl(i,j)
      else
          a=8.0
          b=7.1447
          c=-4.125
          d=116.5
      ujl(i,j)=b
      vjl(i,j)=c
      pjl(i,j)=d
      rjl(i,j)=a
      end if
 
      end do
      
      j=jm
      
! upper boundary 时间相关的边界条件（先判断激波位置，再决定用激波前或者激波后的参数）
      
      do i=ib,im-1
      
         xx=0.5*(x(i,j)+x(i+1,j))
         yy=0.5*(y(i,j)+y(i+1,j))
         x0=1/6.+1.73205/3.+ttime*20*1.73205/3.
         flag=1.73205*(xx-x0)-yy+1.
         
         if(flag.gt.0) then

             a=1.4
             b=0.
             c=0.
             d=1.0

           else

             a=8.0
             b=7.1447
             c=-4.125
             d=116.5

           endif
            
      
      ujr(i,j)=b
      vjr(i,j)=c
      pjr(i,j)=d
      rjr(i,j)=a
 
      
      end do

       
      return
      end
     
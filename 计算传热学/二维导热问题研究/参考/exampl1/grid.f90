!=======计算网格参数=========!

subroutine grid()
    use variable
    implicit none
    integer::i,j
	real::hx,hy

	hx=xl/real(im-2)   
	hy=yl/real(jm-2)

!=======计算控制容积界面坐标(xu,yv)，节点坐标(x,y)，节点间距(dx,dy)，控制容积间距(dxv,dyv)=========!

    xu(2)=0
	yv(2)=0    
 	do i=3,im
	   xu(i)=xu(i-1)+hx
	end do
	do j=3,jm
	   yv(j)=yv(j-1)+hy
	end do
	x(1)=xu(2)
	do i=2,imm
	   x(i)=(xu(i+1)+xu(i))/2.0
	   dx(i)=x(i)-x(i-1)
	   dxv(i)=xu(i+1)-xu(i)
	   xuw(i)=xu(i)-x(i-1)
	   xue(i)=x(i)-xu(i)
	end do
	x(im)=xu(im)
	dx(im)=x(im)-x(imm)
	xuw(2)=0
	xuw(im)=xu(im)-x(imm)
	xue(im)=0

	y(1)=yv(2)
	do j=2,jmm
	   y(j)=(yv(j+1)+yv(j))/2.0
	   dy(j)=y(j)-y(j-1)
       dyv(j)=yv(j+1)-yv(j)
	   yvs(j)=yv(j)-y(j-1)
	   yvn(j)=y(j)-yv(j)
	end do
    y(jm)=yv(jm)
	dy(jm)=y(jm)-y(jmm)
    yvs(2)=0
	yvs(jm)=yv(jm)-y(jmm)
	yvn(jm)=0

!=======统一处理直角坐标系与柱坐标系=========!

    if (mode==1) then
	   R(1:jm)=1
	   Ryv(2:jm)=1
	else 
	   Ryv(2)=R0
	   do j=3,jm
	      Ryv(j)=Ryv(j-1)+dyv(j-1)
	   end do 
	   R(1)=R0
	   do j=2,jmm
	      R(j)=(Ryv(j+1)+Ryv(j))/2.0
	   end do
	   R(jm)=Ryv(jm)
	end if
    return
end subroutine
         
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
     

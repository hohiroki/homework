



!================划分网格==================================
	subroutine grid
	use Varieble 
	implicit none 
	real(8)		::	ax,ay
		
	ax=(x2-x1)/M
	ay=(y2-y1)/N


	do i=1,M+1
		x(i)=(i-1)*ax+x1
	enddo
	do j=1,N+1
		y(j)=(j-1)*ay+y1
	enddo
	Call Cal_dxy
	
	end	subroutine
!=========================================================


!=================计算delta_x(y),dx(y)===================================
	subroutine	Cal_dxy
	use varieble
	implicit none 
	
	do i=1,M
		del_x(i)=x(i+1)-x(i)
	enddo 
	do j=1,N
		del_y(j)=y(j+1)-y(j)
	enddo
	
	dx(1)=0.5*del_x(1)
	dx(M+1)=0.5*del_x(M)
	dy(1)=0.5*del_y(1)
	dy(N+1)=0.5*del_y(N)

	do i=2,M
		dx(i)=0.5*(del_x(i)+del_x(i-1))
	enddo 
	do j=2,N
		dy(j)=0.5*(del_y(j)+del_y(j-1))
	enddo
	end subroutine
!======================================================================


!=================全场lambda赋值，计算界面lambda===================================
	subroutine	Cal_lambda
	use varieble
	implicit none 
	
	lambda=400
!	lambda=40
	if(zoom_bcon==2) then
		Call zoom_size
		Call zoom_con2
	endif
	do j=1,N+1
		do i=2,M-1
			lambda_we(i,j)=2*lambda(i+1,j)*lambda(i,j)/(dx(i)*lambda(i+1,j)+dx(i+1)*lambda(i,j))*del_x(i)
		enddo
	enddo
	do j=1,N+1
		lambda_we(1,j)=2*lambda(2,j)*lambda(1,j)/(2*dx(1)*lambda(2,j)+dx(2)*lambda(1,j))*del_x(1)
		lambda_we(M,j)=2*lambda(M+1,j)*lambda(M,j)/(dx(M)*lambda(M+1,j)+2*dx(M+1)*lambda(M,j))*del_x(M)
	enddo

	do j=2,N-1
		do i=1,M+1
			lambda_ns(i,j)=2*lambda(i,j+1)*lambda(i,j)/(dy(j)*lambda(i,j+1)+dy(j+1)*lambda(i,j))*del_y(j)
		enddo
	enddo
	do i=1,M+1
		lambda_ns(i,1)=2*lambda(i,2)*lambda(i,1)/(2*dy(1)*lambda(i,2)+dy(2)*lambda(i,1))*del_y(1)
		lambda_ns(i,N)=2*lambda(i,N+1)*lambda(i,N)/(dy(N)*lambda(i,N+1)+2*dy(N+1)*lambda(i,N))*del_y(N)
	enddo
	


	end subroutine
!======================================================================



!===========================空心区域绝热边界===========
	subroutine zoom_con2
	use varieble
	implicit none
	do i=2,M
		do  j=2,N
			if(zoom(i,j)==1) lambda(i,j)=1e-30
		enddo
	enddo
	end	subroutine
!======================================================


!===========================空心区域等温边界定义===========
	subroutine zoom_con1
	use varieble
	implicit none
	do i=2,M
		do  j=2,N
			if(zoom(i,j)==1) then
					Sc(i,j)=Sc(i,j)+1e30*zoom_TB
					Sp(i,j)=Sp(i,j)-1e30
			endif
		enddo
	enddo
	end	subroutine
!======================================================


!===================空心区域定义===========
	subroutine zoom_size
	use Varieble
	implicit none
	integer(4)	::	i1,j1,i2,j2
	Call grid
	i1=50
	j1=30
	i2=80
	j2=90
!------------------------------椭圆------------------------------
!	do i=2,M
!		do j=2,N
!			if(((x(i)-30)**2/169.+(y(j)-35)**2/100) .le. 1.0) zoom(i,j)=1
!		enddo
!	enddo
!------------------------------------------------------------------


!------------------------------方形------------------------------
	do i=i1,i2
		do j=j1,j2
			zoom(i,j)=1
		enddo
	enddo
!------------------------------------------------------------------

	end	subroutine
!======================================================
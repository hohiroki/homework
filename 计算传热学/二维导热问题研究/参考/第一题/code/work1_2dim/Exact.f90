!

!============Case4精确解，2 dimension==========================================================
	subroutine Exact4_2dim(T)
	use Varieble
	implicit none
	real(8)		::	T(M+1,N+1)
	real(8)		::	x0,y0,a,pi,L1,L2

	L1=x2-x1
	L2=y2-y1

	ET=0
	pi=3.141592653
    do j=1,N+1
		y0=y(j)-y1
		do i=1,M+1
			x0=x(i)-x1
			do k=0,200
				a=(2.0*k+1.0)*pi/2.0/L1
				ET(i,j)=ET(i,j)+((-1)**k)*4/(2.0*k+1.0)/pi&
				&*(exp(a*y0)-exp(-a*y0))/(exp(a*L2)-exp(-a*L2))*cos(a*x0)
			enddo
		enddo
	enddo

	ET=ET*300+300
	

	ET(1,1)=T(1,1)
	ET(1,N+1)=T(1,N+1)
	ET(M+1,1)=T(M+1,1)
	ET(M+1,N+1)=T(M+1,N+1)

	Call output(ET,fexact)
	
	Call output(abs(T-ET)/ET,ferror)
	

	end subroutine

!=========================================================================================


!=================Case3精确解，unsteady，no source=======================================================================
	subroutine Exact3_un_1dim(T)
	use Varieble
	implicit none
	real(8)		::	T(M+1,N+1)
	real(8)		::	a,pi,L

	
	L=y2-y1
	pi=3.141592653
    do j=1,N+1
		do k=0,200
			a=(2.0*k+1.0)*pi/2.0/L
			ET(1,j)=ET(1,j)+4.0/pi/(2*k+1)*((-1)**k)*cos(a*(y(j)-y1))*exp(-a*a*lambda(1,1)/rho/Cp*time)
		enddo
	enddo

	ET(1,:)=ET(1,:)*300+300
	


	Call output_1dim(ET(1,:),fexact)
	Call output_1dim(abs(T(5,:)-ET(1,:))/ET(1,:),ferror)	
	Call output(T,'Cal3_2dim_5.plt')
    print*,sum(abs(T(5,:)-ET(1,:))/ET(1,:))/(N+1)
	end subroutine

!=========================================================================================

!================Case2:Exact_cy_1dim,source=================================================
	subroutine Exact2_cy_1dim(T)
	use Varieble
	implicit none
	real(8)		::	T(M+1,N+1)
	real(8)		::	C1,C2,C3,Tf0,S0

	S0=100
	Tf0=Tf(3)
	C1=S0/4/lambda(1,1)
	C2=y1*y1*C1*2-qB(1)*y1/lambda(1,1)
	C3=Tf0+C1*y2*y2+S0*y2/2/h(3)-C2*(lambda(1,1)/h(3)/y2+log(y2))

	print*,C1,C2,C3
    do j=1,N+1
	
		ET(1,j)=-C1*y(j)*y(j)+C2*log(y(j))+C3
	enddo

	Call output_1dim(ET(1,:),fexact)
	Call output_1dim(abs(T(5,:)-ET(1,:))/ET(1,:),ferror)
	Call output(T,'Cal2_cy_2dim.plt')
	end subroutine
!=============================================================================

!==================Case2:Exact_car_1dim,source=======================================================
	subroutine Exact2_car_1dim(T)
	use Varieble
	implicit none
	real(8)		::	T(M+1,N+1)
	real(8)		::	C1,C2,S0

	S0=Sc(5,5)
	C1=S0/2/lambda(1,1)
	C2=C1*2*y1-qB(1)/lambda(1,1)
	print*,C1,C2
    do j=1,N+1
		ET(1,j)=C1*(y2*y2-y(j)*y(j))+C2*(y(j)-y2)+TB(3)
	enddo

	Call output_1dim(ET(1,:),fexact)
	Call output_1dim(abs(T(5,:)-ET(1,:))/ET(1,:),ferror)
	Call output(T,'Cal2_car_2dim.plt')
	end subroutine
!=============================================================================


!==================Case1精确解,圆柱坐标系======================================================================
	subroutine Exact1_cy_1dim(T)
	use Varieble
	implicit none
	real(8)		::	T(M+1,N+1)
	real(8)		::	C1


	C1=(TB(3)-TB(1))/(log(y2)-log(y1))
	print*,C1
    do j=1,N+1
		ET(1,j)=C1*(log(y(j))-log(y1))+TB(1)
	enddo

	Call output_1dim(ET(1,:),fexact)
	Call output_1dim(abs(T(5,:)-ET(1,:))/ET(1,:),ferror)
	Call output(T,'Cal1_cy_2dim.plt')

	end subroutine
!=============================================================================

!================Case1精确解，直角坐标========================================================================
	subroutine Exact1_car_1dim(T)
	use Varieble
	implicit none
	real(8)		::	T(M+1,N+1)
	real(8)		::	C1


	C1=(TB(3)-TB(1))/(y2-y1)
	print*,C1
    do j=1,N+1
		ET(1,j)=C1*(y(j)-y1)+TB(1)
	enddo

	Call output_1dim(ET(1,:),fexact)
	Call output_1dim(abs(T(5,:)-ET(1,:))/ET(1,:),ferror)
	Call output(T,'Cal1_car_2dim.plt')
	end subroutine
!=============================================================================


!================Case7精确解，直角坐标========================================================================
	subroutine Exact7_car_1dim(T)
	use Varieble
	implicit none
	real(8)		::	T(M+1,N+1)
	real(8)		::	C1

	qB(1)=q0
	C1=qB(1)/4./40
	print*,C1
    do j=1,N+1
		ET(1,j)=C1*(-y(j)+delta)+qB(1)/4/h(1)+tf(1)
	enddo

	Call output_1dim(ET(1,:),'Case7_exact.plt')
	Call output_1dim(abs(T(M/2+1,:)-ET(1,:))/ET(1,:),'Case7_error.plt')
	
	do j=N+1,1,-1
		if(abs(T(M/2+1,j)-ET(M/2+1,j))/ET(M/2+1,j)>0.0001)then
			print*,dy*(j-1)
			exit
		endif
	enddo

	end subroutine
!=============================================================================
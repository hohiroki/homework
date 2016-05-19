


!===============Steady Solver===================================================
	subroutine Steady(T)
	use Varieble
	implicit none
	real(8)			::	T(M+1,N+1)
	real(8)			::	temp(M+1,N+1),a(M+1),b(M+1),c(M+1),d(M+1),residual
	integer(4)		::	lo
	lo=1


	!---------------Initial：迭代初值-----------------
	do i=1,M+1
		do j=1,N+1
			T(i,j)=50
		enddo
	enddo	

	!-------计算差分方程系数--------------------------------
	
	call cal_a


	!------------------Iteration--------------
	do while(lo==1)
		temp=T
	!--------j=1---------------
		a=ap(:,1)
		b=ae(:,1)
		c=aw(:,1)
		d=an(:,1)*T(:,2)+Sc(:,1)
		call TDMA(T(:,1),a,b,c,d)
	!--------j=2-N--------------
		do j=2,N
			a=ap(:,j)
			b=ae(:,j)
			c=aw(:,j)
			d=an(:,j)*T(:,j+1)+as(:,j)*T(:,j-1)+Sc(:,j)
			call TDMA(T(:,j),a,b,c,d)
		enddo
	!--------j=N+1---------------
		a=ap(:,N+1)
		b=ae(:,N+1)
		c=aw(:,N+1)
		d=as(:,N+1)*T(:,N)+Sc(:,N+1)
		call TDMA(T(:,N+1),a,b,c,d)
	!---------------------------
		residual=sum(abs(Temp-T)/T)
		if(residual<=epsilon) lo=0
	enddo
	if(dimen==1)then
		Call output_1dim(T(5,:),fname)
	else
		call output(T,fname)
	endif

	end subroutine
!============================================================================================






!===========Unsteady Solver(explicit)==============================================================
	

	subroutine Unsteady_explic(T)
	use Varieble
	implicit none
	real(8)			::	T(M+1,N+1)
	real(8)			::	T0(M+1,N+1),ep(M+1,N+1)
	Nt=time/dt+0.5


	T=TT0			!初始化
	Call cal_a		!计算差分方程系数
	T0=0
	ep=0

	!--------------时间推进-----------------------------
	do k=1,Nt
		T0=T


	!--------直接推进-------------------------------
		do j=2,N
			do i=2,M
				T(i,j)=(Sc(i,j)+(2*ap0(i,j)-ap(i,j))*T0(i,j)+ae(i,j)*T0(i+1,j)+aw(i,j)*T0(i-1,j)+an(i,j)*T0(i,j+1)+as(i,j)*T0(i,j-1))/ap0(i,j)
			enddo
		enddo	
	
	
	!-----------边界处理-----------------------------
		if(bcon(1)==1)then
			do i=1,M+1
				T(i,1)=TB(1)
			enddo
		else
			do i=2,M
				T(i,1)=(Sc(i,1)+(2*ap0(i,1)-ap(i,1))*T0(i,1)+ae(i,1)*T0(i+1,1)+aw(i,1)*T0(i-1,1)+an(i,1)*T0(i,2))/ap0(i,1)
			enddo
		endif

		if(bcon(2)==1)then
			do j=1,N+1
				T(M+1,j)=TB(2)
			enddo
		else
			do j=2,N
				T(M+1,j)=(Sc(M+1,j)+(2*ap0(M+1,j)-ap(M+1,j))*T0(M+1,j)+aw(M+1,j)*T0(M,j)+an(M+1,j)*T0(M+1,j+1)+as(M+1,j)*T0(M+1,j-1))/ap0(M+1,j)
			enddo
		endif

		if(bcon(3)==1)then
			do i=1,M+1
				T(i,N+1)=TB(3)
			enddo
		else
			do i=2,M
				T(i,N+1)=(Sc(i,N+1)+(2*ap0(i,N+1)-ap(i,N+1))*T0(i,N+1)+ae(i,N+1)*T0(i+1,N+1)+aw(i,N+1)*T0(i-1,N+1)+as(i,N+1)*T0(i,N))/ap0(i,N+1)
			enddo		
		endif

		if(bcon(4)==1)then
			do j=1,N+1
				T(1,j)=TB(4)
			enddo
		else
			do j=2,N
				T(1,j)=(Sc(1,j)+(2*ap0(1,j)-ap(1,j))*T0(1,j)+ae(1,j)*T0(2,j)+an(1,j)*T0(1,j+1)+as(1,j)*T0(1,j-1))/ap0(1,j)
			enddo
		endif
	!-----------------角点处理---------------------------------------------------------
		if((bcon(4) .ne.1) .and. (bcon(1) .ne.1)) then
			T(1,1)=(Sc(1,1)+(2*ap0(1,1)-ap(1,1))*T0(1,1)+ae(1,1)*T0(2,1)+an(1,1)*T0(1,2))/ap0(1,1)
		else if((bcon(1) .ne.1) .and.(bcon(2) .ne.1)) then
			T(M+1,1)=(Sc(M+1,1)+(2*ap0(M+1,1)-ap(M+1,1))*T0(M+1,1)+aw(M+1,1)*T(M,1)+an(M+1,1)*T0(M+1,2))/ap0(M+1,1)
		else if((bcon(2) .ne.1) .and. (bcon(3) .ne.1)) then
			T(M+1,N+1)=(Sc(M+1,N+1)+(2*ap0(M+1,N+1)-ap(M+1,N+1))*T0(M+1,N+1)+aw(M+1,N+1)*T(M,N+1)+as(M+1,N+1)*T0(M+1,N))/ap0(M+1,N+1)
		else if((bcon(3) .ne.1) .and. (bcon(4) .ne.1)) then
			T(1,N+1)=(Sc(1,N+1)+(2*ap0(1,N+1)-ap(1,N+1))*T0(1,N+1)+ae(1,N+1)*T(2,N+1)+as(1,N+1)*T0(1,N))/ap0(1,N+1)
		else
		endif

	enddo
    print* ,Nt
	if(dimen==1)then
		Call output_1dim(T(5,:),fname)	!输出一维数据
	else
		call output(T,fname)			!输出二维数据
	endif
	end subroutine
!============================================================================================







!===========Unsteady Solver(implicit)==============================================================
	
	subroutine Unsteady_implic(T)
	use Varieble
	implicit none
	real(8)			::	T(M+1,N+1)
	real(8)			::	T0(M+1,N+1)
	real(8)			::	temp(M+1,N+1),a(M+1),b(M+1),c(M+1),d(M+1),residual
	integer(4)		::	lo

	lo=1
	Nt=time/dt+0.5
	T=TT0			!初始化
	Call cal_a		!计算系数

	do k=1,Nt
		T0=T
	!--------隐式迭代---------------------------------
	  do while(lo==1)
		temp=T
	!--------j=1---------------
		a=ap(:,1)
		b=ae(:,1)
		c=aw(:,1)
		d=an(:,1)*T(:,2)+Sc(:,1)+ap0(:,1)*T0(:,1)
		call TDMA(T(:,1),a,b,c,d)
	!--------j=2-N--------------
		do j=2,N
			a=ap(:,j)
			b=ae(:,j)
			c=aw(:,j)
			d=an(:,j)*T(:,j+1)+as(:,j)*T(:,j-1)+Sc(:,j)+ap0(:,j)*T0(:,j)
			call TDMA(T(:,j),a,b,c,d)
		enddo
	!--------j=N+1---------------
		a=ap(:,N+1)
		b=ae(:,N+1)
		c=aw(:,N+1)
		d=as(:,N+1)*T(:,N)+Sc(:,N+1)+ap0(:,N+1)*T0(:,N+1)
		call TDMA(T(:,N+1),a,b,c,d)
	!---------------------------
		residual=sum(abs(Temp-T)/T)
		if(residual<=epsilon) lo=0
	  enddo
	  lo=1
	enddo
	print *, Nt

	if(dimen==1)then
		Call output_1dim(T(5,:),fname)  !输出一维数据
	else
		call output(T,fname)			!输出二维数据
	endif


	end subroutine
!============================================================================================


!===========Unsteady Solver(semi-implicit):加权隐式==============================================================
	
	subroutine Unsteady_semi_im(T)
	use Varieble
	implicit none
	real(8)			::	T(M+1,N+1)
	real(8)			::	T0(M+1,N+1)
	real(8)			::	temp(M+1,N+1),a(M+1),b(M+1),c(M+1),d(M+1),residual,ftheta
	real(8)			::	t_an(M+1,N+1),t_as(M+1,N+1),t_ae(M+1,N+1),t_aw(M+1,N+1),t_ap(M+1,N+1)
	integer(4)		::	lo

	lo=1
	Nt=time/dt+0.5
	T=TT0			!初始化

!---------------计算系数---------------------------
	Call cal_a	
	ftheta=1-theta
	t_an=an
	t_as=as
	t_aw=aw
	t_ae=ae
	ae=theta*t_ae
	aw=theta*t_aw
	an=theta*t_an
	as=theta*t_as
	t_ap=ap-ap0
	ap=t_ap*theta+ap0

	!-------边界处理------------------------------
	if(bcon(1)==1)	then
		ap(:,1)=1e30
		Sc(:,1)=1e30*TB(1)
		t_ap(:,1)=0
	else
	endif

	if(bcon(2)==1)	then
		ap(M+1,:)=1e30
		Sc(M+1,:)=1e30*TB(2)
		t_ap(M+1,:)=0
	else
	endif

	if(bcon(3)==1)	then
		ap(:,N+1)=1e30
		Sc(:,N+1)=1e30*TB(3)
		t_ap(:,N+1)=0
	else
	endif

	if(bcon(4)==1)	then
		ap(1,:)=1e30
		Sc(1,:)=1e30*TB(4)
		t_ap(1,:)=0
	else
	endif
	!------------------------------------------

	do k=1,Nt
		T0=T
	!----------隐式迭代------------------------------------------
	  do while(lo==1)
		temp=T
	!--------j=1---------------
		a=ap(:,1)
		b=ae(:,1)
		c=aw(:,1)
				d(1)=an(1,1)*T(1,2)+Sc(1,1)+(ap0(1,1)-ftheta*t_ap(1,1))*T0(1,1)+ftheta*&
			&	(t_an(1,1)*T0(1,2)+t_ae(1,1)*T0(2,1))
				
				d(M+1)=an(M+1,1)*T(M+1,2)+Sc(M+1,1)+(ap0(M+1,1)-ftheta*t_ap(M+1,1))*T0(M+1,1)+ftheta*&
			&	(t_an(M+1,1)*T0(M+1,2)+t_aw(M+1,1)*T0(M,1))
			
			do i=2,M
				d(i)=an(i,1)*T(i,2)+Sc(i,1)+(ap0(i,1)-ftheta*t_ap(i,1))*T0(i,1)+ftheta*&
			&	(t_an(i,1)*T0(i,2)+t_ae(i,1)*T0(i+1,1)+t_aw(i,1)*T0(i-1,1))
			enddo

		call TDMA(T(:,1),a,b,c,d)
	!--------j=2-N--------------
		do j=2,N
			a=ap(:,j)
			b=ae(:,j)
			c=aw(:,j)

				d(1)=an(1,j)*T(1,j+1)+as(1,j)*T(1,j-1)+Sc(1,j)+(ap0(1,j)-ftheta*t_ap(1,j))*T0(1,j)+ftheta*&
			&	(t_an(1,j)*T0(1,j+1)+t_as(1,j)*T0(1,j-1)+t_ae(1,j)*T0(2,j))
				
				d(M+1)=an(M+1,j)*T(M+1,j+1)+as(M+1,j)*T(M+1,j-1)+Sc(M+1,j)+(ap0(M+1,j)-ftheta*t_ap(M+1,j))*T0(M+1,j)+ftheta* &
			&	(t_an(M+1,j)*T0(M+1,j+1)+t_as(M+1,j)*T0(M+1,j-1)+t_aw(M+1,j)*T0(M,j))
			
			do i=2,M
				d(i)=an(i,j)*T(i,j+1)+as(i,j)*T(i,j-1)+Sc(i,j)+(ap0(i,j)-ftheta*t_ap(i,j))*T0(i,j)+ftheta*&
			&	(t_an(i,j)*T0(i,j+1)+t_as(i,j)*T0(i,j-1)+t_ae(i,j)*T0(i+1,j)+t_aw(i,j)*T0(i-1,j))
			enddo
			call TDMA(T(:,j),a,b,c,d)
		enddo
	!--------j=N+1---------------
		a=ap(:,N+1)
		b=ae(:,N+1)
		c=aw(:,N+1)

				d(1)=as(1,N+1)*T(1,N)+Sc(1,N+1)+(ap0(1,N+1)-ftheta*t_ap(1,N+1))*T0(1,N+1)+ftheta*&
			&	(t_as(1,N+1)*T0(1,N)+t_ae(1,N+1)*T0(2,N+1))
				
				d(M+1)=as(M+1,N+1)*T(M+1,N)+Sc(M+1,N+1)+(ap0(M+1,N+1)-ftheta*t_ap(M+1,N+1))*T0(M+1,N+1)+ftheta*&
			&	(t_as(M+1,N+1)*T0(M+1,N)+t_aw(M+1,N+1)*T0(M,N+1))
			
			do i=2,M
				d(i)=as(i,N+1)*T(i,N)+Sc(i,N+1)+(ap0(i,N+1)-ftheta*t_ap(i,N+1))*T0(i,N+1)+ftheta*&
			&	(t_as(i,N+1)*T0(i,N)+t_ae(i,N+1)*T0(i+1,N+1)+t_aw(i,N+1)*T0(i-1,N+1))
			enddo

		call TDMA(T(:,N+1),a,b,c,d)
	!---------------------------
		residual=sum(abs(Temp-T)/T)
		if(residual<=epsilon) lo=0
	  enddo
	  lo=1
	enddo
	print *, Nt

	if(dimen==1)then
		Call output_1dim(T(5,:),fname)	!输出一维数据
	else
		call output(T,fname)			!输出二维数据
	endif


	end subroutine
!============================================================================================





!===========Unsteady Solver(explicit,animate)==============================================================
	subroutine Unsteady_explic_ani(T)
	use Varieble
	implicit none
	real(8)			::	T(M+1,N+1)
	real(8)			::	T0(M+1,N+1),ep(M+1,N+1)

	Nt=time/dt+0.5

	open(15,file='unsteady_ani.plt')
	40 format('title="Varieble1"')
	10 format('Variable= "x","y","varieble1"')
	20 format('Zone T="variable1", i=',i5,',j=',i5)

	T=TT0

	write(15,40)  
	write(15,10) 
	write(15,20) M+1,N+1
	do j=1,N+1
		do i=1,M+1
			 write(15,"(f8.4,3x,f8.4,3x,f18.4)") dx*(i-1),dy*(j-1),T(i,j)
		enddo
	enddo
	write(15,"()")

	Call cal_a
	T0=0
	ep=0
	do k=1,Nt
		T0=T
	!---------------------------------------
		do j=2,N
			do i=2,M
				T(i,j)=(Sc(i,j)+(2*ap0(i,j)-ap(i,j))*T0(i,j)+ae(i,j)*T0(i+1,j)+aw(i,j)*T0(i-1,j)+an(i,j)*T0(i,j+1)+as(i,j)*T0(i,j-1))/ap0(i,j)
			enddo
		enddo	
	!------------边界处理--------------------------	
		
		if(bcon(1)==1)then
			do i=1,M+1
				T(i,1)=TB(1)
			enddo
		else
			do i=2,M
				T(i,1)=(Sc(i,1)+(2*ap0(i,1)-ap(i,1))*T0(i,1)+ae(i,1)*T0(i+1,1)+aw(i,1)*T0(i-1,1)+an(i,1)*T0(i,2))/ap0(i,1)
			enddo
		endif

		if(bcon(2)==1)then
			do j=1,N+1
				T(M+1,j)=TB(2)
			enddo
		else
			do j=2,N
				T(M+1,j)=(Sc(M+1,j)+(2*ap0(M+1,j)-ap(M+1,j))*T0(M+1,j)+aw(M+1,j)*T0(M,j)+an(M+1,j)*T0(M+1,j+1)+as(M+1,j)*T0(M+1,j-1))/ap0(M+1,j)
			enddo
		endif

		if(bcon(3)==1)then
			do i=1,M+1
				T(i,N+1)=TB(3)
			enddo
		else
			do i=2,M
				T(i,N+1)=(Sc(i,N+1)+(2*ap0(i,N+1)-ap(i,N+1))*T0(i,N+1)+ae(i,N+1)*T0(i+1,N+1)+aw(i,N+1)*T0(i-1,N+1)+as(i,N+1)*T0(i,N))/ap0(i,N+1)
			enddo		
		endif

		if(bcon(4)==1)then
			do j=1,N+1
				T(1,j)=TB(4)
			enddo
		else
			do j=2,N
				T(1,j)=(Sc(1,j)+(2*ap0(1,j)-ap(1,j))*T0(1,j)+ae(1,j)*T0(2,j)+an(1,j)*T0(1,j+1)+as(1,j)*T0(1,j-1))/ap0(1,j)
			enddo
		endif
	!----------角点处理----------------------------------
		if((bcon(4) .ne.1) .and. (bcon(1) .ne.1)) then
			T(1,1)=(Sc(1,1)+(2*ap0(1,1)-ap(1,1))*T0(1,1)+ae(1,1)*T0(2,1)+an(1,1)*T0(1,2))/ap0(1,1)
		else if((bcon(1) .ne.1) .and.(bcon(2) .ne.1)) then
			T(M+1,1)=(Sc(M+1,1)+(2*ap0(M+1,1)-ap(M+1,1))*T0(M+1,1)+aw(M+1,1)*T(M,1)+an(M+1,1)*T0(M+1,2))/ap0(M+1,1)
		else if((bcon(2) .ne.1) .and. (bcon(3) .ne.1)) then
			T(M+1,N+1)=(Sc(M+1,N+1)+(2*ap0(M+1,N+1)-ap(M+1,N+1))*T0(M+1,N+1)+aw(M+1,N+1)*T(M,N+1)+as(M+1,N+1)*T0(M+1,N))/ap0(M+1,N+1)
		else if((bcon(3) .ne.1) .and. (bcon(4) .ne.1)) then
			T(1,N+1)=(Sc(1,N+1)+(2*ap0(1,N+1)-ap(1,N+1))*T0(1,N+1)+ae(1,N+1)*T(2,N+1)+as(1,N+1)*T0(1,N))/ap0(1,N+1)
		else
		endif
	!----------20步输出一组数据------------------------------------------
	if(mod(k,20)==0) then
		write(15,40)  
		write(15,10) 
		write(15,20) M+1,N+1
		do j=1,N+1
			do i=1,M+1
				 write(15,"(f8.4,3x,f8.4,3x,f18.4)") dx*(i-1),dy*(j-1),T(i,j)
			enddo
		enddo
		write(15,"()")
	else
	endif
	enddo
	!-------------------------------------------------
	close(15)	
    print* ,Nt
	Call output(T,fname)	!输出二维数据
	end subroutine
!============================================================================================
	





!===========Unsteady Solver(implicit,Case7)==============================================================
	
	subroutine Unsteady_implic_7(T)
	use Varieble
	implicit none
	real(8)			::	T(M+1,N+1)
	real(8)			::	T0(M+1,N+1)
	real(8)			::	temp(M+1,N+1),a(M+1),b(M+1),c(M+1),d(M+1),residual,tstep
	integer(4)		::	lo,count

	lo=1
	Nt=time/dt+0.5
	T=TT0			!初始化
	Call cal_a		!计算系数

	open(15,file='result_1dim.plt')
	open(5,file='boundy1_-4.plt')
	open(6,file='boundy2_-4.plt')

	write(15,10)  
	write(5,10) 
	write(6,10)
	tstep=0
	count=0
	print *, Nt
!	pause
	do k=1,Nt
		T0=T
!		print*,k
	!--------隐式迭代---------------------------------


	  do while(lo==1)
		temp=T

		tstep=tstep+dt
		if(tstep<=0.015)then
			qB(1)=3e5
		elseif(tstep<0.06)then
			qB(1)=0
			if(tstep-dt<=0.015)then
				Sc=0
				Sp=0
				Call cal_a
			endif
		elseif(tstep>=0.06)then
			qB(1)=3e5
			tstep=0
				Sc=0
				Sp=0
				Call cal_a
		else
			print*,'Erorr!'
		endif


	!--------j=1---------------
		a=ap(:,1)
		b=ae(:,1)
		c=aw(:,1)
		d=an(:,1)*T(:,2)+Sc(:,1)+ap0(:,1)*T0(:,1)
		call TDMA(T(:,1),a,b,c,d)
	!--------j=2-N--------------
		do j=2,N
			a=ap(:,j)
			b=ae(:,j)
			c=aw(:,j)
			d=an(:,j)*T(:,j+1)+as(:,j)*T(:,j-1)+Sc(:,j)+ap0(:,j)*T0(:,j)
			call TDMA(T(:,j),a,b,c,d)
		enddo
	!--------j=N+1---------------
		a=ap(:,N+1)
		b=ae(:,N+1)
		c=aw(:,N+1)
		d=as(:,N+1)*T(:,N)+Sc(:,N+1)+ap0(:,N+1)*T0(:,N+1)
		call TDMA(T(:,N+1),a,b,c,d)
	!---------------------------
		residual=sum(abs(Temp-T)/T)
		if(residual<=epsilon) lo=0
	  enddo
	  lo=1
	  	if(mod(k,1000)==0)then
			write(5,"(f12.7,3x,f18.4)"), dt*k,T(M/2+1,1)
			write(6,"(f12.7,3x,f18.4)"), dt*k,T(M/2+1,N+1)
!			count=count+1
!			write(15,20),count,N+1
!			do j=1,N+1
!				write(15,"(f12.7,3x,f18.4)") (j-1)*dy,T(M/2+1,j)
!			enddo
!			write(15,"()")
		endif
	enddo



10 format('title="temperature"')
20 format('Zone T="n=', i5,'",j=',i5)
30 format('Zone T="boundry"',',j=',i5)

		close(15)
		close(5)
		close(6)
		Call output_1dim(T(M/2+1,:),fname)  !输出一维数据
		call output(T,'Case7_2dim.plt')			!输出二维数据



	end subroutine
!============================================================================================
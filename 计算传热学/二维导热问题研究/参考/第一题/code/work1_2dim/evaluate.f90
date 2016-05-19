
!====================================================================
	subroutine evaluate
	use Varieble
	implicit none

!------------Varieble----------------------

	rho=5 !2719    
	Cp=5 !871		

	M=60
	N=50
	time=5
	dt=0.01
!-----------------------
	x1=0.0
	x2=60.0
	y1=10.0
	y2=60.0



	epsilon=1e-8


	exact_symbol=0	
	theta=2				!默认为稳态问题
	dimen=2				!默认为二维问题
	zoom_bcon=-1		!默认为实心
	animate=0			!默认为不输出多组结果

!	Call Case1_car_1dim
!	Call Case1_cy_1dim
!	Call Case2_car_1dim
!	Call Case2_cy_1dim
!	Call Case3_un_1dim
!	Call Case4_car_2dim
!	Call Case5_un_2dim
!	Call Case5_steady_2dim
!	Call Case6_un_2dim
!	Call Case7_un_1dim
	Call Case8_zoom

!分配内存
	allocate(ET(M+1,N+1))
	allocate(rn(N+1),rs(N+1))
	allocate(ae(M+1,N+1),aw(M+1,N+1),an(M+1,N+1),as(M+1,N+1),ap(M+1,N+1),ap0(M+1,N+1))
	allocate(Sc(M+1,N+1),Sp(M+1,N+1),lambda_we(M,N+1),lambda_ns(M+1,N),lambda(M+1,N+1))
	allocate(del_x(M),del_y(N),dx(M+1),dy(N+1),x(M+1),y(N+1))
	allocate(zoom(M+1,N+1))
	zoom=0
!-------------------------------------------
	if(zoom_bcon==1) then
		Call zoom_size		!计算空心区域大小
		Call zoom_con1		!空心区域边界为等温边界
	endif

!	theta=1.0			!隐式格式
!	theta=0.0			!显式格式
!	theta=0.5			!c-n格式
!	theta=2				!稳态问题

	



!	Call Test_un


!	Call Symm
!	Call Unsteady(1)
!	Call Sample2(1)
!	call test(1)
	Call cal_r			!计算r

!----------Source-----------------------------

	do i=2,M
		Sp(i,1)=Sp(i,1)*dx(i)*dy(1)*(rn(1)+rs(1))/2
		Sc(i,1)=Sc(i,1)*dx(i)*dy(1)*(rn(1)+rs(1))/2
		Sp(i,N+1)=Sp(i,N+1)*dx(i)*dy(1)*(rn(N+1)+rs(N+1))/2
		Sc(i,N+1)=Sc(i,N+1)*dx(i)*dy(1)*(rn(N+1)+rs(N+1))/2
		do j=2,N
			Sp(i,j)=Sp(i,j)*dx(i)*dy(j)*(rn(j)+rs(j))/2
			Sc(i,j)=Sc(i,j)*dx(i)*dy(j)*(rn(j)+rs(j))/2
		enddo
	enddo
	do j=2,N
		Sp(1,j)=Sp(1,j)*dx(1)*dy(j)*(rn(j)+rs(j))/2
		Sc(1,j)=Sc(1,j)*dx(1)*dy(j)*(rn(j)+rs(j))/2
		Sp(M+1,j)=Sp(M+1,j)*dx(M+1)*dy(j)*(rn(j)+rs(j))/2
		Sc(M+1,j)=Sc(M+1,j)*dx(M+1)*dy(j)*(rn(j)+rs(j))/2
	enddo
		Sp(1,1)=Sp(1,1)*dx(1)*dy(1)*(rn(1)+rs(1))/2
		Sc(1,1)=Sc(1,1)*dx(1)*dy(1)*(rn(1)+rs(1))/2
		Sp(1,N+1)=Sp(1,N+1)*dx(1)*dy(N+1)*(rn(N+1)+rs(N+1))/2
		Sc(1,N+1)=Sc(1,N+1)*dx(1)*dy(N+1)*(rn(N+1)+rs(N+1))/2
		Sp(M+1,1)=Sp(M+1,1)*dx(M+1)*dy(1)*(rn(1)+rs(1))/2
		Sc(M+1,1)=Sc(M+1,1)*dx(M+1)*dy(1)*(rn(1)+rs(1))/2
		Sp(M+1,N+1)=Sp(M+1,N+1)*dx(M+1)*dy(N+1)*(rn(N+1)+rs(N+1))/2
		Sc(M+1,N+1)=Sc(M+1,N+1)*dx(M+1)*dy(N+1)*(rn(N+1)+rs(N+1))/2
end subroutine
!================================================================================================




!==================Case1:Exact_car_1dim,no source=========================
	subroutine Case1_car_1dim
	use Varieble
	implicit none
	coordinate=1
	fname='Cal1_car_1dim.plt'
	fexact='Exact1_car_1dim.plt'
	ferror='Error1_car_1dim.plt'
	dimen=1
	exact_symbol=1
	theta=2				!稳态问题
!--------Boundary conditions----------

	TB(1)=500
	TB(3)=300
	qB=0
	bcon=2
	bcon(3)=1
	bcon(1)=1
!------------source---------------------
	Sp=0
	Sc=0
	end subroutine
!========================================================



!=================Case1:Exact_cy_1dim,no source==========================
	subroutine Case1_Cy_1dim
	use Varieble
	implicit none

	coordinate=2
	fname='Cal1_cy_1dim.plt'
	fexact='Exact1_cy_1dim.plt'
	ferror='Error1_cy_1dim.plt'
	dimen=1
	exact_symbol=1
	theta=2				!稳态问题
!--------Boundary conditions----------

	TB(1)=500
	TB(3)=300
	qB=0
	bcon=2
	bcon(3)=1
	bcon(1)=1

!------------source---------------------
	Sp=0
	Sc=0
	end subroutine
!===========================================================================

!==================Case2:Exact_car_1dim,source=========================
	subroutine Case2_car_1dim
	use Varieble
	implicit none
	coordinate=1

	fname='Cal2_car_1dim.plt'
	fexact='Exact2_car_1dim.plt'
	ferror='Error2_car_1dim.plt'
	dimen=1
	exact_symbol=2
	theta=2				!稳态问题
!--------Boundary conditions----------

	TB(3)=300
	qB=0
	qB(1)=600

	bcon=2
	bcon(3)=1

!------------source---------------------
	Sp=0
	Sc=100
	end subroutine
!========================================================



!=================Case2:Exact_cy_1dim,source==========================
	subroutine Case2_Cy_1dim
	use Varieble
	implicit none

	coordinate=2

	fname='Cal2_cy_1dim.plt'
	fexact='Exact2_cy_1dim.plt'
	ferror='Error2_cy_1dim.plt'
	dimen=1
	exact_symbol=2
	theta=2				!稳态问题
!--------Boundary conditions----------

	h(3)=150
	Tf=300
	qB=0
	qB(1)=400
	bcon=2
	bcon(3)=3

!------------source---------------------
	Sp=0
	Sc=100
	end subroutine
!===========================================================================


!============Case3:Unsteady_1dim,no source,T0=600======================================
	subroutine Case3_Un_1dim
	use Varieble
	implicit none
	
	TT0=600
	coordinate=1
	if(theta==0.)then
		fname='Exp3_1dim.plt'
		ferror='Ex3_Error.plt'
	elseif(theta==1.)then
		fname='Imp3_1dim.plt'
		ferror='Im3_Error.plt'
	elseif(theta==0.5)then
		fname='CN3_1dim.plt'
		ferror='CN3_Error.plt'
	else
		fname='Case3_1dim.plt'
		ferror='Cas3_Error.plt'
	endif
	fexact='Exact3_un.plt'
	exact_symbol=3
	dimen=1
	!	theta=1.0			!隐式格式
	theta=0.0			!显式格式
!	theta=0.5			!c-n格式

!	animate=1
!--------Boundary conditions----------
	TB(3)=300
!	TB(1)=500
	qB=0
	bcon=2
	bcon(3)=1
!	bcon(1)=1

!------------source---------------------
	Sp=0
	Sc=0
	end subroutine
!========================================================


!==============Case4:steady_2dim,no source==========================
	subroutine Case4_car_2dim
	use Varieble
	implicit none
	coordinate=1
	fname='Case4_2dim.plt'
	fexact='Exact4_2dim.plt'
	ferror='Error4_2dim.plt'
	exact_symbol=4
	theta=2				!稳态问题
	rho=5 !2719    
	Cp=5 !871		

	M=60
	N=50
!-----------------------
	x1=0.0
	x2=60.0
	y1=10.0
	y2=60.0

!--------Boundary conditions----------
	TB(1)=300
	TB(2)=300
	TB(3)=600
	qB(4)=0
	bcon=1
	bcon(4)=2
!------------source---------------------
	Sp=0
	Sc=0

	end subroutine
!========================================================

!===================Case5==========================
	subroutine Case5_un_2dim
	use Varieble
	implicit none


	TT0=300

	coordinate=1
	if(theta==0.)then
		fname='Exp5_2dim.plt'
		ferror='Ex5_Error.plt'
	elseif(theta==1.)then
		fname='Imp5_2dim.plt'
		ferror='Im5_Error.plt'
	elseif(theta==0.5)then
		fname='CN5_2dim.plt'
		ferror='CN5_Error.plt'
	else
		fname='Case5_2dim.plt'
		ferror='Cas5_Error.plt'
	endif
	fexact='Exact5_2dim.plt'

	!	theta=1.0			!隐式格式
	theta=0.0			!显式格式
!	theta=0.5			!c-n格式

!	animate=1
!--------Boundary conditions----------
	TB(1)=300
	qB(3)=2000
	TB(4)=500
	Tf=400
	h(2)=150
	bcon=1
	bcon(2)=3
	bcon(3)=2
!------------source---------------------
	Sp=0
	Sc=100
	end subroutine
!========================================================
!==============Case5:steady_2dim,source==========================
	subroutine Case5_steady_2dim
	use Varieble
	implicit none
	coordinate=1
	fname='Case5_ste_2dim.plt'
	fexact='Exact5_2dim.plt'
	ferror='Error5_2dim.plt'
!--------Boundary conditions----------
	TB(1)=300
	qB(3)=2000
	TB(4)=500
	Tf=400
	h(2)=150
	bcon=1
	bcon(2)=3
	bcon(3)=2
!------------source---------------------
	Sp=0
	Sc=100

	end subroutine
!========================================================


!===================Case6==========================
	subroutine Case6_un_2dim
	use Varieble
	implicit none
	
	TT0=300
	coordinate=1
	if(theta==0.)then
		fname='Exp6_2dim.plt'
		ferror='Ex6_Error.plt'
	elseif(theta==1.)then
		fname='Imp6_2dim.plt'
		ferror='Im6_Error.plt'
	elseif(theta==0.5)then
		fname='CN6_2dim.plt'
		ferror='CN6_Error.plt'
	else
		fname='Case6_2dim.plt'
		ferror='Cas6_Error.plt'
	endif
	fexact='Exact6_2dim.plt'


	!	theta=1.0			!隐式格式
	theta=0.0			!显式格式
!	theta=0.5			!c-n格式

!	animate=1
!--------Boundary conditions----------
	TB=500

	bcon=1


!------------source---------------------
	Sp=0
	Sc=100
	end subroutine
!========================================================



!============Case7:Unsteady_1dim,周期性边界,T0=500======================================
	subroutine Case7_Un_1dim
	use Varieble
	implicit none
	
	TT0=500
	coordinate=1

	fname='Case7_1dim.plt'
	theta=-1

	!---------Case7--------
	rho=7700
	Cp=460

	M=20
	N=25
	num=5000
	time=0.06*num
	dt=1e-4
	!----------------------

	!-------Case7--------------------
	delta=0.008
	x1=0.0
	x2=delta
	y1=0.0
	y2=delta
	!--------------------------------

!	animate=1
!--------Boundary conditions----------
	h=600
	tf=260
	qB=0
	q0=3e5
	qB(1)=q0
	bcon=2
	bcon(3)=3

!------------source---------------------
	Sp=0
	Sc=0
	end subroutine
!========================================================



!==============Case8：椭圆区域,等温边界(1)和绝热边界(2)==========================
	subroutine Case8_zoom
	use Varieble
	implicit none
	coordinate=2
	fname='Case8_zoom1_1_cy.plt'
	theta=2				!稳态问题

	rho=5 !2719    
	Cp=5 !871		
	x1=0.0
	x2=60.0
	y1=0.0
	y2=50.0
	M=100
	N=100
	zoom_bcon=1			!等温边界(1)和绝热边界(2)
	zoom_TB=400
!--------Boundary conditions----------
	TB(1)=300
	TB(2)=300
	TB(3)=600
	qB(4)=0
	qB(1)=0
	bcon=1
	bcon(1)=2
	bcon(4)=2
!------------source---------------------

	Sp=0
	Sc=0

	end subroutine
!========================================================



!============test======================================
	subroutine Test_un
	use Varieble
	implicit none
	theta=0.2
	TT0=600
	coordinate=1
	if(theta==0.)then
		fname='Expte_1dim.plt'
		ferror='Exte_Error.plt'
	elseif(theta==1.)then
		fname='Impte_1dim.plt'
		ferror='Imte_Error.plt'
	elseif(theta==0.5)then
		fname='CNte_1dim.plt'
		ferror='CNte_Error.plt'
	else
		fname='Casete_1dim.plt'
		ferror='Caste_Error.plt'
	endif
	fexact='Exactte_un.plt'
!	exact_symbol=3
	dimen=1
!	animate=1
!--------Boundary conditions----------
	TB(3)=300
!	TB(1)=500
	qB=0
	bcon=2
	bcon(3)=1
!	bcon(1)=1

!------------source---------------------
	Sp=0
	Sc=0
	end subroutine
!========================================================






































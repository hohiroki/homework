





!================main��������========================================
	program project1_2dim

	use Varieble

	implicit none
	real(8),allocatable	::	T(:,:)
	real(8),allocatable	::	a(:),b(:),c(:),d(:)
	real(8)				::	residual


!-------������ֵ-------------------
	Call evaluate

!----------��������ͼ����ȵ���-------------	
	Call grid
	Call Cal_lambda
!-------------------------------------
!-------���������ڴ�-------------------
	allocate(a(M+1),b(M+1),c(M+1),d(M+1))
	allocate(T(M+1,N+1))
	

!---------�жϲ�����һ�������---------------------------
	if(theta==-1)then
		print*,'case7'
		Call Unsteady_implic_7(T)
		Call Exact7_car_1dim(T)
	elseif(theta==2)then
		Call Steady(T)
	elseif(theta==0)then
		if(animate==0) then
			Call Unsteady_explic(T)
		else
			Call Unsteady_explic_ani(T)
		endif
	elseif(theta==1)then
		Call Unsteady_implic(T)
	elseif(theta>0 .and. theta<1)then
		Call Unsteady_semi_im(T)
	else
	endif
!-------------------------------------------------------



!-----------�ж��Ƿ������ȷ��---------------------------
	if(dimen==1)then
	if(exact_symbol==1 .and. coordinate==1)then
		Call Exact1_car_1dim(T)
	else if(exact_symbol==1 .and. coordinate==2)then
		Call Exact1_cy_1dim(T)
	else if(exact_symbol==2 .and. coordinate==1)then
		Call Exact2_car_1dim(T)
	else if(exact_symbol==2 .and. coordinate==2)then
		Call Exact2_cy_1dim(T)
	else if(exact_symbol==3)then
		Call Exact3_un_1dim(T)
	else
	endif
	else if(exact_symbol==4)then
		Call Exact4_2dim(T)
	else
	endif


	!----------------�ͷ��ڴ�--------------------------------
	deallocate(a,b,c,d)
	deallocate(T,ET)
	deallocate(rn,rs)
	deallocate(ae,aw,an,as,ap,ap0)
	deallocate(Sc,Sp,lambda,lambda_we,lambda_ns)
	deallocate(del_x,del_y,dx,dy,x,y)
	deallocate(zoom)
	end program
!============================================================



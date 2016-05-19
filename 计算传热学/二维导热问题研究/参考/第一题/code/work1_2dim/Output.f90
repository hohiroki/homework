


!===========================二维数据的输出====================
!T：一维数据；fname1：文件名

subroutine output(T,fname1)
use Varieble
implicit none
real(8)			::	T(M+1,N+1)
character(20)	::	fname1

if(zoom_bcon .ne. -1)then
	do i=2,M
		do  j=2,N
			if(zoom(i,j)==1) T(i,j)=0
		enddo
	enddo
endif

open(15,file=fname1)
	write(15,40)  
	40 format('title="Varieble1"')
     

	10 format('Variable= "x","y","varieble1"')
	20 format('Zone T="variable1", i=',i5,',j=',i5)
	write(15,10) 
	write(15,20) M+1,N+1
	do j=1,N+1
		do i=1,M+1
			 write(15,"(f12.8,3x,f12.8,3x,f18.8)") x(i),y(j),T(i,j)
		enddo
	enddo

	close(15)		
end subroutine

!================================================================================





!=========一维数据的输出==================================================
!T：一维数据；fname1：文件名
subroutine output_1dim(T,fname1)
use Varieble
implicit none
real(8)			::	T(N+1)
character(20)	::	fname1

open(15,file=fname1)
	write(15,40)  
	40 format('title="Varieble1"')
     

	do j=1,N+1
		write(15,"(f12.8,3x,f18.8)") y(j),T(j)
	enddo

	close(15)		
end subroutine

!================================================================================



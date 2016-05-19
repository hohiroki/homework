!This is a program for calculationg the heat conduction of a wall.
! Formulations in detail can be found in Book < Numerical Transfer & Heat Conduction> 
! by Patankar.
! In this program,we use uniform grid to compute unsteady changing temperature 
! of a wall.
! We compare two numerical methods : point iteration and block iteration.

program steady_heat_transfer_boundary_conditon_1
    use variable
    implicit none
	integer::i,j,n 
    last=10 
    call initial()
	call grid()
    call start()  
    do while(.not.Converged_all) 
	    time=time+dt
		! 边界赋值或更新源项内容.
        call boundary()
        ! Compute coefficient ap,ae.... for each nodes.
        call compute()
        ! Compute temperature field at each time step.  
        call solve()
        ! Check whether achieve a steady state ( Converged_all )
        call check_all() 
        tlevel=tlevel+1          
    end do 
	if (Converged_all) then
	    call output()
	else 
	write(*,*) "The calculation does not converged please check!"
	end if
end program

subroutine initial()
    use variable
    implicit none
    integer::i,j
    xl=1 
    yl=2
	mode=1
	Cp=1
	rou=0
	last=10
	T0=0
	k0=1
    R0=1 ! 直角坐标R0=1，柱坐标R0为r方向初始半径
    eps=1E-05
    dt=6000.
    Converged_in=.false.
    Converged_all=.false.
    time=0.0
    tlevel=0
    return
end subroutine

!=======节点赋初值========!
subroutine start()
    use variable
    implicit none
    integer::i,j
	do j=1,jm
	    do i=1,im
		   if (i==1.or.i==im.or.j==1.or.j==jm) then
		   	   TN(i,j)=x(i)+y(j)+x(i)*y(j)
               TN_all(i,j)=x(i)+y(j)+x(i)*y(j)
	           T(i,j)=x(i)+y(j)+x(i)*y(j)
			else
			   TN(i,j)=T0
               TN_all(i,j)=T0
		  	   T(i,j)=T0
			end if
		end do
	end do
return
end subroutine


subroutine boundary()
    use variable
    implicit none
    integer::i,j
	do j=1,jm
	   do i=1,im
		  k(i,j)=k0  !第(i,j)控制容积导热系数
	   end do
	end do

!===========================================
!  界面扩散系数节点扩散系数调和平均进行计算 
!===========================================
   !=========x方向界面导热系数==========!
	do j=2,jmm
	   kx(2,j)=k(1,j)
	   do i=3,imm
        kx(i,j)=k(i,j)*k(i-1,j)*dx(i)/(k(i,j)*xuw(i)+k(i-1,j)*xue(i))
	   end do
	   kx(im,j)=k(im,j)
	end do
   !=========y方向界面导热系数==========!
	do j=3,jmm
	   do i=2,imm
	      ky(i,2)=k(i,1)
	      ky(i,jm)=k(i,jm)
	      ky(i,j)=k(i,j)*k(i,j-1)*dy(j)/(k(i,j)*yvs(j)+k(i,j-1)*yvn(j))   
	   end do
	end do   
 return
 end subroutine


 subroutine output()
    use variable
    implicit none
    !character(len=20)::Tnode
    integer::i,j
    open(unit=20,file='Tnode.dat')
	!open(unit=10,file='axx.dat')
 
    write(20,*) 'TITLE="line"'
    write(20,*) 'VARIABLES="x","y","T"'
    write(20,*)'ZONE x=',im,'y=',jm,' T'
	!write(10,*)'ZONE x=',im,'y=',jm,' ae',' aw',' an',' as'
    do j=1,jm
	   do i=1,im
		  write(20,*)x(i),y(j),T(i,j)
          !write(10,*)x(i),y(j),node(i,j)%ae,node(i,j)%aw,node(i,j)%an,node(i,j)%as
	   end do
	end do
  return
end subroutine

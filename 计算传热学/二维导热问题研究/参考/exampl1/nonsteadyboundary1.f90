!This is a program for calculationg the heat conduction of a wall.
! Formulations in detail can be found in Book < Numerical Transfer & Heat Conduction> 
! by Patankar.
! In this program,we use uniform grid to compute unsteady changing temperature 
! of a wall.
! We compare two numerical methods : point iteration and block iteration.

program nonsteady_heat_transfer_boundary_conditon_1
    use variable
    implicit none
	integer::i,j,n 
    call initial()
	call grid()
    call start()  
    do while(.not.Converged_all) 
		! 边界赋值或更新源项内容.
        call boundary()
        ! Compute coefficient ap,ae.... for each nodes.
        call compute()
        ! Compute temperature field at each time step.  
        call solve()
        ! Check whether achieve a steady state ( Converged_all )
        call check_all() 
  
      if (dt<10) then
	   	 if (Converged_all) then
		    if (abs(time-timeend).le.dt) then
		       goto 110
		    else
	          Converged_all=.false.
	          time=time+dt
		      tlevel=tlevel+1
			 ! write(filename,"(i5,'.dat')")tlevel
		     ! call output()
		    end if   
		  end if
	  else
	      if (Converged_all) then
	          filename="Tnode.dat"
		      call output()
			 write(*,*)'After ',tlevel,' iterations achieve converged '
   		  else
		     tlevel=tlevel+1	  
	      end if
	  end if
    end do 
110 filename="T15.dat"
    call output()
end program



subroutine initial()
    use variable
    implicit none
    integer::i,j
    xl=1 
    yl=1
	mode=1
	Cp=800
	rou=8000.0
	timeend=1.0
	Tfw=50
	Tfe=100
	T0=0
	hw=5
	he=10
	k0=20
    R0=1 ! 直角坐标R0=1，柱坐标R0为r方向初始半径
    eps=1E-5
    dt=0.5
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
	     	k(i,j)=k0  !第(i,j)控制容积导热系数
		   if (i==im) then
		   	   TN(i,j)=100.0*(1.0+y(j))
               TN_old(i,j)=100.0*(1.0+y(j))
	           T(i,j)=100.0*(1.0+y(j))
			else
			   TN(i,j)=T0
               TN_old(i,j)=T0
		  	   T(i,j)=T0
			end if
		end do
	end do
   write(filename,"(i5,'.dat')")time
   call output()
return
end subroutine


subroutine boundary()
    use variable
    implicit none
    integer::i,j
	real::temph,Spad,Scad,areax,areays
	
!========控制容积内源项处理========!
	do j=2,jmm
	    do i=2,imm
		    Sc(i,j)=0
			Sp(i,j)=0
			b(i,j)=0;  !控制容积内部源项
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

!=========绝热边界条件上下左边界节点值更新==========!
	do i=2,imm
	   T(i,1)=T(i,2)
	   T(i,jm)=T(i,jmm)
	end do
    do j=2,jmm
	   i=1
	   T(i,j)=T(i+1,j)
	end do

 return
 end subroutine

 subroutine output()
    use variable
    implicit none
    integer::i,j
	character(len=20)::name
	name=filename

   ! T(1,1)=T(1,2)+T(2,1)-T(2,2)
	!T(1,jm)=T(1,jmm)+T(2,jm)-T(jmm,jmm)
!	T(im,1)=T(imm,1)+T(im,2)-T(imm,imm)
	!T(im,jm)=T(im,jmm)+T(imm,jm)-T(imm,jmm)

    open(unit=20,file=name)
    write(20,'(2Xa)') 'TITLE="Animate"'
    write(20,100) 'VARIABLES="x","y","T",'
    write(20,200) 'ZONE I=',im, 'J=',jm,'F=POINT'
100 format(1xa22)
200 format(a8,i3,6xa3,i3,3xa10)
300 format(3xf6.4,6xf6.4,8xf8.4)
    do j=1,jm
	   do i=1,im
		  write(20,300)x(i),y(j),T(i,j)
	   end do
	end do
  return
 ! close(20)
end subroutine

!****************************************************************************
!
!  PROGRAM: Second_heat_transfer
!
!  PURstringE:   数值计算第二题_稳态，无源项混合流动
!****************************************************************************
!一阶格式的
program Second_heat_transfer
use head
implicit none
!---------------------------定义外部函数
real,external::A
real,external::D
real,external::F
real,external::Den_ve
real,external::APP
integer:: i,j
!-----------------------------
M=60                      !x方向网格数量
N=60                     !y方向网格数量
dx=L/M
dy=H/N

!---------------------------分配动态数组
allocate(xf(M+2),yf(N+2))
allocate(x(M+1,N+1),y(M+1,N+1),u(M+1,N+1),v(M+1,N+1),temp(M+1,N+1),temp2(M+1,N+1),Px(M+1,N+1),Py(M+1,N+1))
allocate(ae(M+1,N+1),aw(M+1,N+1),as(M+1,N+1),an(M+1,N+1),ap(M+1,N+1))

allocate(temp1(M+1,N+1),temp5(M+1,N+1),tempfalse(M+1,N+1))

temp=300                      !温度场初始化
!----------------------------


       call coord()         !初始化节点坐标和控制面坐标
       call velocity()      !初始化节点速度
       call boundary()      !设置温度边界条件
       call coef()          !求解系数
!       call coef2()          !求解二阶系数
       call intera()        !迭代求解
!	   call temp_false()     !假扩散

!       Px1=density*u(1,N/2)*dx/lamda			!对于层流的还适用
!	   Px2=density*u(1,2*N/3+3)*dx/lamda
!	   write(*,*)Px1,Px2
		
		Pxmax=0;Pymax=0;Pxjun=0;Pyjun=0
		do i=1,M+1
			do j=1,N+1
				Px(i,j)=density*u(i,j)*dx/lamda
				Py(i,j)=density*v(i,j)*dy/lamda
				Pxjun=Pxjun+abs(Px(i,j))
				Pyjun=Pyjun+abs(Py(i,j))
				if(abs(Px(i,j))>Pxmax)then
					Pxmax=Px(i,j)
				end if
				if(abs(Py(i,j))>Pymax)then
					Pymax=Py(i,j)
				end if
			end do
		end do
		Pxjun=Pxjun/((M+1)*(N+1))
		Pyjun=Pyjun/((M+1)*(N+1))
		write(*,*)Pxmax,Pymax,Pxjun,Pyjun
 

deallocate(xf,yf)
deallocate(x,y,u,v,temp,temp2,ae,aw,as,an,ap,Px,Py)
deallocate(temp1,temp5,tempfalse)


end program Second_heat_transfer


!***************************初始化坐标
subroutine coord()
  use head
  integer::i,j                         !循环变量
!--------------------初始化节点坐标----------------- 
 do i=1,M+1
     do j=1,N+1
	   x(i,j)=(i-1)*dx
	   y(i,j)=(j-1)*dy 
     enddo
 enddo

!--------------------初始化控制面坐标---------------
 xf(1)=0.0
 xf(2)=(x(2,1)+x(1,1))*0.5
 xf(M+2)=x(M+1,1)
 xf(M+1)=(x(M+1,1)+x(M,1))*0.5
 yf(1)=0.0
 yf(2)=(y(1,2)+y(1,1))*0.5
 yf(N+2)=y(1,N+1)
 yf(N+1)=(y(1,N+1)+y(1,N))*0.5
!----------------------------------------------
  do i=3,M
      xf(i)=(i-1)*dx-0.5*dx
  enddo
  do j=3,N
      yf(j)=(j-1)*dy-0.5*dy
  enddo
!----------------------------------------------------------
end subroutine

!**************************初始化速度
subroutine velocity()
 use head
  integer::i,j                         !循环变量
    open(unit=10,file='velocity-mesh1.txt')
      open(unit=102,file='velocity-mesh1.dat')
         write(102,*)"Variables= x,y,u,v"  
	      write(102,*)''
           write(102,*) "ZONE I=",N+1,"J=",M+1,"F=point"
  	  do i=1,M+1
          do j=1,N+1
 	        read(10,*) x1,y1,u(i,j),v(i,j)
			write(*,*)i,j
 		    write(102,*)  x1,y1,u(i,j),v(i,j)
 	     enddo
       enddo
  close(10)
	close(102)

 end subroutine

!subroutine velocity()
! use head
!  integer::i,j 
!   open(unit=102,file='velocity.txt')
!	do i=1,M+1
!		do j=1,N/3
!			u(i,j)=1
!			v(i,j)=0
!			write(102,*)  x(i,j),y(i,j),u(i,j),v(i,j)
!		end do 
!		do j=N/3+1,2*N/3+1
!			u(i,j)=2
!			v(i,j)=0
!			write(102,*)  x(i,j),y(i,j),u(i,j),v(i,j)
!		end do 

!		do j=2*N/3+2,N+1
!			u(i,j)=1
!			v(i,j)=0
!			write(102,*)  x(i,j),y(i,j),u(i,j),v(i,j)
!		end do 
!	end do

!	close(102)
!end subroutine

!******************************设置温度边界条件
subroutine boundary()
  use head
 integer::i,j                         !循环变量
  do i=2,M+1
     temp(i,1)= 300 !27     !300
	 temp(i,N+1)=300 !57   !300	 
  enddo
	 

  do j=1,N+1
	 temp(1,j)= 450 !57     !330                !温度的单位可能要调整一下		
  enddo 

end subroutine

!*******************************求解系数
subroutine coef()
 use head
 integer::i,j                         !循环变量
	ae=0.0;aw=0.0;as=0.0;an=0.0
	  do i=2,M
	      do j=2,N
		      ae(i,j)=A(i,j,'e') 		 
			  as(i,j)=A(i,j,'s')
		      aw(i,j)=A(i,j,'w')
			  an(i,j)=A(i,j,'n') 
			  ap(i,j)=(ae(i,j)+as(i,j)+aw(i,j)+an(i,j))+(F(i,j,'e')-F(i,j,'w'))+(F(i,j,'n')-F(i,j,'s'))
		 enddo
	  enddo
	 
	 do j=2,N
	    ae(M+1,j)=0
		as(M+1,j)=A(M+1,j,'s')
		aw(M+1,j)=A(M+1,j,'w')
		an(M+1,j)=A(M+1,j,'n') 
		ap(M+1,j)=(ae(M+1,j)+as(M+1,j)+aw(M+1,j)+an(M+1,j))+(F(M+1,j,'e')-F(M+1,j,'w'))+(F(M+1,j,'n')-F(M+1,j,'s'))	

	enddo       
 
end subroutine

!***********函数A
 function A(i,j,string)

  use head
  real::A
  integer::i,j   
  character::string
  ff=F(i, j, string)
  dd=D(i, j, string)
 
 	if (string == 'e'.OR. string == 'n')then 	
	   A =dd* APP(ff/dd )+max(-ff, 0.0);
	elseif (string == 'w'.OR. string == 's')then
 	   A =dd* APP(ff/dd)+max(ff, 0.0);
    endif
    return;
 end
!************函数D
function D(i,j,string)
  use head
  real::D
  integer::i,j   
  character::string    
	if(string == 'n')then
	   
	   D= lamda*(xf(i+1) - xf(i)) / (y(i,j+1) - y(i,j));
	elseif(string == 's')then
	   D= lamda*(xf(i+1) - xf(i)) / (y(i,j) - y(i,j-1));
	elseif(string== 'w')then
	   D= lamda*(yf(j+1) - yf(j)) / (x(i,j) - x(i-1,j));
	elseif(string == 'e')then
	   D= lamda* (yf(j+1) - yf(j)) / (x(i+1,j) - x(i,j));
	endif 

end
!************函数F
function F(i,j,string)
  use head
  real::F
  integer::i,j   
  character::string
   
	if(string == 'n') then
	   F= Den_ve (i, j, 'n') * (xf(i+1) - xf(i));
	elseif(string == 's')then
	   F= Den_ve (i, j, 's') * (xf(i+1) - xf(i));
	elseif(string == 'w')then
	   F= Den_ve (i, j, 'w') * (yf(j+1) - yf(j));
	elseif(string == 'e')then
	   F= Den_ve (i, j, 'e') * (yf(j+1) - yf(j));
    endif
end
!************函数Den_ve density*velocity
 function Den_ve(i,j,string)
  use head
  real::Den_ve
  integer::i,j   
  character::string
  !print*,i,j
    if(string == 'n')then	 
		Den_ve=(density* v(i,j) * (y(i,j+1)-yf(j+1) ) + density * v(i,j+1)*(yf(j+1) -y(i,j)))/ (y(i,j+1)-y(i,j));	 
	elseif(string == 's')then 
	    Den_ve=(density* v(i,j) * (yf(j)-y(i,j-1) ) + density* v(i,j-1) * (y(i,j) - yf(j)))/ (y(i,j)-y(i,j-1)); 
	elseif(string == 'w')then 
		Den_ve=(density* u(i,j) * (xf(i)-x(i-1,j) ) + density* u(i-1,j) * (x(i,j) - xf(i)))/ (x(i,j)-x(i-1,j)); 
	elseif(string == 'e')then 
		if(i==M+1)then
		Den_ve=density* u(M+1,j)
		else
		Den_ve=(density* u(i,j) * (x(i+1,j)-xf(i+1) ) + density* u(i+1,j) * (xf(i+1) -x(i,j)))/ (x(i+1,j)-x(i,j));
        endif
	endif
 
end
!********************APP--计算A(|p|)---根据不同的choose
 function APP(p)
   use head
   real::APP
   real::p
	if( choose == 1)then      ! 中心差分格式
	    APP=1 - 0.5 * abs(p);       
	elseif(choose == 2)then   ! 迎风格式
	    APP= 1;                   
	elseif(choose == 3)then   ! 指数格式    
	    APP=(abs(p)+1e-8)/(exp(abs(p) + 1e-8) - 1.0);
	elseif(choose == 4) then  ! 混合格式
	    APP=max(1-0.5*abs(p),0.0);
	elseif(choose == 5)then ! 乘方格式
        APP=max(0.0,(1.0-0.1*(abs(p)))**5)
	endif
	return	
 end


subroutine intera()    !用G-S迭代法求解ap*Tp=ae*Te+aw*Tw+an*Tn+as*Ts
 use head
 	integer i, j
    real:: maxTemp = 0.0, maxDelta = 0.0
	real:: tempTemp, deltaTemp
	real:: value
	integer::time=0

do while(epsilon1.GT.epsilon.AND.time.LT.mostinter.AND.epsilon1.LT.1e5)
    time=time+1
	maxdelta=0.0
    temp2=temp                                 !更新温度场
  do i=2,M+1
      do j=2,N
 	    if(temp(i,j).GT.maxTemp)then
			 maxTemp = temp(i,j);              !   //温度最大值
		endif
	enddo
  enddo

	do j=2,N
	   do i=2,M	   
	    value = ae(i,j) * temp(i+1,j) + aw(i,j) * temp(i-1,j)+(an(i,j) * temp(i,j+1) + as(i,j) * temp(i,j-1))		
		temp(i,j) = value / ap(i,j)			   !  //TP=value/ap ,新的到的Tp,一次循环过后，所有温度都将被更新一遍	
	   enddo

       value = 0 + aw(M+1,j) * temp(M,j)+(an(M+1,j) * temp(M+1,j+1) + as(M+1,j) * temp(M+1,j-1))
	   temp(M+1,j) = value / ap(M+1,j)
	end do

	do i=2,M+1
		do j=2,N
				deltaTemp =abs(temp(i,j)- temp2(i,j));  	 
				if(deltaTemp.GT.maxDelta)then !   // 如果两次温度差大于最大温度差，则将本次循环得到的温度差赋给最大温度差
				maxDelta = deltaTemp;   
 				endif  
		enddo
	enddo

!	print*,maxDelta

	epsilon1 = maxDelta ! / maxTemp;                              ! //精度为：两次温度差/最大温度
    residual(time)=epsilon1                     !残差
	
enddo

print*,time

   open(unit=101,file='result-lamda-e-1-geshi5-mesh1.dat')
      write(101,*)"Variables= x,y,u,v,temp"
	  write(101,*)''
    write(101,*) "ZONE I=",N+1,"J=",M+1,"F=point"
       do i=1,M+1
          do j=1,N+1
	 		write(101,*) x(i,j),y(i,j),u(i,j),v(i,j),temp(i,j) 
 	      enddo
      enddo
  close(101)

open(unit=111,file='residual-lamda-e-1-geshi5-mesh1.dat')          !残差
       do i=1,time
         write(111,*) i,residual(i) 
      enddo
  close(111)


open(unit=201,file='temp-lamda-e-1-geshi5-mesh1.txt')                   !温度
    do i=1,M+1
	  do j=1,N+1
	     write(201,*) x(i,j),y(i,j),temp(i,j)
	  enddo
	enddo
close(201)


end subroutine
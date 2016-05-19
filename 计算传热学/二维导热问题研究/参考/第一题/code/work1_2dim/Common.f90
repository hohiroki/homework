Module Varieble



real(8)				::	TB(4),qB(4),h(4),tf(4)	!边界条件取值
integer(4)			::	i,j,k,M,N,Nt			!循环变量和网格数


!-------------------------------------
!四个边界条件的类型:1：第一类边界；2：第二类边界；3：第三类边界
integer(4)			::	bcon(4)					
real(8)				::	x1,x2,y1,y2,epsilon	!区域范围，迭代残差
real(8),allocatable	::	x(:),y(:),dx(:),dy(:),del_x(:),del_y(:),lambda(:,:),lambda_we(:,:),lambda_ns(:,:)	!空间步长，热导率
!------------------------------------------------------------------------------------------
!ET:精确解；ae,aw,an,as,ap,ap0:差分方程的系数；S：源项
real(8),allocatable	::	ET(:,:),ae(:,:),aw(:,:),an(:,:),as(:,:),ap(:,:),ap0(:,:)
!--------------------------------------------------------------------------------
!rn,rs：圆柱坐标的控制变量；Sc,Sp：源项线性化的系数；Scad,Spad：附加源项的系数
real(8),allocatable	::	rn(:),rs(:),Sc(:,:),Sp(:,:),Scad(:,:),Spad(:,:)

!-------------------------------------------------------------------
!time：计算时间,  dt：时间步长,   rho：密度,  Cp比热,   TT0：初始温度
!theta：显式格式（0）；隐式格式（1）；加权隐式（0—1）；定常问题（2）
real(8)				::	time,dt,rho,Cp,TT0,theta
!character(50)		::	cond


!-----------------------------------------------------------
!fname：数值解；    fexact：精确解         ferror：与精确解的误差
character(50)		::	fname,fexact,ferror


!----------------------------------------------------------------
!exact_symbol：精确解控制变量；animate：输出多组数据的控制变量；
!dimen:问题的维度；coordinate：坐标系类型
integer(4)			::	exact_symbol,animate,dimen,coordinate	
	

!------------Case7--------------------------------------------
real(8)		::	delta,q0	!Case7参数：气缸壁厚，热流大小
integer(4)	::	num

!------------空心区域范围----------------------------------------
integer(4),allocatable		::	zoom(:,:)
integer(4)	::	zoom_bcon
real(8)		::	zoom_TB

end module 
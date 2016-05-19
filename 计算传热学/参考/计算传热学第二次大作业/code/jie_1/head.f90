module head 
!--------------------------------------
	real,allocatable::xf(:),yf(:),x(:,:),y(:,:),u(:,:),v(:,:),temp(:,:),temp2(:,:)
	real,allocatable::ae(:,:),aw(:,:),as(:,:),an(:,:),ap(:,:),Px(:,:),Py(:,:)

	real,allocatable::temp1(:,:),temp5(:,:),tempfalse(:,:)

	integer::M,N
 !--------------------------------------- 
	real::epsilon=1E-4,epsilon1=5.0       !迭代精度
	real::residual(1e6)                   !残差
	integer::mostinter=200000             !最大迭代步数

	real::H=0.3,L=0.6
	real::dx,dy,Pxmax,Pymax,Pxjun,Pyjun

	real::density=1.225          !kg/m3
	real::lamda=2.5E-2    !0.58   !2.5E-2     !   4 摄氏度

	real::choose=5              !离散方程类型1~5         !1-中心差分;2-迎风格式;3-指数格式;4-混合格式;5-乘方格式

	real::ff,dd
	real::x1,y1
	end
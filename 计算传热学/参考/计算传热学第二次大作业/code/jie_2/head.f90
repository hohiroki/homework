module head 
!--------------------------------------
	real,allocatable::xf(:),yf(:),x(:,:),y(:,:),u(:,:),v(:,:),temp(:,:),temp2(:,:),Px(:,:),Py(:,:)
	real,allocatable::ae(:,:),aw(:,:),as(:,:),an(:,:),ap(:,:)

	real,allocatable::temp1(:,:),temp5(:,:),tempfalse(:,:)

	integer::M,N
 !--------------------------------------- 
	real::epsilon=1E-4,epsilon1=5.0       !��������
	real::residual(1e6)                   !�в�
	integer::mostinter=200000             !����������

	real::H=0.3,L=0.6
	real::dx,dy,Pxmax,Pymax

	real::density=1.225          !kg/m3
	real::lamda=2.5E-2    !0.58   !2.5E-2     !   4 ���϶�

	real::choose=2               !��ɢ��������1~5         !1-���Ĳ��;2-ӭ���ʽ;3-ָ����ʽ;4-��ϸ�ʽ;5-�˷���ʽ

	real::ff,dd
	real::x1,y1
	end

	!2��ӭ��
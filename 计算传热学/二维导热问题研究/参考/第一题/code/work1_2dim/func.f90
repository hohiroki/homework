



!===============Boundary：确定边界差分方程系数==============================
	subroutine boundary(ap1,ae1,aw1,an1,as1,ap10,tSc,tSp,bn,tbcon,rn1,rs1,dxy)
	use Varieble
	implicit none
	integer(4)		::	bn,tbcon
	real(8)			::	ae1,ap1,aw1,an1,as1,ap10,tSc,tSp,rn1,rs1,dxy
	real(8)			::	Sca,Spa
	call con(Sca,Spa,bn,tbcon)


	if(bn==1)then
		tSc=tSc+Sca*dxy*rs1
		ap1=ae1+aw1+an1-Spa*dxy*rs1+ap10-tSp
	elseif(bn==2)then
		tSc=tSc+Sca*dxy*(rn1+rs1)/2
		ap1=as1+aw1+an1-Spa*dxy*(rn1+rs1)/2+ap10-tSp
	elseif(bn==3)then
		tSc=tSc+Sca*dxy*rn1
		ap1=ae1+aw1+as1-Spa*dxy*rn1+ap10-tSp
	elseif(bn==4)then
		tSc=tSc+Sca*dxy*(rn1+rs1)/2
		ap1=as1+ae1+an1-Spa*dxy*(rn1+rs1)/2+ap10-tSp
	else
		print*,"Error in Boundary(bn)!"
	endif		
	end subroutine
!===============================================================


!===========conditions：不同边界对应的附加源项==================================
	subroutine con(Sc1,Sp1,bn,tbcon)
	use Varieble
	implicit none
	real(8)		::	Sc1,Sp1
	integer(4)	::	bn,tbcon
	if (tbcon==1) then
		Sc1=1e30*TB(bn)
		Sp1=-1e30
	elseif(tbcon==2) then
		Sc1=qB(bn)
		Sp1=0
	elseif(tbcon==3) then
		Sc1=h(bn)*Tf(bn)
		Sp1=-h(bn)
	else
		print *, "Error in con(bn)!"
	endif
	end subroutine
!===============================================================
	
	
!===============corner：角点的差分方程的系数=========================================
	subroutine corner(n1)
	use Varieble
	implicit none
	integer(4)	::	n1
	real(8)		::	Sc1,Sp1,Sc2,Sp2

	if(n1==1)then
		Call con(Sc1,Sp1,4,bcon(4))
		Call con(Sc2,Sp2,1,bcon(1))
		ae(1,1)=dy(1)/del_x(1)*lambda_we(1,1)*(rn(1)+rs(1))/2
		an(1,1)=dx(1)/del_y(1)*lambda_ns(1,1)*rn(1)
		aw(1,1)=0
		as(1,1)=0
		ap(1,1)=ae(1,1)+an(1,1)-Sp1*dy(1)*(rn(1)+rs(1))/2-Sp2*dx(1)*rs(1)+ap0(1,1)-Sp(1,1)
		Sc(1,1)=Sc(1,1)+Sc1*dy(1)*(rn(1)+rs(1))/2.+Sc2*dx(1)*rs(1)
	elseif(n1==2)then
		Call con(Sc1,Sp1,1,bcon(1))
		Call con(Sc2,Sp2,2,bcon(2))
		aw(M+1,1)=dy(1)/del_x(M)*lambda_we(M,1)*(rn(1)+rs(1))/2
		an(M+1,1)=dx(M+1)/del_y(1)*lambda_ns(M+1,1)*rn(1)
		ae(M+1,1)=0
		as(M+1,1)=0
		ap(M+1,1)=aw(M+1,1)+an(M+1,1)-Sp1*dy(1)*(rn(1)+rs(1))/2-Sp2*dx(M+1)*rs(1)+ap0(M+1,1)-Sp(M+1,1)
		Sc(M+1,1)=Sc(M+1,1)+Sc1*dy(1)*(rn(1)+rs(1))/2+Sc2*dx(M+1)*rs(1)
	elseif(n1==3)then
		Call con(Sc1,Sp1,2,bcon(2))
		Call con(Sc2,Sp2,3,bcon(3))
		aw(M+1,N+1)=dy(N+1)/del_x(M)*lambda_we(M,N+1)*(rn(N+1)+rs(N+1))/2
		as(M+1,N+1)=dx(M+1)/del_y(N)*lambda_ns(M+1,N)*rs(N+1)
		ae(M+1,N+1)=0
		an(M+1,N+1)=0
		ap(M+1,N+1)=aw(M+1,N+1)+as(M+1,N+1)-Sp1*dy(N+1)*(rn(N+1)+rs(N+1))/2-Sp2*dx(M+1)*rn(N+1)+ap0(M+1,N+1)-Sp(M+1,N+1)
		Sc(M+1,N+1)=Sc(M+1,N+1)+Sc1*dy(N+1)*(rn(N+1)+rs(N+1))/2+Sc2*dx(M+1)*rn(N+1)
	elseif(n1==4)then
		Call con(Sc1,Sp1,3,bcon(3))
		Call con(Sc2,Sp2,4,bcon(4))
		ae(1,N+1)=dy(N+1)/del_x(1)*lambda_we(1,N+1)*(rn(N+1)+rs(N+1))/2
		as(1,N+1)=dx(1)/del_y(N)*lambda_we(1,N)*rn(N+1)
		aw(1,N+1)=0
		an(1,N+1)=0
		ap(1,N+1)=ae(1,N+1)+as(1,N+1)-Sp1*dy(N+1)*(rn(N+1)+rs(N+1))/2-Sp2*dx(1)*rn(N+1)+ap0(1,N+1)-Sp(1,N+1)
		Sc(1,N+1)=Sc(1,N+1)+Sc1*dy(N+1)*(rn(N+1)+rs(N+1))/2+Sc2*dx(1)*rn(N+1)
	else
		print*,"Error in corner(n)!"
	endif		
	end subroutine
!==============================================================


	
	
!======================Cal_a：计算差分方程系数================================
	subroutine cal_a
	use Varieble

	if(theta==2)then
		ap0=0
	else 
	do j=1,N+1
		do i=1,M+1
			ap0(i,j)=rho*Cp*dx(i)*dy(j)/dt*(rn(j)+rs(j))/2
		enddo
	enddo
	endif

	do j=2,N
		do i=2,M
			ae(i,j)=dy(j)/del_x(i)*lambda_we(i,j)*(rn(j)+rs(j))/2
			aw(i,j)=dy(j)/del_x(i-1)*lambda_we(i-1,j)*(rn(j)+rs(j))/2
			an(i,j)=dx(i)/del_y(j)*lambda_ns(i,j)*rn(j)
			as(i,j)=dx(i)/del_y(j-1)*lambda_ns(i,j-1)*rs(j)
			ap(i,j)=ae(i,j)+aw(i,j)+an(i,j)+as(i,j)+ap0(i,j)-Sp(i,j)
		enddo
	enddo

		do i=2,M
			ae(i,1)=dy(1)/del_x(i)*lambda_we(i,1)*(rn(1)+rs(1))/2
			aw(i,1)=dy(1)/del_x(i-1)*lambda_we(i-1,1)*(rn(1)+rs(1))/2
			an(i,1)=dx(i)/del_y(1)*lambda_ns(i,1)*rn(1)
			as(i,1)=0
		enddo

		do i=2,M
			ae(i,N+1)=dy(N+1)/del_x(i)*lambda_we(i,N+1)*(rn(N+1)+rs(N+1))/2
			aw(i,N+1)=dy(N+1)/del_x(i-1)*lambda_we(i-1,N+1)*(rn(N+1)+rs(N+1))/2
			as(i,N+1)=dx(i)/del_y(N)*lambda_ns(i,N)*rn(N+1)
			an(i,N+1)=0
		enddo

	do j=2,N
			ae(1,j)=dy(j)/del_x(1)*lambda_we(1,j)*(rn(j)+rs(j))/2
			an(1,j)=dx(1)/del_y(j)*lambda_ns(1,j)*rn(j)
			as(1,j)=dx(1)/del_y(j-1)*lambda_ns(1,j-1)*rs(j)
			aw(1,j)=0
	enddo

	do j=2,N
			aw(M+1,j)=dy(j)/del_x(M)*lambda_we(M,j)*(rn(j)+rs(j))/2
			an(M+1,j)=dx(M+1)/del_y(j)*lambda_ns(M+1,j)*rn(j)
			as(M+1,j)=dx(M+1)/del_y(j-1)*lambda_ns(M+1,j-1)*rs(j)
			ae(M+1,j)=0
	enddo




!	S=0
	!-----------j=1---------------------------------------------
	Call corner(1)
	Call corner(2)
	do i=2,M
		Call boundary(ap(i,1),ae(i,1),aw(i,1),an(i,1),as(i,1),ap0(i,1),Sc(i,1),Sp(i,1),1,bcon(1),rn(1),rs(1),dx(i))
	enddo
	!------------j=2-N-----------------
		do j=2,N
			Call boundary(ap(1,j),ae(1,j),aw(1,j),an(1,j),as(1,j),ap0(1,j),Sc(1,j),Sp(1,j),4,bcon(4),rn(j),rs(j),dy(j))
			Call boundary(ap(M+1,j),ae(M+1,j),aw(M+1,j),an(M+1,j),as(M+1,j),ap0(M+1,j),Sc(M+1,j),Sp(M+1,j),2,bcon(2),rn(j),rs(j),dy(j))
		enddo
		!------------j=N+1-----------------
		Call corner(4)
		Call corner(3)
		do i=2,M
			Call boundary(ap(i,N+1),ae(i,N+1),aw(i,N+1),an(i,N+1),as(i,N+1),ap0(i,N+1),Sc(i,N+1),Sp(i,N+1),3,bcon(3),rn(N+1),rs(N+1),dx(i))
		enddo
		
!		print*,ap(2,2),aw(2,2),ae(2,2),an(2,2),as(2,2),ap0(2,2)

!		Call output(ap,'ap.txt')
!		pause
	end subroutine
!=========================================================================================================


!======================TDMA算法=====================================
	subroutine TDMA(T,a,b,c,d)
	use Varieble
	implicit none
	real(8)			::	T(M+1),a(M+1),b(M+1),c(M+1),d(M+1)
	real(8)			::	P(M+1),Q(M+1)
	p(1)=b(1)/a(1)
	Q(1)=d(1)/a(1)
	do i=2,M+1
		P(i)=b(i)/(a(i)-c(i)*P(i-1))
		Q(i)=(c(i)*Q(i-1)+d(i))/(a(i)-c(i)*P(i-1))
	enddo
	T(M+1)=Q(M+1)
	do i=M,1,-1
         T(i)=P(i)*T(i+1)+Q(i)
	enddo
	end subroutine	
!===================================================================


!======================Calulate r：计算圆柱坐标的r===============================
	subroutine cal_r
	use Varieble
	implicit none
	if(coordinate==1) then   !直角坐标
		rn=1
		rs=1
	elseif(coordinate==2) then !圆柱坐标
		rn(1)=y(1)+dy(1)
		rs(1)=y(1)
		rn(N+1)=y(N+1)
		rs(N+1)=y(N+1)-dy(N+1)
		do j=2,N
			rs(j)=rn(j-1)
			rn(j)=rs(j)+dy(j)
		enddo
	else
		print*, 'Error in cal_r!'
	endif
	end subroutine
!======================================================================
program te
use test
implicit none
real(8)		::	a=0.1
i=1
print*,i
call printf
end program



!===============Boundary==============================
	subroutine boundary(ap1,ae1,aw1,an1,as1,ap10,b,bn,tbcon)
	use Varieble
	implicit none
	integer(4)		::	bn,tbcon
	real(8)			::	ae1,ap1,aw1,an1,as1,ap10,b
	real(8)			::	Sca,Spa
	call con(Sca,Spa,bn,tbcon)
	if(cond=='Steady') then
		ap10=0
	else if(cond=='Explicit') then
		ap10=rho*Cp*dx*dy/dt/2
	else
		print*, 'Error in con!'
	endif
	if(bn==1)then
		ae1=lambda*dy/dx/2
	    aw1=ae1
	    an1=lambda*dx/dy
	    as1=0
		b=Sca*dx
		ap1=ae1+aw1+an1-Spa*dx+ap10
	elseif(bn==2)then
		aw1=lambda*dy/dx
	    ae1=0
	    an1=lambda*dx/dy/2
	    as1=an1
		b=Sca*dy
		ap1=as1+aw1+an1-Spa*dy+ap10
	elseif(bn==3)then
		ae1=lambda*dy/dx/2
	    aw1=ae1
	    as1=lambda*dx/dy
	    an1=0
		b=Sca*dx
		ap1=ae1+aw1+as1-Spa*dx+ap10
	elseif(bn==4)then
		ae1=lambda*dy/dx
	    aw1=0
	    an1=lambda*dx/dy/2
	    as1=an1
		b=Sca*dy
		ap1=as1+ae1+an1-Spa*dy+ap10
	else
		print*,"Error in Boundary(bn)!"
	endif		
	end subroutine
!===============================================================


!===================conditions==================================
	subroutine con(Sc,Sp,bn,tbcon)
	use Varieble
	implicit none
	real(8)		::	Sc,Sp
	integer(4)	::	bn,tbcon
	if (tbcon==1) then
		Sc=1e30*TB(bn)
		Sp=-1e30
	elseif(tbcon==2) then
		Sc=qB(bn)
		Sp=0
	elseif(tbcon==3) then
		Sc=h(bn)*Tf(bn)
		Sp=-h(bn)
	else
		print *, "Error in con(bn)!"
	endif
	end subroutine
	!===============================================================
	
	
	!===============corner=========================================
	subroutine corner(ap1,ae1,aw1,an1,as1,ap10,b,n1)
	use Varieble
	implicit none
	real(8)		::	ap1,ae1,aw1,an1,as1,ap10,b
	integer(4)	::	n1
	real(8)		::	Sc1,Sp1,Sc2,Sp2


	if(cond=='Steady')then
		ap10=0
	else if(cond=='Explicit') then
		ap10=rho*Cp*dx*dy/dt/4
	else
		print*, 'Error in con!'
	endif

	if(n1==1)then
		Call con(Sc1,Sp1,4,bcon(4))
		Call con(Sc2,Sp2,1,bcon(1))
		ae1=lambda*dy/dx/2
		aw1=0
		an1=lambda*dx/dy/2
		as1=0
	elseif(n1==2)then
		Call con(Sc1,Sp1,1,bcon(1))
		Call con(Sc2,Sp2,2,bcon(2))
		ae1=0
		aw1=lambda*dy/dx/2
		an1=lambda*dx/dy/2
		as1=0
	elseif(n1==3)then
		Call con(Sc1,Sp1,2,bcon(2))
		Call con(Sc2,Sp2,3,bcon(3))
		ae1=0
		aw1=lambda*dy/dx/2
		an1=0
		as1=lambda*dx/dy/2
	elseif(n1==4)then
		Call con(Sc1,Sp1,3,bcon(3))
		Call con(Sc2,Sp2,4,bcon(4))
		ae1=lambda*dy/dx/2
		aw1=0
		an1=0
		as1=lambda*dx/dy/2
	else
		print*,"Error in corner(n)!"
	endif		
	ap1=ae1+aw1+an1+as1-Sp1*dy/2-Sp2*dx/2+ap10
	b=Sc1*dy/2+Sc2*dx/2
	end subroutine
!==============================================================


	
	
!======================Cal_a================================
	subroutine cal_a
	use Varieble

	if(cond=='Steady')then
		ap0=0
	else if(cond=='Explicit') then
		ap0=rho*Cp*dx*dy/dt
	else
		print*, 'Error in con!'
	endif

	ae=dy/dx*lambda
	aw=ae
	an=dx/dy*lambda
	as=an
	ap=ae+aw+an+as+ap0

	S=0
	!-----------j=1---------------------------------------------
	Call corner(ap(1,1),ae(1,1),aw(1,1),an(1,1),as(1,1),ap0(1,1),S(1,1),1)
	Call corner(ap(M+1,1),ae(M+1,1),aw(M+1,1),an(M+1,1),as(M+1,1),ap0(M+1,1),S(M+1,1),2)
	do i=2,M
		Call boundary(ap(i,1),ae(i,1),aw(i,1),an(i,1),as(i,1),ap0(i,1),S(i,1),1,bcon(1))
	enddo
	!------------j=2-N-----------------
		do j=2,N
			Call boundary(ap(1,j),ae(1,j),aw(1,j),an(1,j),as(1,j),ap0(1,j),S(1,j),4,bcon(4))
			Call boundary(ap(M+1,j),ae(M+1,j),aw(M+1,j),an(M+1,j),as(M+1,j),ap0(M+1,j),S(M+1,j),2,bcon(2))
		enddo
		!------------j=N+1-----------------
		Call corner(ap(1,N+1),ae(1,N+1),aw(1,N+1),an(1,N+1),as(1,N+1),ap0(1,N+1),S(1,N+1),4)
		Call corner(ap(M+1,N+1),ae(M+1,N+1),aw(M+1,N+1),an(M+1,N+1),as(M+1,N+1),ap0(M+1,N+1),S(M+1,N+1),3)
		do i=2,M
			Call boundary(ap(i,N+1),ae(i,N+1),aw(i,N+1),an(i,N+1),as(i,N+1),ap(i,N+1),S(i,N+1),3,bcon(3))
		enddo
	end subroutine
!=========================================================================================================


!======================TDMA=====================================
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

	



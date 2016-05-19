module variable
    implicit none
    integer,parameter::im=7,jm=7
	integer,parameter::imm=im-1,jmm=jm-1
	integer::solve_type
    real,dimension(2:imm,2:jmm)::ap,ae,aw,an,as,b,ap0,apt
    real,dimension(1:im,1:jm)::T,TN,TN_old,k
    real,dimension(2:imm,2:jmm)::Sc,Sp
    real,dimension(2:im,2:jmm)::kx
    real,dimension(2:imm,2:jm)::ky
	real,dimension(1:im)::x,P,Q
	real,dimension(2:im)::xu,dx,xue,xuw
	real,dimension(2:imm)::dxv
	real,dimension(1:jm)::y,R,yr
	real,dimension(2:jm)::yv,dy,Ryv,yvn,yvs
	real,dimension(2:jmm)::dyv
    real::hw,he,Tfw,Tfe
    real::xl,yl,R0
    real::rou,Cp,K0
    real::T0,QB       !T0代表迭代初值,QB热流密度
	integer::mode,last
    real::eps
    real::time,timeend
    real::dt
    integer::tlevel,nprint                       
    logical::Converged_in
    logical::Converged_all
    character(len=20)::filename
end module   
 

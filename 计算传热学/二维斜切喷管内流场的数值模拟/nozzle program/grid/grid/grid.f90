program main
  implicit none
  real(8) pi
  integer i,j,ib,jb,imc,jmc,im,jm,il,in
  real(8) li,lt,lb,lx,lm
  real(8) hi,ht,hb,hx
  real(8) xx,yy,xu,xm,yu,ym
  real(8) alpha,theta,phi

  jmc=80
  imc=150
  il=jmc/2
  pi=atan(1.d0)*4.d0

  alpha=pi/9.d0
  theta=pi/12.d0
!    phi=pi/6.d0
!    phi=pi/4.d0
!    phi=pi/3.d0
!    phi=pi*5.0/12.d0
     phi=pi/2.d0

  hi=1.50
  li=1.0
  
  ht=1.0
  lt=(hi-ht)/tan(alpha)

  hb=2.0
  lb=(hb-ht)/tan(theta)

  lm=li+lt+lb
  
!   lx=lm+4.0*tan(phi)/(1-tan(phi)*tan(theta))
    if (phi==pi/2.d0) then
    lx=lm
	else
    lx=lm+4.0/(tan(phi)-tan(theta))
	end if
!         open (5,file='../grid.in',status='unknown')
        open (5,file='grid.in',status='unknown')
        write(5,*) imc,jmc
        ib=3
        jb=3
        im=imc+ib
        jm=jmc+jb
! inlet rectanglar section x=0+lxin/nx  y=-hi+2*hi/ny
        do i=1,il
		    xx=(i-1)*li/il
        do j=1,jmc+1
		    yy=(j-1)*2*hi/jmc-hi
            write(5,*) xx,yy
        end do
        end do
! convegaration section 
		do i=1,il+1
		    xx=(i-1)*lt/il+li
            hx=hi+(i-1)*(ht-hi)/il
        do j=1,jmc+1
		    yy=(j-1)*2*hx/jmc-hx
            write(5,*) xx,yy
        end do
        end do

		in=imc-2*il	
		do i=2,imc-2*il+1
		    xu=(i-1)*lb/in+li+lt
			xm=(i-1)*(lx-li-lt)/in+li+lt
            yu=ht+(i-1)*(hb-ht)/in
			ym=-ht-(i-1)*(lx-li-lt)*tan(theta)/in
        do j=1,jmc+1
			xx=(j-1)*(xu-xm)/jmc+xm
		    yy=(j-1)*(yu-ym)/jmc+ym
            write(5,*) xx,yy
        end do
        end do
end
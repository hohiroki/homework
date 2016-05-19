subroutine compute()
    use variable
    implicit none
	real::areax,areays,areayn !x,y方向导热面积
    integer::i,j

    do j=2,jmm
        do i=2,imm
		    areax=dyv(j)*R(j)        !x方向导热面积
			areays=dxv(i)*Ryv(j)     !y方向下界面导热面积
			areayn=dxv(i)*Ryv(j+1)   !y方向上界面导热面积
            ae(i,j)=kx(i+1,j)*areax/dx(i+1)
            aw(i,j)=kx(i,j)*areax/dx(i)
            an(i,j)=ky(i,j+1)*areayn/dy(j+1)
            as(i,j)=ky(i,j)*areays/dy(j)
            ap0(i,j)=2.0*rou*Cp*R(j)*dxv(i)*dyv(j)/dt  
            b(i,j)=b(i,j)+ap0(i,j)*TN(i,j) 
            ap(i,j)=aw(i,j)+ae(i,j)+an(i,j)+as(i,j)+ap0(i,j)-Sp(i,j)*R(j)*dxv(i)*dyv(j)
        end do
    end do
      
       return
     end subroutine       

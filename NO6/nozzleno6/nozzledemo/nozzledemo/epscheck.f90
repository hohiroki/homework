subroutine epscheck
!check if eps alright
use main
 eps = 0.0
 eps=abs((rovxml(1,1)-rovx(1,1))/rovx(1,1))
 do i=ib,im
      do j=jb,jm
           eps=eps+ abs((rovxml(i,j)-rovx(i,j))/(rouvx(i,j)+0.00000000001))
      end do
 end do
 
       return
      end
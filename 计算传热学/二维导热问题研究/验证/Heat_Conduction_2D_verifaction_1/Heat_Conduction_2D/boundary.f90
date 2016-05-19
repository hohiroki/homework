subroutine boundary

use variable

implicit none

do j=1,y_node_number
T2(1,j)=T2(2,j)
T2(x_node_number,j)=300
end do

do i=1,x_node_number
T2(i,1)=300
T2(i,y_node_number)=600
end do

return

end subroutine
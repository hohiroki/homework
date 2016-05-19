subroutine boundary

use variable

implicit none

do j=1,y_node_number
    T2(1,j)=0
    T2(x_node_number,j)=0
end do

do i=1,x_node_number
    T2(i,1)=0
    T2(i,y_node_number)=0
end do

return

end subroutine
subroutine boundary

use variable

implicit none

do j=1,y_node_number
    T2(1,j)=exp((0.0+y_node(j))/2.0-time)
    T2(x_node_number,j)=exp((x_length+y_node(j))/2.0-time)
end do

do i=1,x_node_number
    T2(i,1)=exp((0.0+x_node(i))/2.0-time)
    T2(i,y_node_number)=exp((y_length+x_node(i))/2.0-time)
end do

return

end subroutine
subroutine boundary

use variable

implicit none

do j=1,y_node_number
    T2(1,j)=500
end do

do j=2,y_node_number
    k(x_node_number,j)=k_boundary
end do

do j=2,(y_node_number-1)
    S_c(x_node_number-1,j)=100.0+1.0/dx_face(x_node_number-1)*T_f/(1.0/h+dx_face_node_W(x_node_number)/k(x_node_number,j))
    S_p(x_node_number-1,j)=0.0-1.0/dx_face(x_node_number-1)*1.0/(1.0/h+dx_face_node_W(x_node_number)/k(x_node_number,j))
end do

do j=2,y_node_number
    k(x_node_number,j)=0.0
end do

do i=1,x_node_number
    T2(i,1)=300.0
end do

do i=2,(x_node_number-1)
    S_c(i,y_node_number-1)=100+q_boundary/dy_face(y_node_number-1)
end do

do i=2,x_node_number
    k(i,y_node_number)=0
end do

return

end subroutine
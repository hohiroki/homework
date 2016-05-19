subroutine boundary

use variable

implicit none

!East face
do j=2,(y_node_number-1)
    S_c(2,j)=S_c_constant+q_boundary/dx_face(2)
    S_p(2,j)=S_p_constant+0
end do

do j=2,y_node_number
    k(2,j)=0.0
end do

!west face
do j=2,y_node_number
    k(x_node_number,j)=k_boundary
end do

do j=2,(y_node_number-1)
    S_c(x_node_number-1,j)=S_c_constant+1.0/dx_face(x_node_number-1)*T_f/(1.0/h+dx_face_node_W(x_node_number)/k(x_node_number,j))
    S_p(x_node_number-1,j)=S_p_constant-1.0/dx_face(x_node_number-1)*1.0/(1.0/h+dx_face_node_W(x_node_number)/k(x_node_number,j))
end do

do j=2,y_node_number
    k(x_node_number,j)=0.0
end do

!south face
do i=2,x_node_number
    k(i,2)=k_boundary
end do

do i=2,(x_node_number-1)
    S_c(i,2)=S_c_constant+1.0/dy_face(2)*T_f/(1.0/h+dy_face_node_N(2)/k(i,2))
    S_p(i,2)=S_p_constant-1.0/dy_face(2)*1.0/(1.0/h+dy_face_node_N(2)/k(i,2))
end do

do i=2,x_node_number
    k(i,2)=0.0
end do

!north face
do i=2,x_node_number
    k(i,y_node_number)=k_boundary
end do

do i=2,(x_node_number-1)
    S_c(i,y_node_number-1)=S_c_constant+1.0/dy_face(y_node_number-1)*T_f/(1.0/h+dy_face_node_S(y_node_number)/k(i,y_node_number))
    S_p(i,y_node_number-1)=S_p_constant-1.0/dy_face(y_node_number-1)*1.0/(1.0/h+dy_face_node_S(y_node_number)/k(i,y_node_number))
end do

do i=2,x_node_number
    k(i,y_node_number)=0
end do

return

end subroutine
program main

use variable

implicit none

integer:: m=1E+3
coordinate_mode=1


call initial()
call mesh()

convergence_index=.false.

do step=1,m

T1=T2

time=time+dt

call boundary()

call influence_coefficient()              


!=====Alternating direction interation=====!

500 do j=2,(y_node_number-1)
        do i=2,(x_node_number-1)
            T2(i,j)=(a_E(i,j)*T1(i+1,j)+a_W(i,j)*T2(i-1,j)+a_S(i,j)*T2(i,j-1)+a_N(i,j)*T1(i,j+1)+b(i,j))/a_P(i,j)
        end do
    end do

    T1=T2

    do i=2,(x_node_number-1)
        do j=2,(y_node_number-1)
            T2(i,j)=(a_E(i,j)*T1(i+1,j)+a_W(i,j)*T2(i-1,j)+a_S(i,j)*T2(i,j-1)+a_N(i,j)*T1(i,j+1)+b(i,j))/a_P(i,j)
        end do
    end do

    call check_convergence()

    if(.not.convergence_index) then
        T1=T2
        goto 500
    end if


end do

    do j=1,y_node_number
        T2(x_node_number,j)=(h*dx_face_node_W(x_node_number)/k_boundary*T_f+T2(x_node_number-1,j))/(h*dx_face_node_W(x_node_number)/k_boundary+1.0)
    end do

    do i=1,x_node_number
        T2(i,y_node_number)=T2(i,y_node_number-1)+q_boundary*dy_face_node_S(y_node_number)/k_boundary
    end do

call output

end program

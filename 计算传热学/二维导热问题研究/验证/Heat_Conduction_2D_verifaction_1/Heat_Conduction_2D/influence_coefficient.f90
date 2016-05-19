subroutine influence_coefficient

use variable

implicit none

real:: x_area   !heat conduction area in the x direction
real:: y_area_S !heat conduction area in the y direction(the south face)
real:: y_area_N !heat conduction area in the y direction(the north face)

do j=2,(y_node_number-1)
    x_area=dy_face(j)*y_R_node(j)

    do i=2,(x_node_number-1)
        y_area_S=dx_face(i)*y_R_face(j)
        y_area_N=dx_face(i)*y_R_face(j+1)

        a_E(i,j)=k(i+1,j)*x_area/dx_node(i+1)
        a_W(i,j)=k(i,j)*x_area/dx_node(i)
        a_S(i,j)=k(i,j)*y_area_S/dy_node(j)
        a_N(i,j)=k(i,j+1)*y_area_N/dy_node(j+1)

        !a_P0(i,j)=rou*Cp*y_R_node(j)*dx_face(i)*dy_face(j)/dt !nonsteady heat conduction
        a_P0(i,j)=0                                           !steady heat conduction

        a_P(i,j)=a_E(i,j)+a_W(i,j)+a_S(i,j)+a_N(i,j)+a_P0(i,j)-S_p(i,j)*y_R_node(j)*dx_face(i)*dy_face(j)
        b(i,j)=S_c(i,j)*y_R_node(j)*dx_face(i)*dy_face(j)+a_P0(i,j)*T1(i,j)

    end do
end do

return

end subroutine
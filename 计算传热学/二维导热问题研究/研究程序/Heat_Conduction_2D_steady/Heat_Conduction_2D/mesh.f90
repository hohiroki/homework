subroutine mesh

use variable

implicit none

real:: x_grid,y_grid !grid length in the x and y length

x_grid=x_length/(x_node_number-2) !uniform mesh
y_grid=y_length/(y_node_number-2) !uniform mesh

!=====control face position:x_face,y_face=====!
x_face(2)=0
y_face(2)=0

do i=3,x_node_number
    x_face(i)=x_face(i-1)+x_grid
end do

do j=3,y_node_number
    y_face(j)=y_face(j-1)+y_grid
end do

!=====control face distance:dx_face,dy_face=====!
do i=2,(x_node_number-1)
    dx_face(i)=x_face(i+1)-x_face(i)
end do

do j=2,(y_node_number-1)
    dy_face(j)=y_face(j+1)-y_face(j)
end do

!=====node position:x_node,y_node=====!
x_node(1)=x_face(2)
y_node(1)=y_face(2)

do i=2,(x_node_number-1)
    x_node(i)=(x_face(i)+x_face(i+1))/2.0
end do

do j=2,(y_node_number-1)
    y_node(j)=(y_face(j)+y_face(j+1))/2.0
end do

x_node(x_node_number)=x_face(x_node_number)
y_node(y_node_number)=y_face(y_node_number)

!=====node distance:dx_node,dy_node=====!
do i=2,x_node_number
    dx_node(i)=x_node(i)-x_node(i-1)
end do

do j=2,y_node_number
    dy_node(j)=y_node(j)-y_node(j-1)
end do

!=====face node distance:dx_face_node_E,dx_face_node_W,dy_face_node_S,dy_face_node_N=====!
!=====this is for interplot of the heat conductivities=====!
dx_face_node_E(2)=x_node(2)-x_face(2)
dx_face_node_W(2)=0

do i=3,(x_node_number-1)
dx_face_node_E(i)=x_node(i)-x_face(i)
dx_face_node_W(i)=x_face(i)-x_node(i-1)
end do

dx_face_node_E(x_node_number)=0
dx_face_node_W(x_node_number)=x_face(x_node_number)-x_node(x_node_number-1)

dy_face_node_N(2)=y_node(2)-y_face(2)
dy_face_node_S(2)=0

do j=3,(y_node_number-1)
dy_face_node_N(j)=y_face(j)-y_node(j)
dy_face_node_S(j)=y_face(j)-y_node(j-1)
end do

dy_face_node_N(y_node_number)=0
dy_face_node_S(y_node_number)=y_face(y_node_number)-y_node(y_node_number-1)

!=====conversion between the Cartesian coordinate and the Cylindrical coordinate=====!
if (coordinate_mode==1)then
    y_R_face(2:y_node_number)=1
    y_R_node(1:y_node_number)=1
else
    y_R_face(2)=y_R0
    do j=3,y_node_number
        y_R_face(j)=y_R_face(j-1)+dy_face(j-1)
    end do

    y_R_node(1)=y_R0
    do j=2,(y_node_number-1)
        y_R_node(j)=(y_R_face(j)+y_R_face(j+1))/2.0
    end do
    y_R_node(y_node_number)=y_R_face(y_node_number)
end if

return

end subroutine

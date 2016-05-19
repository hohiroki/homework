program main

use variable

implicit none

integer:: m=500
coordinate_mode=1

call initial()
call mesh()

!temperation reintial
do j=1,y_node_number
    do i=1,x_node_number
        T1(i,j)=600.0*sin(x_node(i))*sin(y_node(j))
        T2(i,j)=600.0*sin(x_node(i))*sin(y_node(j))
    end do
end do

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

call output

end program

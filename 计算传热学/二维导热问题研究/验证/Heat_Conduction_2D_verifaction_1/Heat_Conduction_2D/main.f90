program main

use variable

implicit none

coordinate_mode=1

call initial()
call mesh()

convergence_index=.false.

500 call boundary()

    call influence_coefficient()              


!=====Alternating direction interation=====!

    do j=2,(y_node_number-1)
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

call output

end program

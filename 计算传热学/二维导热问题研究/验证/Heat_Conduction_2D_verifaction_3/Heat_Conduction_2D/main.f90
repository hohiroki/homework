program main

use variable

implicit none

integer:: m=1E+3
coordinate_mode=1

call initial()
call mesh()

do j=1,y_node_number
    do i=1,x_node_number
        T1(i,j)=exp((x_node(i)+y_node(j))/2.0)
        T2(i,j)=exp((x_node(i)+y_node(j))/2.0)
    end do
end do

convergence_index=.false.

do step=1,m

T1=T2

time=time+dt

call boundary()

do j=2,(y_node_number-1)
    do i=2,(x_node_number-1)
        S_c(i,j)=-3.0/2.0*exp((x_node(i)+y_node(j))/2.0-time)
        S_p(i,j)=0.0
    end do
end do

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

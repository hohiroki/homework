subroutine initial

use variable

implicit none

real:: temp

x_length=0.2
y_length=0.1

k=130.0
rou=19350.0
Cp=134.0

S_c_constant=0.0
S_p_constant=0.0

open(200,file='T_Wu.dat')

do i=1,3
    read(200,*)
end do

do j=1,101
    do i=1,101
        read(200,*) temp,temp,T1(i,j)
    end do
end do

T2=T1

k_boundary=130.0
T_f=900.0
h=460.0
q_boundary=0.0

dt=1E-2
time=0
error_requirement=1E-5

return

end subroutine
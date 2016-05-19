subroutine initial

use variable

implicit none

x_length=0.2
y_length=0.1

k=50.0
rou=100.0
Cp=50.0

S_c_constant=0.0
S_p_constant=0.0

T1=300.0
T2=300.0

k_boundary=50.0
T_f=900.0
h=460.0
q_boundary=0.0

dt=1E-3
time=0
error_requirement=1E-5

return

end subroutine
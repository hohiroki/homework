subroutine initial

use variable

implicit none

x_length=1.0
y_length=1.0

k=1.0
rou=1.0
Cp=1.0

S_c=100.0
S_p=0

T1=300.0
T2=300.0

k_boundary=1.0
T_f=400.0
h=150.0
q_boundary=200.0

dt=1E-3
time=0
error_requirement=1E-6

return

end subroutine
subroutine check_convergence()

use variable

implicit none

real:: temp

error=0

do j=1,y_node_number
    do i=1,x_node_number
        temp=abs((T2(i,j)-T1(i,j))/(T1(i,j)+1E-30))
        if (temp.gt.error) then
            error=temp
        else
            error=error
        end if
    end do
end do

if (error.gt.error_requirement) then
    convergence_index=.false.
else
    convergence_index=.true.
end if

return

end subroutine
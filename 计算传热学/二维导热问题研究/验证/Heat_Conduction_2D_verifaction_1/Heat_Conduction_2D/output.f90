subroutine output

use variable

implicit none

open(15,file='T.dat')

10 format(2xa)
20 format(1xa22)
30 format(1xa8,i3,6xa3,i3,3xa10)
40 format(2xf6.4,6xf6.4,8xf8.4)

write(15,10) 'Title="2D heat conduction problem"'
write(15,20) 'Variables="x","y","T"'
write(15,30) 'ZONE I=',x_node_number,'J=',y_node_number,'F=POINT'

do j=1,y_node_number
    do i=1,x_node_number
        write(15,40) x_node(i),y_node(j),T2(i,j)
    end do
end do

return

end subroutine
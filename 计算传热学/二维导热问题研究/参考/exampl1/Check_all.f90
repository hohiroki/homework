subroutine Check_all()
    use variable
    implicit none
    integer::i,j
    real::dT_max,temp
    dT_max=0
    do i=1,im
        do j=1,jm
             temp=abs((TN_old(i,j)-T(i,j))/(TN_old(i,j)+1E-20))
             if(temp.gt.dT_max) then 
			     dT_max=temp
			 end if
	     TN_old(i,j)=T(i,j)
         end do
     end do  

     if(dT_max.le.eps)then
         Converged_all=.true.
     else
         Converged_all=.false.
     end if           
       
end subroutine
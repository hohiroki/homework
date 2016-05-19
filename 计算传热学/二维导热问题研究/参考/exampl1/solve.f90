

!==========交替方向线迭代求解================!
subroutine solve()
    use variable
    implicit none
    integer::i,j
	integer*4::step
	real,dimension(1:im,1:jm)::Ttemp
    step=0
    Converged_in=.false.
    do while(.not.Converged_in)
        ! Kernel of the block iteration.
		Converged_in=.false.
		step=step+1
        ! 求解x方向
        call update_Tx()
		! 求解y方向
        call update_Ty()
       ! Check whether converged at this time level or not.
        call check_in()
         do i=1,im
             do j=1,jm
                TN(i,j)=T(i,j)
             end do
         end do   
		if (step==5) then
		   Converged_in=.true.	
		end if
    end do
end subroutine    

subroutine check_in()
          use variable
          implicit none
          integer::i,j
          real::dT_max,temp
          dT_max=abs(TN(1,1)-T(1,1))
          ! Find the max dT
          do j=1,jm
              do i=1,im
                  temp=abs((TN(i,j)-T(i,j))/(T(i,j)+0.1**10))
                  if(temp.gt.dT_max) then 
				    dT_max=temp
				   end if				  
              end do
          end do
          if(abs(dT_max).le.eps)then
              Converged_in=.true.
          else
              Converged_in=.false.
          end if           
end subroutine


subroutine update_Tx()
    use variable
    implicit none
    integer::i,ii,j,jj
	real::temp1,temp2
	!*************x方向往复扫描一次************
     do J=2,jmm
	    P(1)=0
		Q(1)=T(1,j)
	    do i=2,imm
		   temp1=ap(i,j)-P(i-1)*aw(i,j)
		   P(i)=ae(i,j)/temp1
           temp2=b(i,j)+an(i,j)*T(i,j+1)+as(i,j)*T(i,j-1)
           Q(I)=(TEMP2+aw(i,j)*Q(i-1))/temp1
        end do
		do ii=2,imm
		   i=imm+2-ii
		   T(i,j)=T(i+1,j)*P(i)+Q(i)
		end do
	end do

    
	!write(*,*)'update x'
end subroutine

subroutine update_Ty()
    use variable
    implicit none
    integer::i,ii,j,jj
	real::temp1,temp2

	!*************y方向往复扫描一次************
    do i=2,imm
	    P(1)=0
		Q(1)=T(i,1)
	    do j=2,jmm
		   temp1=ap(i,j)-P(j-1)*as(i,j)
		   P(j)=an(i,j)/temp1
           temp2=b(i,j)+ae(i,j)*T(i+1,j)+aw(i,j)*T(i-1,j)
           Q(j)=(TEMP2+as(i,j)*Q(j-1))/temp1
        end do
		do jj=2,jmm
		   j=jmm+2-jj
		   T(i,j)=T(i,j+1)*P(j)+Q(j)
		end do
	end do

   	!write(*,*)'update y'
end subroutine
    
    

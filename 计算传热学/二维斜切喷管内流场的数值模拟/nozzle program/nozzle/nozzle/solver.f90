          subroutine solver(nrk)

!
      use main
!*****************************************************************************
      CALL CPU_TIME ( time1 )
!     compute the area-vectors and interpolation coordinates  
        
!     ==========================      
!      CALL CPU_TIME ( time4 )
!*****************************************************************************
!     set ghost-cell at boundary for various boundary conditions
      call bc_ghostcell_value
!      CALL CPU_TIME ( time5 )
!     ======================

!     interpolate the primitive variables to cell interfaces
      call interpolation_to_interface
!      CALL CPU_TIME ( time6 )
!     ===============================

!     MUSCL interpolation at cell interface except boundaries
      call muscl_interpolation 
!      CALL CPU_TIME ( time10 )
!     ========================

!     set muscl_interpolated value at the cell interfaces on the boundary
      call bc_muscl_interpolation
!      CALL CPU_TIME ( time11 )
!     ===========================

!     evaluate the inviscid fluxes 
!      if(irsolver==0) then
!      call inviscid_fluxes_roe
!      else
!      call inviscid_fluxes_lax
!      end if

	  select case (irsolver)
	  case (0) 
	   call inviscid_fluxes_roe
	  case (1)
	   call inviscid_fluxes_lax
	  case (2)
	   call inviscid_fluxes_rtroe
	  end select
!     ====================  
!      CALL CPU_TIME ( time12 )
!*****************************************************************************
!     Runge-Kutta time stepping
      call Runge_Kutta(nrk)
!     =====================                
!      CALL CPU_TIME ( time13 )
      
!*****************************************************************************
!     update variables
      call update_variables
!     ===================== 
      
      CALL CPU_TIME ( time2 )
      
      write(*,*) 'time spent'
      write(*,*) time2-time1
    !  pause
      return
      end
      
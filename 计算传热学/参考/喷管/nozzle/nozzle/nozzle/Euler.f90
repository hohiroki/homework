!==================================================================!
!======= 2d Euler solver ==========!
!==================================================================!
!                  V 1.0 by REN Yuxin, Oct. 12, 2007
!==================================================================!
!
!
      program Euler
     !======================
     

     use main

!
      open (2,file='Euler.his',status='unknown')
!     bl.his:file that store the computation history
      
!
!     initialization
!     --------------
      write(*,*) '1'
!
      call startup
      call outputuns

      write(*,*) '2'
! 	  write(*,*) maxval(maxval(vy))
! 	  write(*,*) maxval(maxval(vx))
!
!     start of main iteration loop
!     ----------------------------
!
      iend=0
!     if iend==1, the computation stops
      
      istart=0
!     istart:store the advancing time steps for every new start up
!
!     start a new time step
!     ---------------------
 

 100  ip=0
!     if ip==1 output the results

      n=n+1
!     n: accumulated time step of the computation
      istart=istart+1
      
      if(mod(n,1).eq.0) then
          write(*,*) 'now, computing the ',n,'th step'
          write(2,*) 'now, computing the ',n,'th step'
      end if
      if(n.eq.nmax) iend=1
      
!      ip=0.
      if(mod(n,iprint).eq.0) ip=1
      
!     compute time step
      call time_step
!     ============= 
!     
      do nunsloop=1,2
!       nunsloop: number of step in Runge-Kutta time stepping scheme 
!
!     solve Euler equations
!     ------------------
!
      call solver(nunsloop)
!     subroutine solver: solve the ns equation for one Runge-Kutta stage      
      
      end do
!
!
      if(.not.steady1) ttime=ttime+step(ib,jb)
!
!     store variables of previous time steps
!     ---------------------------------------
!

      rms1=0.
      rms2=0.
      rms3=0.
      rms4=0.
      
      
      do i=ib,im-1
      do j=jb,jm-1
      
      tinv=1./step(i,j)
      
      r1=abs(ro  (i,j)-rom1  (i,j))*tinv
      r2=abs(rovx(i,j)-rovxm1(i,j))*tinv
      r3=abs(rovy(i,j)-rovym1(i,j))*tinv
      r4=abs(roe (i,j)-roem1 (i,j))*tinv
      
      if(r1>rms1) then
        rsm1=r1
        imax1=i
        jmax1=j
      end if 
      
      if(r2>rms2) then
        rsm2=r2
        imax2=i
        jmax2=j
      end if 
      
      if(r3>rms3) then
        rsm3=r3
        imax3=i
        jmax3=j
      end if
      
      if(r4>rms4) then
        rsm4=r4
        imax4=i
        jmax4=j
      end if  
      
      end do
      end do
      
      if(mod(n,1).eq.0) then
          write(*,*) 'time', ttime*basel/basev
          write(*,*) rsm1,imax1,jmax1
          write(*,*) rsm2,imax2,jmax2
          write(*,*) rsm3,imax3,jmax3
          write(*,*) rsm4,imax4,jmax4
          write(*,*) rtflk
          write(2,*) rsm1,imax1,jmax1
          write(2,*) rsm2,imax2,jmax2
          write(2,*) rsm3,imax3,jmax3
          write(2,*) rsm4,imax4,jmax4
		  write(2,*) rtflk
      end if
      
      do i=ib,im
      do j=jb,jm 
!
      rom1  (i,j)=ro  (i,j)
      rovxm1(i,j)=rovx(i,j)
      rovym1(i,j)=rovy(i,j)
      roem1 (i,j)=roe (i,j)
!
      end do
      end do

      call output
      if(.not.steady1) call outputuns
      
!     ==============
!
      if(iend.eq.1) then
      close (2)
      write(*,*) 'Normal Termination'
     
      stop
      endif
      goto 100
      end
!
!     end of main loop
!     ----------------
!

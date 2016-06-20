      subroutine area
!     ***************
!     area vectors
!     ***************
      
!
      use main
!
!     this subrouitne computes cell face areas and interpolation coefficients.
!

!=====compute the coordinates of cell centers      
      
 
      do i=ib,im-1
      do j=jb,jm-1
      
          xc(i,j)=0.25*(x    (i,j)+x    (i+1,j)&
     &                      +x  (i,j+1)+x  (i+1,j+1))
          yc(i,j)=0.25*(y    (i,j)+y    (i+1,j)&
     &                      +y  (i,j+1)+y  (i+1,j+1))       
      end do
      end do
     
!
!     i face: surface vector components, grid velocity, shifted volume
!
      do i=ib,im
      do j=jb,jm-1
      	
          aix(i,j)=(y(i,j+1)-y(i,j))  !Ãæ»ýÊ¸Á¿
          aiy(i,j)=-(x(i,j+1)-x(i,j))
          vib(i,j)=0.
          
          aim(i,j)=sqrt(aix(i,j)**2+aiy(i,j)**2)
     
      end do
      end do
!
!     j face: surface vector components, grid velocity, shifted volume
!
      do i=ib,im-1
      do j=jb,jm
          
          ajx(i,j)=-(y(i+1,j)-y(i,j))
          ajy(i,j)=(x(i+1,j)-x(i,j))
          vjb(i,j)=0.
          
          ajm(i,j)=sqrt(ajx(i,j)**2+ajy(i,j)**2)
      end do
      end do
      
      
!     interpolation coefficient

      do i=ib,im-1
      do j=jb,jm-1
      
      xw=0.5*(x(i,j)+x(i,j+1))
      yw=0.5*(y(i,j)+y(i,j+1))
     
      xe=0.5*(x(i+1,j)+x(i+1,j+1))
      ye=0.5*(y(i+1,j)+y(i+1,j+1))
     
      xs=0.5*(x(i,j)+x(i+1,j))
      ys=0.5*(y(i,j)+y(i+1,j))
     
      xn=0.5*(x(i,j+1)+x(i+1,j+1)) 
      yn=0.5*(y(i,j+1)+y(i+1,j+1)) 
     
     
      alagmx(i,j)=sqrt((xw-xe)**2+(yw-ye)**2)
       
      alagmy(i,j)=sqrt((xs-xn)**2+(ys-yn)**2)
       
      
      end do
      end do
     
      return
      end
     
!==============================================================================
      subroutine volume_static
!     **********************************
!     compute volume through coordinates
!     **********************************      
!

      use main

      do i=ib,im-1
      do j=jb,jm-1
      
          vol(i,j)=0.5*((x(i+1,j+1)-x(i,j))*&
     &                    (y(i,j+1)-y(i+1,j))-&
     &                    (y(i+1,j+1)-y(i,j))*&
     &                    (x(i,j+1)-x(i+1,j)))
          
          if(vol(i,j).lt.0.) then
                 write(*,*) 'volume.lt.0',i,j !,k
                 stop
          end if
          
      end do
      end do
!
      return
      end
      

      
      subroutine outputuns
!     ********************
      use main
      character*7 nam
        if(mod(n,iprint)==0) then
        k=(n/iprint)
        if(k<1000) then
         
            if(k<10) then
               write(nam,'(a,i1,a)') '00',k,'.dat'
             else if(k<100) then
               write(nam,'(a,i2,a)') '0',k,'.dat'
             else
               write(nam,'(i3,a)') k,'.dat'
            end if
            
            open(30,file=nam)
            write(30,*) 'title="contour"'
            write(30,*) 'variables="x","y","u","v",'
            write(30,*)'"p","den","mach","entropy","vmu"'
     
            write(30,*) 'zone i=',im-ib,' j=',jm-jb &
                   &,' f=point'
       
          do j=jb,jm-1
          do i=ib,im-1
             xcar=xc(i,j)*basel
             ycar=yc(i,j)*basel
             
             ucar=vx(i,j)*basev
             vcar=vy(i,j)*basev
             
             pcar=p(i,j)*baser*basev*basev
             rcar=ro(i,j)*baser
             entr=pcar/rcar**ga
             amach=sqrt(ucar*ucar+vcar*vcar)&
     &             /sqrt(ga*pcar/rcar)
             vmur=vmu(i,j)*baser*basev*basel
             
             write(30,*) xcar,ycar&
     &                  ,ucar,vcar&
     &                  ,pcar,rcar,amach,entr,vmur
           end do
           end do
           close(30)
        
        
        end if
        end if
      
      return
      end
      
      
      subroutine output
!     *****************
      use main
!
!     save data for restart
!     ---------------------
!
      if(ip.eq.1) then
      write(*,*) 'now,backup the intermediate results'
      open(1,file='bl.sav',form='unformatted')
      write(1) n,ttime
      
      do  i=ib,im
      do  j=jb,jm
          write(1)    vx(i,j),vy(i,j)&
     &               ,p(i,j),ro(i,j)&
     &               ,vmul(i,j),vmu(i,j)&
     &               ,rom1(i,j),roem1(i,j)&
     &               ,rovxm1(i,j),rovym1(i,j)
      end do
      end do
     
      close(1)
      write(*,*) 'backup finish'
      end if
!
      return
      end

      

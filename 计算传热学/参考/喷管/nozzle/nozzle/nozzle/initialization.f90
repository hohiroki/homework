!==============================================================================================      
      subroutine startup
      !=================
      
      use main
!      implicit none

	
      open (4,file='euler.in',status='unknown')
!
!     read in main integer/logical control variables
!     ----------------------------------------------
!
      read(4,*) nmax,irest,iprint 
!         nmax:  maximum advancing steps;
!         irest: start from very beginning(0) or start from a previous computation(1);
!         iprint: output every iprint steps
!

      read(4,*) steady1
!         steady1: logical variable, ==.true. compute steady flow; ==.false. the flow is unsteady

      read(4,*) icgl, irsolver  
!         icgl=1: gloal time step; =0 local time step
!         irsolver=0 Roe; =1 Lax -riemann solver
!     read constants and CFL number 
!     ------------------------------
!
      read(4,*) cp,ga,prl
!         cp: specfic heat at constant pressure
!         ga: ratio of specific heat
!         prl:Prandt.al number
      prte=prl
      prt=prl
      
      read(4,*) cfl
!         cfl: cfl number


!     read inlet/outlet parameters: uniform inlet, need to be modified.
!     -----------------------------------------------------------------
!
!      read(4,*) vxin,vyin,pin,tpin
!      read(4,*) pout
!     
      read(4,*) weno_od
	  call weno2d_init(weno_od)

      close(4)
!===========end of the input file bl.in=========================
!
!     read grids 
!     ----------
!
        open (5,file='grid.in',status='unknown')
        read(5,*) imc,jmc
        ib=3
        jb=3
        im=imc+ib
        jm=jmc+jb

        do i=ib,im
        do j=jb,jm
            read(5,*) x(i,j),y(i,j)
!			print *, i,j,x(i,j),y(i,j)
        end do
        end do
        
      close(5)
!
!     end of data input section
!     -------------------------
!

      open (30,file='grid.dat')
      write(30,*) 'title="contour"'
      write(30,*) 'variables="x","y"'
     
          write(30,*) 'zone i=',im-ib+1,' j=',jm-jb+1&
                    &,' f=point'
       
          do j=jb,jm
          do i=ib,im
             xcar=x(i,j)
             ycar=y(i,j) 
             write(30,*) xcar,ycar
          end do
          end do
        
       close(30)
            
!============================================================================      
!
!      non-dimensionalization
!
      
       basel=1. !x(im,jb)-x(ib,jb)
       baser=1. !pin/(tpin*(ga-1.)*cp)
       basev=1. !sqrt(pin/baser)
       baset=1. !tpin 
       basec=1. !sqrt((ga-1.)*cp*baset)
       basem=1. !sqrt(vxin*vxin+vyin*vyin)/basec
       
       
       write(*,*) 'baseparameters', basel,baser,basev,baset,basem
       
!      non-dim parameters

       cp=cp*baset/(basev*basev)
       
!      non-dim inlet/outlet parameters 
       
       vxin=vxin/basev
       vyin=vyin/basev
       
       pin=pin/(baser*basev*basev)
       tpin=tpin/baset
       pout=pout/(baser*basev*basev)
       
!      non-dim coordinates       
       
        do i=ib,im
        do j=jb,jm
           x(i,j)=x(i,j)/basel
           y(i,j)=y(i,j)/basel
        end do
        end do



! geometry
  li=1.0
  
  lt=(hi-ht)/tan(alpha)

  lb=(hb-ht)/tan(theta)

  lm=li+lt+lb
  lx=lm+4*tan(phi)/(1-tan(phi)*tan(theta))
!
!
!     set various constants
!     ---------------------
!
      ga1=ga-1
      fga=(ga-1.0)/ga
      rfga=1.0/fga
      rcp=1.0/cp
      cv=cp/ga
      rcpcv=cp-cv

!
!     compute aera vectors
!     --------------------------------------
!

      call area
      call volume_static    
       
!     ==========
!
!
!
!     initialization
!     --------------
!
       ror0=1.4
       ur0=0.
       vr0=0.
       pr0=10.0

       rol0=5.625
       ul0=0.0
       vl0=0.0
       pl0=60

	  p0=pl0
	  rho0=rol0
      tp0=p0/rho0/rcpcv

      if(irest.eq.0) n=0
!
!     initial condition
!     -----------------
!
!
!       do i=ib,im-1
!       do j=jb,jm-1
!               xx=0.25*(x(i,j)+x(i,j+1)+x(i+1,j)+x(i+1,j+1))
!               yy=0.25*(y(i,j)+y(i,j+1)+y(i+1,j)+y(i+1,j+1))
!               flag=1.73205*(xx-0.1666667)-yy
!               if(flag.gt.0) then
!                    ro(i,j)=1.4
!                    vx(i,j)=0.
!                    vy(i,j)=0.
!                     p(i,j)=1.0
!               else
!                     ro(i,j)=8.0
!                     vx(i,j)=7.1447
!                     vy(i,j)=-4.125
!                      p(i,j)=116.5
!               endif
!       end do
!       end do
      il=jmc/20
      do i=ib,il+ib-1
      do j=jb,jm-1

                   ro(i,j)=rol0
                   vx(i,j)=ul0
                   vy(i,j)=vl0
                    p(i,j)=pl0
! 		   vx(i,j)=1.578172536102060
!            vy(i,j)=0.0
! 		   e      =.5*(vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j))
! 		   tp(i,j)=tp0-e/cp
! 		   p (i,j)=p0*(tp(i,j)/tp0)**rfga
! 		   ro(i,j)=p(i,j)/tp(i,j)/rcpcv

      end do
      end do

      do i=il+ib,im-1
      do j=jb,jm-1

                   ro(i,j)=ror0
                   vx(i,j)=ur0
                   vy(i,j)=vr0
                    p(i,j)=pr0

      end do
      end do
      
      do  i=ib,im-1
      do  j=jb,jm-1
      
      e=.5*(vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j))
      roe(i,j)=p(i,j)/(ga-1)+(e)*ro(i,j)
      ho(i,j)=ga*(roe(i,j)/ro(i,j)-e)+e
      ros=ro(i,j)
      rovx(i,j)=ros*vx(i,j)
      rovy(i,j)=ros*vy(i,j)
      tp(i,j)=p(i,j)/(ro(i,j)*rcpcv)

      rom1  (i,j)=ro  (i,j)
      roem1 (i,j)=roe (i,j)
      rovxm1(i,j)=rovx(i,j)
      rovym1(i,j)=rovy(i,j)
      
      end do
      end do
           
     
!
!     restart from backup intermediate results
!     --------------------------
!
      if(irest.eq.0) goto 10000
!
!      read:
!      velocities, static pressure, density,viscosity
!      ----------------------------------------------
!
      write(*,*) 'read  backuped intermediate results'
      open(5,file='bl.rd', form='unformatted')  !
      read(5) n,ttime
      
      do  i=ib,im
      do  j=jb,jm
             read(5) vx(i,j),vy(i,j)&
     &               ,p(i,j),ro(i,j)&
     &               ,vmul(i,j),vmu(i,j)&
     &               ,rom1(i,j),roem1(i,j)&
     &               ,rovxm1(i,j),rovym1(i,j)
     
      end do
      end do
!
      close(5)
      
      do i=ib,im-1
      do j=jb,jm-1
      temini=tpin*baset
      vmul(i,j)=1.458*abs(temini)**1.5/(temini+110.4)*1.0d-6/&
         (baser*basev*basel)
      vmu(i,j)=vmul(i,j)
      end do
      end do
      
      write(*,*) 'read  backuped intermediate results finish'
!
!      other variables
!      ---------------
!
      do  i=ib,im-1
      do  j=jb,jm-1
      
      vsqh=.5*(vx(i,j)*vx(i,j)+vy(i,j)*vy(i,j))
      roe(i,j)=p(i,j)/ga1+vsqh*ro(i,j)
      ho(i,j)=ga*(roe(i,j)/ro(i,j)-vsqh)+vsqh
      ros=ro(i,j)
      rovx(i,j)=ros*vx(i,j)
      rovy(i,j)=ros*vy(i,j)
      tp(i,j)=p(i,j)/(ro(i,j)*rcpcv)
      
      end do
      end do

10000 continue
!
      return
      end

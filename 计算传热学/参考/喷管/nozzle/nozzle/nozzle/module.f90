      module main
      
      parameter(iq=109,jq=59,ngmax=max(iq,jq)+9)
      parameter(small=1.0d-20,small1=1.0d-10,alarge=1.0d20)
      parameter(one=1.0,zero=0.0)  
	  parameter( pi=3.1415926535897932384626)
	  parameter(hi=1.5,ht=1.0,hb=2.0)
	  parameter(alpha=pi/9.0,theta=pi/12.0,phi=pi/6.0)
      
      logical :: steady1
      
      
      integer :: &
     & ib,jb &
     &,im,jm,irsolver
     
      real :: &
     & basel,basev,baser,baset,basem

	 real :: rol0,pl0,ul0,vl0&
	 &      ,ror0,pr0,ur0,vr0
      
     real :: li,lt,lb,lm,lx

     real :: &
     & ro(iq,jq),p(iq,jq)&
     &,roe(iq,jq),tp(iq,jq)&
     &,rovx(iq,jq),rovy(iq,jq)&
     &,vx(iq,jq),vy(iq,jq)&
     &,ho(iq,jq)&
     &,hro(iq,jq),hre(iq,jq)&
     &,hrx(iq,jq),hry(iq,jq)&
     &,rom1(iq,jq),roem1(iq,jq)&
     &,rovxm1(iq,jq),rovym1(iq,jq)
     
     
      real :: &
     & x(iq,jq),y(iq,jq)&
     &,xc(iq,jq),yc(iq,jq),step(iq,jq)&
     &,vol(iq,jq)

!     
      real :: &
     &            aix(iq,jq),aiy(iq,jq)&
     &           ,ajx(iq,jq),ajy(iq,jq)&
     &           ,aim(iq,jq),ajm(iq,jq)&
     &           ,alagmx(iq,jq),alagmy(iq,jq)
!   
      real :: &
!
     & vib(iq,jq),vjb(iq,jq)
!       
      real :: &
!      
     &uil (iq,jq),uir(iq,jq),&
     &vil (iq,jq),vir(iq,jq),&
     &pil (iq,jq),pir(iq,jq),&
     &ril (iq,jq),rir(iq,jq),&

     &ujl (iq,jq),ujr(iq,jq),&
     &vjl (iq,jq),vjr(iq,jq),&
     &pjl (iq,jq),pjr(iq,jq),&
     &rjl (iq,jq),rjr(iq,jq)

     
    
      integer :: &
     &nmax,n,irest,&
     &ioup,iinp,icgl,iimplicit,idgstart&
     &,nts1p,nunsloop,ninnerloop,ip,iprint,istart,nb

      real :: &
!
     & cfl,ttime,tstep,epsall
!     
      real :: &
!
     & cp,ga,fga,ga1,rcp,rcpcv,rfga,cv,prl,prt,prte
!     
!
      real :: &
!
     &vmul(iq,jq),vmu(iq,jq)


	  real :: &
    & vxin,vyin,tpin,pin,pout
    
     real :: &
     &  rms1,rms2,rms3,rms4
     
     integer :: imax1,jmax1, &
    &            imax2,jmax2, &
	  &            imax3,jmax3, &
    &            imax4,jmax4

     real :: du2,adu,rtnc(3,2)
     integer :: rtm,rtfli,rtflj,rtflk

!   weno addin
     integer :: weno_od
	 integer,parameter :: weno_max=5
	 real(8) :: Binv(weno_max,weno_max,weno_max)!,Pint(weno_max,weno_max),&
!    &           BPBt(weno_max,weno_max)
     real(8) :: crb(2,weno_max,weno_max),&
!     
	&             crg(2,weno_max,weno_max),&
!
	&              drb(2,weno_max),  &
	&              drg(2,weno_max),  &
!
	&             cowr(weno_max,weno_max)

!     real    :: stencil(weno_max),pol(weno_max),kxi

	 real    :: rownj(2,iq,jq),rowni(2,iq,jq),&
    &           vxwnj(2,iq,jq),vxwni(2,iq,jq),&
    &           vywnj(2,iq,jq),vywni(2,iq,jq),&
    &            pwnj(2,iq,jq), pwni(2,iq,jq),&
!
!
	&            rilg(2,iq,jq),  rirg(2,iq,jq),&
	&            uilg(2,iq,jq),  uirg(2,iq,jq),&
	&            vilg(2,iq,jq),  virg(2,iq,jq),&
	&            pilg(2,iq,jq),  pirg(2,iq,jq),&
!
!
	&            rjlg(2,iq,jq), rjrg(2,iq,jq),&
	&            ujlg(2,iq,jq), ujrg(2,iq,jq),&
	&            vjlg(2,iq,jq), vjrg(2,iq,jq),&
	&            pjlg(2,iq,jq), pjrg(2,iq,jq)



     end module main
      
!

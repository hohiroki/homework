      SUBROUTINE CALCU
	INCLUDE 'HEAD'
C-----------Interior Points
	DO 10 J=2,NJM1
	DO 10 I=3,NIM1
C-----------Calculate Areas and Volume	
	AREAS=SWEU(I)*RV(J)
	AREAN=SWEU(I)*RV(J+1)
	AREAW=SSN(J)*RCV(J)
	AREAE=SSN(J)*RCV(J)
	VOL=SWEU(I)*SSN(J)*RCV(J)
C-----------Calculate Convection Coefficients
      CSW=(DEN(I-1,J-1)*XICVS(J)+DEN(I-1,J)*XICVN(J))*V(I-1,J)
	CSE=(DEN(I,J-1)*XICVS(J)+DEN(I,J)*XICVN(J))*V(I,J)
      CS=(CSW*XICUW(I)+CSE*XICUE(I))*AREAS
	CNW=(DEN(I-1,J)*XICVS(J+1)+DEN(I-1,J+1)*XICVN(J+1))
     1	*V(I-1,J+1)
	CNE=(DEN(I,J)*XICVS(J+1)+DEN(I,J+1)*XICVN(J+1))*V(I,J+1)
	CN=(CNW*XICUW(I)+CNE*XICUE(I))*AREAN
	CW=DEN(I-1,J)*(U(I-1,J)*XICCW(I-1)+U(I,J)*XICCE(I-1))*AREAW
	IF(I.EQ.3)CW=DEN(I-2,J)*U(I-1,J)*AREAW
	CE=DEN(I,J)*(U(I,J)*XICCW(I)+U(I+1,J)*XICCE(I))*AREAE
	IF(I.EQ.NIM1)CE=DEN(I+1,J)*U(I+1,J)*AREAE
	REMASS=CN-CS+CE-CW
C---------Calculate Diffusion Coefficients
      VISSW=(VIS(I-1,J-1)*XICVS(J)+VIS(I-1,J)*XICVN(J))
      VISSE=(VIS(I,J-1)*XICVS(J)+VIS(I,J)*XICVN(J))
      VISS=(VISSW*XICUW(I)+VISSE*XICUE(I))
      VISNW=(VIS(I-1,J)*XICVS(J+1)+VIS(I-1,J+1)*XICVN(J+1))
      VISNE=(VIS(I,J)*XICVS(J+1)+VIS(I,J+1)*XICVN(J+1))
	VISN=(VISNW*XICUW(I)+VISNE*XICUE(I))
	VISW=VIS(I-1,J)
	IF(I.EQ.3)VISW=VIS(I-2,J)
	VISE=VIS(I,J)
	IF(I.EQ.NIM1)VISE=VIS(I+1,J)	
	DS=VISS*AREAS/DYSP(J)+1E-30
	DN=VISN*AREAN/DYPN(J)+1E-30
	DW=VISW*AREAW/DXWPU(I)+1E-30
	DE=VISE*AREAE/DXPEU(I)+1E-30
C---------Calculate the Peclet Number	
      PS=CS/DS
	PN=CN/DN
	PW=CW/DW
	PE=CE/DE
C---------Assemble Main Coefficients
      AS(I,J)=DS*AMAX1(0.,(1.-.1*ABS(PS))**5)+AMAX1(CS,.0)
	AN(I,J)=DN*AMAX1(0.,(1.-.1*ABS(PN))**5)+AMAX1(-CN,.0)
	AW(I,J)=DW*AMAX1(0.,(1.-.1*ABS(PW))**5)+AMAX1(CW,.0)
	AE(I,J)=DE*AMAX1(0.,(1.-.1*ABS(PE))**5)+AMAX1(-CE,.0)
c     AS(I,J)=AMAX1(ABS(0.5*CS),DS)+0.5*CS
c	AN(I,J)=AMAX1(ABS(0.5*CN),DN)-0.5*CN
c	AW(I,J)=AMAX1(ABS(0.5*CW),DW)+0.5*CW
c	AE(I,J)=AMAX1(ABS(0.5*CE),DE)-0.5*CE
      IF(REMASS.LE.0.0)THEN
	SP(I,J)=0.0
	SU(I,J)=-REMASS*U(I,J)
	ELSE
	SP(I,J)=-REMASS
	SU(I,J)=0.0
	END IF
C---------Calculate the Source
	DUDXE=VISE*(U(I+1,J)-U(I,J))/DXPEU(I)*AREAE
	DUDXW=VISW*(U(I,J)-U(I-1,J))/DXWPU(I)*AREAW
	SU(I,J)=SU(I,J)+(DUDXE-DUDXW)
	DVDXN=VISN*(V(I,J+1)-V(I-1,J+1))/SWEU(I)*AREAN
	DVDXS=VISS*(V(I,J)-V(I-1,J))/SWEU(I)*AREAS
	SU(I,J)=SU(I,J)+(DVDXN-DVDXS)
C---------Callculate the Source of Pressure
      DDU(I,J)=(0.5*AREAW+0.5*AREAE)
     	SU(I,J)=SU(I,J)+DDU(I,J)*(P(I-1,J)-P(I,J))
10    CONTINUE
C---------Problem Modification
      CALL PRMOD(1)
C---------Final AP(I,J,K) and Residual Source Calculation
	RESORU=0
	DO 20 J=2,NJM1
	DO 20 I=3,NIM1
	AP(I,J)=AN(I,J)+AS(I,J)+AW(I,J)+AE(I,J)-SP(I,J)
	RESOR=AS(I,J)*U(I,J-1)+AN(I,J)*U(I,J+1)+
     1      AW(I,J)*U(I-1,J)+AE(I,J)*U(I+1,J)+
     1      SU(I,J)-AP(I,J)*U(I,J)
      RESORU=RESORU+ABS(RESOR)
20    CONTINUE
C---------Under-Relaxation
	DO 30 J=2,NJM1
	DO 30 I=3,NIM1
	AP(I,J)=AP(I,J)/URFU
      SU(I,J)=SU(I,J)+(1.-URFU)*AP(I,J)*U(I,J)
30    CONTINUE
C---------SIMPLEC Method
	DO 40 J=2,NJM1
	DO 40 I=3,NIM1
C          TEMP=AP(I,J,K)
C          DDU(I,J,K)=DDU(I,J,K)/TEMP	
	TEMP=AP(I,J)-AN(I,J)-AS(I,J)-AW(I,J)-AE(I,J)+1E-30
	DDU(I,J)=DDU(I,J)/TEMP
40    CONTINUE
C---------Solution of Difference Equation
      DO 50 N=1,NSWPU
	CALL LISOLV(3,NIM1,2,NJM1,IT,JT,U,1)
50    CONTINUE
      RETURN
	END
      SUBROUTINE CALCED
      INCLUDE 'HEAD'
C-----------Interior points	
	DO 10 J=2,NJM1
	DO 10 I=2,NIM1
C-----------Calculate Areas and Volume
      AREAS=SWE(I)*RV(J)
	AREAN=SWE(I)*RV(J+1)
	AREAW=SSN(J)*RCV(J)
	AREAE=SSN(J)*RCV(J)
	VOL=SWE(I)*SSN(J)*RCV(J)
C---------Calculate Convection Coefficients
      CS=(DEN(I,J-1)*XICVS(J)+DEN(I,J)*XICVN(J))*V(I,J)*AREAS
	CN=(DEN(I,J)*XICVS(J+1)+DEN(I,J+1)*XICVN(J+1))*V(I,J+1)*AREAN
	CW=(DEN(I-1,J)*XICUW(I)+DEN(I,J)*XICUE(I))*U(I,J)*AREAW
	CE=(DEN(I,J)*XICUW(I+1)+DEN(I+1,J)*XICUE(I+1))*U(I+1,J)*AREAE
	REMASS=CN-CS+CE-CW
C---------Calculate Diffusion Coefficients
      VISS=(VIS(I,J-1)*XICVS(J)+VIS(I,J)*XICVN(J))/PRED
	VISN=(VIS(I,J)*XICVS(J+1)+VIS(I,J+1)*XICVN(J+1))/PRED
	VISW=(VIS(I-1,J)*XICUW(I)+VIS(I,J)*XICUE(I))/PRED
	VISE=(VIS(I,J)*XICUW(I+1)+VIS(I+1,J)*XICUE(I+1))/PRED
	DS=VISS*AREAS/DYSP(J)+1E-30
	DN=VISN*AREAN/DYPN(J)+1E-30
	DW=VISW*AREAW/DXWP(I)+1E-30
	DE=VISE*AREAE/DXPE(I)+1E-30
C---------Calculate the Peclet Number
      PS=CS/DS
	PN=CN/DN
	PW=CW/DW
	PE=CE/DE
C---------Assemble Main Cofficients
      AS(I,J)=DS*AMAX1(0.,(1.-.1*ABS(PS))**5)+AMAX1(CS,0.)
	AN(I,J)=DN*AMAX1(0.,(1.-.1*ABS(PN))**5)+AMAX1(-CN,.0)
	AW(I,J)=DW*AMAX1(0.,(1.-.1*ABS(PW))**5)+AMAX1(CW,0.)
	AE(I,J)=DE*AMAX1(0.,(1.-.1*ABS(PE))**5)+AMAX1(-CE,.0)
c     AS(I,J)=AMAX1(ABS(0.5*CS),DS)+0.5*CS
c	AN(I,J)=AMAX1(ABS(0.5*CN),DN)-0.5*CN
c	AW(I,J)=AMAX1(ABS(0.5*CW),DW)+0.5*CW
c	AE(I,J)=AMAX1(ABS(0.5*CE),DE)-0.5*CE
      IF(REMASS.LE.0.0)THEN
	SP(I,J)=0.0
	SU(I,J)=-REMASS*ED(I,J)
	ELSE
	SP(I,J)=-REMASS
	SU(I,J)=0.0
	END IF
C---------Calculate the Source
      SP(I,J)=SP(I,J)-C2*DEN(I,J)*ED(I,J)/(TE(I,J)+1E-30)*VOL
      SU(I,J)=SU(I,J)+C1*CMU*DEN(I,J)*TE(I,J)*GENT(I,J)*VOL
10    CONTINUE
C---------Problem Modification
      CALL PRMOD(5)
C---------Final AP(I,J,K) and Residual Source Calculation
      RESORE=0.0
      DO 20 J=2,NJM1
	DO 20 I=2,NIM1
	AP(I,J)=AN(I,J)+AS(I,J)+AW(I,J)+AE(I,J)-SP(I,J)	
	RESOR=AS(I,J)*ED(I,J-1)+AN(I,J)*ED(I,J+1)+
     1      AW(I,J)*ED(I-1,J)+AE(I,J)*ED(I+1,J)+
     1      SU(I,J)-AP(I,J)*ED(I,J)
      IF(ABS(SP(I,J)).GE.1E25)RESOR=0.0     
      RESORE=RESORE+ABS(RESOR)                	
20    CONTINUE
C---------Under-Relaxation
	DO 30 J=2,NJM1
	DO 30 I=2,NIM1
	AP(I,J)=AP(I,J)/URFE
      SU(I,J)=SU(I,J)+(1.-URFE)*AP(I,J)*ED(I,J)
30    CONTINUE
C---------Solution of Difference Equation
      DO 40 N=1,NSWPD
	CALL LISOLV(2,NIM1,2,NJM1,IT,JT,ED,1)
40    CONTINUE
      DO 1110 I=1,NI
	ED(I,1)=ED(I,2)
1110	CONTINUE
      DO 1120 J=1,NJ
	ED(NI,J)=ED(NIM1,J)
1120  CONTINUE	
C************30/10
      RETURN
	END
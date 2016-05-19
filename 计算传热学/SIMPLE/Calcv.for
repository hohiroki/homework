      SUBROUTINE CALCV
	INCLUDE 'HEAD'
C-----------Interior Points
	DO 10 J=3,NJM1
	DO 10 I=2,NIM1
C-----------Calculate Areas and Volume
      AREAS=SWE(I)*R(J-1)
	AREAN=SWE(I)*R(J)
	AREAW=SSNV(J)*RV(J)
	AREAE=SSNV(J)*RV(J)
	VOL=SWE(I)*SSNV(J)*RV(J)
C---------Calculate Convection Coefficients
      CS=DEN(I,J-1)*(V(I,J)*XICCN(J-1)+V(I,J-1)*XICCS(J-1))*AREAS
	IF(J.EQ.3)CS=DEN(I,J-2)*V(I,J-1)*AREAS
	CN=DEN(I,J)*(V(I,J+1)*XICCN(J)+V(I,J)*XICCS(J))*AREAN
	IF(J.EQ.NJM1)CN=DEN(I,J+1)*V(I,J+1)*AREAN
	CWN=(DEN(I-1,J)*XICUW(I)+DEN(I,J)*XICUE(I))*U(I,J)
	CWS=(DEN(I-1,J-1)*XICUW(I)+DEN(I,J-1)*XICUE(I))*U(I,J-1)
	CW=(CWN*XICVN(J)+CWS*XICVS(J))*AREAW
	CEN=(DEN(I,J)*XICUW(I+1)+DEN(I+1,J)*XICUE(I+1))*U(I+1,J)
	CES=(DEN(I,J-1)*XICUW(I+1)+DEN(I+1,J-1)*XICUE(I+1))
     1	*U(I+1,J-1)
	CE=(CEN*XICVN(J)+CES*XICVS(J))*AREAE
	REMASS=CN-CS+CE-CW
C---------Calculate Diffusion Coefficients
	VISN=VIS(I,J)
	IF(J.EQ.NJM1)VISN=VIS(I,J+1)
	VISS=VIS(I,J-1)
	IF(J.EQ.3)VISS=VIS(I,J-2)		
	VISWN=(VIS(I-1,J)*XICUW(I)+VIS(I,J)*XICUE(I))
	VISWS=(VIS(I-1,J-1)*XICUW(I)+VIS(I,J-1)*XICUE(I))
	VISW=(VISWN*XICVN(J)+VISWS*XICVS(J))
	VISEN=(VIS(I,J)*XICUW(I+1)+VIS(I+1,J)*XICUE(I+1))
	VISES=(VIS(I,J-1)*XICUW(I+1)+VIS(I+1,J-1)*XICUE(I+1))
	VISE=(VISEN*XICVN(J)+VISES*XICVS(J))
      DS=VISS*AREAS/DYSPV(J)+1E-30
	DN=VISN*AREAN/DYPNV(J)+1E-30
	DW=VISW*AREAW/DXWP(I)+1E-30
	DE=VISE*AREAE/DXPE(I)+1E-30
C---------Calculate the Peclet Number
      PS=CS/DS
	PN=CN/DN
	PW=CW/DW
	PE=CE/DE
C---------Assemble Main Coefficients
      AS(I,J)=DS*AMAX1(0.,(1.-.1*ABS(PS))**5)+AMAX1(CS,0.)
	AN(I,J)=DN*AMAX1(0.,(1.-.1*ABS(PN))**5)+AMAX1(-CN,0.)
	AW(I,J)=DW*AMAX1(0.,(1.-.1*ABS(PW))**5)+AMAX1(CW,0.)
	AE(I,J)=DE*AMAX1(0.,(1.-.1*ABS(PE))**5)+AMAX1(-CE,0.)
c     AS(I,J)=AMAX1(ABS(0.5*CS),DS)+0.5*CS
c	AN(I,J)=AMAX1(ABS(0.5*CN),DN)-0.5*CN
c	AW(I,J)=AMAX1(ABS(0.5*CW),DW)+0.5*CW
c	AE(I,J)=AMAX1(ABS(0.5*CE),DE)-0.5*CE
      IF(REMASS.LE.0.0)THEN
	SP(I,J)=0.0
	SU(I,J)=-REMASS*V(I,J)
	ELSE
	SP(I,J)=-REMASS
	SU(I,J)=0.0
	END IF
C---------Calculate the Source
      DUDYW=VISW*(U(I,J)-U(I,J-1))/SSNV(J)*AREAW
	DUDYE=VISE*(U(I+1,J)-U(I+1,J-1))/SSNV(J)*AREAE
	SU(I,J)=SU(I,J)+(DUDYE-DUDYW)
	DVDYS=VISS*(V(I,J)-V(I,J-1))/DYSPV(J)*AREAS
	DVDYN=VISN*(V(I,J+1)-V(I,J))/DYPNV(J)*AREAN
	SU(I,J)=SU(I,J)+(DVDYN-DVDYS)
	IF(INDCOS.EQ.1)THEN
c	DWDYL=VISL*(W(I,J,K)-W(I,J-1,K))/SSNV(J)*AREAL
c      DWDYR=VISR*(W(I,J,K+1)-W(I,J-1,K+1))/SSNV(J)*AREAR	
      ELSE	
c	DWDYL=VISL*(W(I,J,K)-W(I,J-1,K))/SSNV(J)*AREAL-
c    1      VISL*(W(I,J,K)*XICVN(J)+W(I,J-1,K)*XICVS(J))/YV(J)*AREAL
c      DWDYR=VISR*(W(I,J,K+1)-W(I,J-1,K+1))/SSNV(J)*AREAR-
c     1      VISR*(W(I,J,K+1)*XICVN(J)+W(I,J-1,K+1)*XICVS(J))/YV(J)
c     1      *AREAR     
      END IF
c	SU(I,J,K)=SU(I,J,K)+(DWDYR-DWDYL)
	IF(INDCOS.NE.1)THEN
c	DWDZ=(VISS*XICVS(J)+VISN*XICVN(J))*
c     1	(W(I,J,K+1)*XICVN(J)+W(I,J-1,K+1)*XICVS(J)-
c     1       W(I,J,K)*XICVN(J)-W(I,J-1,K)*XICVS(J))/YV(J)*
c     1      (0.5*AREAL+0.5*AREAR)
c      SU(I,J,K)=SU(I,J,K)-2*DWDZ
c      SU(I,J,K)=SU(I,J,K)+(DEN(I,J,K)*XICVN(J)+DEN(I,J-1,K)*XICVS(J))*
c     1         ((W(I,J,K)*XICVN(J)+W(I,J-1,K)*XICVS(J))*XICCL(K)+
c     1          (W(I,J,K+1)*XICVN(J)+W(I,J-1,K+1)*XICVS(J))*XICCR(K))
c     1         **2/YV(J)*VOL
      SP(I,J)=SP(I,J)-2*(VISS*XICVS(J)+VISN*XICVN(J))
     1          /YV(J)/YV(J)*VOL     
      END IF
C---------Callculate the Source of Pressure
      DDV(I,J)=(0.5*AREAS+0.5*AREAN)
     	SU(I,J)=SU(I,J)+DDV(I,J)*(P(I,J-1)-P(I,J))     
10    CONTINUE
C---------Problem Modification
      CALL PRMOD(2)
C---------Final AP(I,J,K) and Residual Source Calculation
	RESORV=0
	DO 20 J=3,NJM1
	DO 20 I=2,NIM1
	AP(I,J)=AN(I,J)+AS(I,J)+AW(I,J)+AE(I,J)-SP(I,J)	
	RESOR=AS(I,J)*V(I,J-1)+AN(I,J)*V(I,J+1)+
     1      AW(I,J)*V(I-1,J)+AE(I,J)*V(I+1,J)+
     1      SU(I,J)-AP(I,J)*V(I,J)
      RESORV=RESORV+ABS(RESOR)                  	
20    CONTINUE
C---------Under-Relaxation
	DO 30 J=3,NJM1
	DO 30 I=2,NIM1
	AP(I,J)=AP(I,J)/URFV
      SU(I,J)=SU(I,J)+(1.-URFV)*AP(I,J)*V(I,J)
30    CONTINUE
C---------SIMPLEC
	DO 40 J=3,NJM1
	DO 40 I=2,NIM1	
CCCCCCCC      TEMP=AP(I,J)
CCCCCCCC      DDV(I,J)=DDV(I,J)/TEMP	
	TEMP=AP(I,J)-AN(I,J)-AS(I,J)-AW(I,J)-AE(I,J)+1e-30
      DDV(I,J)=DDV(I,J)/TEMP
40    CONTINUE
C---------Solution of Difference Equation
      DO 50 N=1,NSWPV
	CALL LISOLV(2,NIM1,3,NJM1,IT,JT,V,1)
50    CONTINUE
      RETURN
	END
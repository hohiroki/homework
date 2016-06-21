      SUBROUTINE CALCP
      INCLUDE 'HEAD'
	DO 710 I=1,NI
	DO 710 J=1,NJ
	PP(I,J)=0.0
710   CONTINUE		
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
	CN=(DEN(I,J)*XICVS(J+1)+DEN(I,J+1)*XICVN(J+1))
     1	*V(I,J+1)*AREAN
	CW=(DEN(I-1,J)*XICUW(I)+DEN(I,J)*XICUE(I))*U(I,J)*AREAW
	CE=(DEN(I,J)*XICUW(I+1)+DEN(I+1,J)*XICUE(I+1))
     1	*U(I+1,J)*AREAE
	REMASS=CN-CS+CE-CW
	SP(I,J)=0.0
      SU(I,J)=-REMASS
C---------ASSEMBLE MAIN COEFFICIENTS
      AS(I,J)=(DEN(I,J-1)*XICVS(J)+DEN(I,J)*XICVN(J))*
     1	    DDV(I,J)*AREAS
      AN(I,J)=(DEN(I,J)*XICVS(J+1)+DEN(I,J+1)*XICVN(J+1))*
     1	    DDV(I,J+1)*AREAN
      AW(I,J)=(DEN(I-1,J)*XICUW(I)+DEN(I,J)*XICUE(I))*
     1          DDU(I,J)*AREAW     
      AE(I,J)=(DEN(I,J)*XICUW(I+1)+DEN(I+1,J)*XICUE(I+1))*
     1	    DDU(I+1,J)*AREAE
10    CONTINUE
C---------Problem Modification
      CALL PRMOD(3)
C---------Final AP(I,J,K) and Residual Source Calculation
      RESORM=0.0
      DO 20 J=2,NJM1
	DO 20 I=2,NIM1
	AP(I,J)=AS(I,J)+AN(I,J)+AW(I,J)+AE(I,J)-SP(I,J)	
	RESOR=SU(I,J)
      RESORM=RESORM+ABS(RESOR)                	
20    CONTINUE
C---------Solution of Difference Equation
      DO 30 N=1,NSWPP
	CALL LISOLV(2,NIM1,2,NJM1,IT,JT,PP,1)
30    CONTINUE
C---------Velocities U(I,J,K): Modified by pp
	DO 40 J=2,NJM1
	DO 40 I=3,NIM1
	U(I,J)=U(I,J)+DDU(I,J)*(PP(I-1,J)-PP(I,J))
40    CONTINUE	
C---------Velocities V(I,J,K): Modified by pp
	DO 50 J=3,NJM1
	DO 50 I=2,NIM1
	V(I,J)=V(I,J)+DDV(I,J)*(PP(I,J-1)-PP(I,J))		
50    CONTINUE	
C---------Pressures(With Provision for Under-Relaxation)
      PPREF=PP(IPREF,JPREF)
	DO 70 J=2,NJM1
	DO 70 I=2,NIM1
	P(I,J)=P(I,J)+URFP*(PP(I,J)-PPREF)
70    CONTINUE
      RETURN
	END
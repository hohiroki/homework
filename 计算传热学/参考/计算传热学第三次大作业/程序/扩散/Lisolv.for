      SUBROUTINE LISOLV(IS,IE,JS,JE,IT,JT,PHI,ISN)
C
CHAPTER  0  0  0  0  0  0  0  0  PRELIMINARIES  0  0  0  0  0  0  0  0
C
	PARAMETER(IX=350,JY=20)
      DIMENSION PHI(IT,JT),P(JY),Q(JY)
      COMMON
     1/COEF/AP(IX,JY),AN(IX,JY),AS(IX,JY),AE(IX,JY),
     1      AW(IX,JY),SU(IX,JY),SP(IX,JY)
      GO TO (1,2),ISN
Cif isn=1, obtain the new phis from north boundary to sorth boundary
Cif isn=2, obtain the new phis from sorth boundary to north boundary
C
CHAPTER**1**1**1**1**N TO S****N TO S****N TO S****N TO S**
C
1     CONTINUE
      JSM1=JS-1
	DO 10 I=IS,IE
	P(JSM1)=0.0
	Q(JSM1)=PHI(I,JSM1)
	DO 20 J=JS,JE
	TEMP=AP(I,J)-AS(I,J)*P(J-1)+1E-30
	P(J)=AN(I,J)/TEMP
	Q(J)=(AW(I,J)*PHI(I-1,J)+AE(I,J)*PHI(I+1,J)+
     1    AS(I,J)*Q(J-1)+SU(I,J))/TEMP
20    CONTINUE
      DO 30 J=JE,JS,-1
	PHI(I,J)=PHI(I,J+1)*P(J)+Q(J)
30    CONTINUE
10    CONTINUE
      RETURN
C
CHAPTER**2**2**2**2**S TO N****S TO N****S TO N****S TO N**
C
2     CONTINUE
      JEP1=JE+1
	DO 40 I=IS,IE
	P(JEP1)=0.0
	Q(JEP1)=PHI(I,JEP1)
	DO 50 J=JE,JS,-1
	TEMP=AP(I,J)-AN(I,J)*P(J+1)+1E-30
	P(J)=AS(I,J)/TEMP
	Q(J)=(AW(I,J)*PHI(I-1,J)+AE(I,J)*PHI(I+1,J)+
     1      AN(I,J)*Q(J+1)+SU(I,J))/TEMP
50    CONTINUE
      DO 60 J=JS,JE
	PHI(I,J)=PHI(I,J-1)*P(J)+Q(J)
60    CONTINUE
40    CONTINUE
      RETURN
	END
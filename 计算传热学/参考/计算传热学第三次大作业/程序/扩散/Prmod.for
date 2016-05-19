      SUBROUTINE PRMOD(NCHAP)
	INCLUDE 'HEAD'
	GO TO(1,2,3,4,5,6,7,8,9,10,11),NCHAP
C
CHAPTER  1  1  1  1  1  1  1  1  U MOMENTUM  1  1  1  1  1  1  1  1  1
C
1     CONTINUE
C-----UUUU--West Wall (the first section in X-Direction)
c      I=3
c	DO 1010 J=11,NJM1
c	U(I-1,J)=0.0
c	AW(I,J)=0.0
1010  CONTINUE
C-----UUUU--East Wall (last section in the X-Direction)
      I=NIM1
      DO 1120 J=1,NJ
      U(I+1,J)=(DEN(I-1,J)*XICUW(I)+DEN(I,J)*XICUE(I))*U(I,J)/DEN(I+1,J)
1120  CONTINUE	
C-----UUUU--Sorth Wall (symmtry axis)(the first section in Y-Direction)
      J=2
	DO 30 I=3,NIM1
	AS(I,J)=0.0
	U(I,1)=U(I,J)
30    CONTINUE	
C-----UUUU--North Wall (last section in Y-Direction)
      J=NJM1
	DO 40 I=3,NIM1
	U(I,NJ)=0.0
40    CONTINUE
      RETURN
C
CHAPTER  2  2  2  2  2  2  2  2  V MOMENTUM  2  2  2  2  2  2  2  2  2
C
2     CONTINUE
C-----VVVV--West Wall (the first section in X-Direction)
c      I=2
c	DO 70 J=12,NJM1
c	V(I-1,J)=0.0
c	AW(I,J)=0.0
70    CONTINUE
C-----VVVV--East Wall (last section in the X-Direction)
      I=NIM1
	DO 80 J=3,NJM1
	V(I+1,J)=0.0
80    CONTINUE	
C-----VVVV--Sorth Wall (symmtry axis)(the first section in Y-Direction)
      J=3
	DO 90 I=2,NIM1
	AS(I,J)=0.0
	V(I,2)=0.0
	V(I,1)=0.0
90    CONTINUE	
C-----VVVV--North Wall (last section in Y-Direction)
      J=NJM1
	DO 100 I=2,NIM1
      AN(I,J)=0.0
	V(I,NJ)=0.0
100    CONTINUE
      RETURN
C
CHAPTER  3  3  3  3  3  3  PRESSURE CORRECTION  3  3  3  3  3  3  3  3
C
3     CONTINUE
C-----PPPP--West Wall (the first section in X-Direction)
      I=2
	DO 190 J=2,NJM1
	AW(I,J)=0.0
190    CONTINUE	
C-----PPPP--East Wall (last section in the X-Direction)
      I=NIM1
	DO 200 J=2,NJM1
	AE(I,J)=0.0
200    CONTINUE	
C-----PPPP--Sorth Wall (symmtry axis)(the first section in Y-Direction)
      J=2
	DO 210 I=2,NIM1
	AS(I,J)=0.0
210    CONTINUE	
C-----PPPP--North Wall (last section in Y-Direction)
      J=NJM1
	DO 220 I=2,NIM1
	AN(I,J)=0.0
220    CONTINUE
      RETURN	
C
CHAPTER  4  4  4  4  4  TURBULENT KINETIC ENERGY  4  4  4  4  4  4  4  4
C
4     CONTINUE
C-----TETE--West Wall (the first section in X-Direction)
c      I=2
c      DO 250 J=11,NJM1
c	TEWALL=CK*((U(I,J)*XICCW(I)+U(I+1,J)*XICCE(I))**2+
c     1           (V(I,J)*XICCS(J)+V(I,J+1)*XICCN(J))**2)
c      SU(I,J)=1.E30*TEWALL
c	SP(I,J)=-1.E30
c      AW(I,J)=0.0	
250   CONTINUE		
C-----TETE--East Wall (last section in the X-Direction)
      I=NIM1
	DO 260 J=1,NJ
      TE(NI,J)=TE(I,J)
260   CONTINUE		
C-----TETE--Sorth Wall (symmtry axis)(the first section in Y-Direction)
      J=2
	DO 270 I=2,NIM1
      AS(I,J)=0.0
	TE(I,1)=TE(I,J)
270   CONTINUE	
C-----TETE--North Wall (last section in Y-Direction)
      J=NJM1
	DO 280 I=2,NIM1
	TEWALL=CK*((U(I,J)*XICCW(I)+U(I+1,J)*XICCE(I))**2+
     1           (V(I,J)*XICCS(J)+V(I,J+1)*XICCN(J))**2)
      SU(I,J)=1.E30*TEWALL
	SP(I,J)=-1.E30
      AN(I,J)=0.0
280   CONTINUE		
      RETURN
C
CHAPTER  5  5  5  5  5  5  5  5  DISSIPATION  5  5  5  5  5  5  5  5  5
C
5     CONTINUE
      CMU75=CMU**0.75
	XK=0.4183
C-----EDED--West Wall (the first section in X-Direction)
c      I=2
c      DO 310 J=11,NJM1
c      EDWALL=CMU75*TE(I,J)*SQRT(TE(I,J))/XK/(X(2)-X(1))
c	SU(I,J)=1E30*EDWALL
c	SP(I,J)=-1E30
c      AW(I,J)=0.0	
310   CONTINUE	
C-----EDED--East Wall (last section in the X-Direction)
      I=NIM1
	DO 320 J=1,NJ
	ED(NI,J)=ED(I,J)
320   CONTINUE	
C-----EDED--Sorth Wall (symmtry axis)(the first section in Y-Direction)
      J=2
	DO 330 I=2,NIM1
      AS(I,J)=0.0
	ED(I,1)=ED(I,J)
330   CONTINUE	
C-----EDED--North Wall (last section in Y-Direction)
      J=NJM1
	DO 340 I=2,NIM1
      EDWALL=CMU75*TE(I,J)*SQRT(TE(I,J))/XK/(Y(NJ)-Y(NJM1))
	SU(I,J)=1E30*EDWALL
	SP(I,J)=-1E30	
	AN(I,J)=0.0
340   CONTINUE	
      RETURN
C
CHAPTER  6  6  6  6  6  6  6  6  ENERGY  6  6  6  6  6  6  6  6  6
C
6     CONTINUE
C-----CPT CPT--West Wall (the first section in X-Direction)
c      I=2
c      DO 370 J=11,NJM1
370   CONTINUE	
C-----CPT CPT--East Wall (last section in the X-Direction)
      I=NIM1
	DO 380 J=1,NJ
	T(NI,J)=T(I,J)
	DEN(NI,J)=DEN(I,J)
380   CONTINUE	
C-----CPT CPT--Sorth Wall (symmtry axis)(the first section in Y-Direction)
      J=2
	DO 390 I=2,NIM1
      AS(I,J)=0.0
	T(I,1)=T(I,J)
	DEN(I,1)=DEN(I,J)
390   CONTINUE	
C-----CPT CPT--North Wall (last section in Y-Direction)
      J=NJM1
	DO 400 I=2,NIM1
	DEN(I,NJ)=DEN(I,J)
400   CONTINUE	
      RETURN
C
CHAPTER  7  7  7  7  7  7  7  7  CH4  7  7  7  7  7  7  7  7  7
C
7     CONTINUE
C-----CH4 CH4--West Wall (the first section in X-Direction)
c      I=2
c      DO 430 J=11,NJM1
c	AW(I,J)=0.0
430   CONTINUE	
C-----CH4 CH4--East Wall (last section in the X-Direction)
      I=NIM1
	DO 440 J=1,NJ
	YCH4(NI,J)=YCH4(I,J)
440   CONTINUE	
C-----CH4 CH4--Sorth Wall (symmtry axis)(the first section in Y-Direction)
      J=2
	DO 450 I=2,NIM1
      AS(I,J)=0.0
	YCH4(I,1)=YCH4(I,J)
450   CONTINUE	
C-----CH4 CH4--North Wall (last section in Y-Direction)
      J=NJM1
	DO 460 I=2,NIM1
	AN(I,J)=0.0
460   CONTINUE	
      RETURN
C
CHAPTER  8  8  8  8  8  8  8  8  O2  8  8  8  8  8  8  8  8  8
C
8     CONTINUE
C-----O2 O2--West Wall (the first section in X-Direction)
c      I=2
c      DO 490 J=11,NJM1
c	AW(I,J)=0.0
490   CONTINUE	
C-----O2 O2--East Wall (last section in the X-Direction)
      I=NIM1
      DO 500 J=1,NJ
      YO2(NI,J)=YO2(I,J)
500   CONTINUE	
C-----O2 O2--Sorth Wall (symmtry axis)(the first section in Y-Direction)
      J=2
      DO 510 I=2,NIM1
      AS(I,J)=0.0
      YO2(I,1)=YO2(I,J)
510   CONTINUE	
C-----O2 O2--North Wall (last section in Y-Direction)
      J=NJM1
      DO 520 I=2,NIM1
      AN(I,J)=0.0
520   CONTINUE	
      RETURN
C
CHAPTER  9  9  9  9  9  9  CO2  9  9  9  9  9  9  9  
C
9     CONTINUE
C-----CO2 CO2--West Wall (the first section in X-Direction)
c      I=2
c      DO 550 J=11,NJM1
c      AW(I,J)=0.0
550   CONTINUE	
C-----CO2 CO2--East Wall (last section in the X-Direction)
      I=NIM1
      DO 560 J=1,NJ
      YCO2(NI,J)=YCO2(I,J)
560   CONTINUE	
C-----CO2 CO2--Sorth Wall (symmtry axis)(the first section in Y-Direction)
      J=2
      DO 570 I=2,NIM1
      AS(I,J)=0.0
      YCO2(I,1)=YCO2(I,J)
570   CONTINUE	
C-----CO2 CO2--North Wall (last section in Y-Direction)
      J=NJM1
      DO 580 I=2,NIM1
      AN(I,J)=0.0
580   CONTINUE	
      RETURN
C
CHAPTER  10  10  10  10  10  10  H2O  10  10  10  10  10  10  10
C
10    CONTINUE
C-----H2O H2O--West Wall (the first section in X-Direction)
c      I=2
c      DO 610 J=11,NJM1
c      AW(I,J)=0.0
610   CONTINUE	
C-----H2O H2O--East Wall (last section in the X-Direction)
      I=NIM1
      DO 620 J=1,NJ
	YH2O(NI,J)=YH2O(I,J)
620   CONTINUE	
C-----H2O H2O --Sorth Wall (symmtry axis)(the first section in Y-Direction)
      J=2
      DO 630 I=2,NIM1
      AS(I,J)=0.0
      YH2O(I,1)=YH2O(I,J)
630   CONTINUE	
C-----H2O H2O--North Wall (last section in Y-Direction)
      J=NJM1
      DO 640 I=2,NIM1
      AN(I,J)=0.0
640   CONTINUE
      RETURN
C
CHAPTER  11  11  11  11  11  11  C  11  11  11  11  11  11  11
C
11    CONTINUE
C-----C C--West Wall (the first section in X-Direction)
c      I=2
c      DO 710 J=11,NJM1
c      AW(I,J)=0.0
710   CONTINUE	
C-----C C--East Wall (last section in the X-Direction)
      I=NIM1
      DO 720 J=1,NJ
	C(NI,J)=C(I,J)
720   CONTINUE	
C-----C C --Sorth Wall (symmtry axis)(the first section in Y-Direction)
      J=2
      DO 730 I=2,NIM1
      AS(I,J)=0.0
      C(I,1)=C(I,J)
730   CONTINUE	
C-----H2O H2O--North Wall (last section in Y-Direction)
      J=NJM1
      DO 740 I=2,NIM1
      AN(I,J)=0.0
740   CONTINUE
      RETURN
      END
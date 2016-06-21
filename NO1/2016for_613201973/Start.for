      SUBROUTINE START
      INCLUDE 'HEAD'
C                                                                       
CHAPTER  1  1  1  1  1  1  INITIAL OPERATIONS  1  1  1  1  1  1  1  1  1
C
C---------INPUT INITIAL FIELDS
      IF(.NOT.IREAD) GO TO 298                                        
      OPEN(11,FILE='OUT',STATUS='OLD')
      READ(11,*)U
      READ(11,*)V
      READ(11,*)P
	READ(11,*)DEN
	READ(11,*)T
	READ(11,*)CPT
	READ(11,*)YCH4
	READ(11,*)YO2
	READ(11,*)YCO2
	READ(11,*)YH2O
	READ(11,*)YN2
	READ(11,*)VIS
	READ(11,*)VSLL
      READ(11,*)PP
      CLOSE(11)
      IF(IREAD) GO TO 228
298   CONTINUE
C---------INITIALIZE U-FIELDS
      DO 100 I=1,NI
	DO 100 J=1,NJ
	U(I,J)=UIN
	IF(I.EQ.1.OR.I.EQ.2)THEN
	IF(J.GE.11)THEN
	U(I,J)=0.0
	END IF 
	END IF
	IF(J.EQ.NJ)THEN
	U(I,J)=0.0
	END IF
100   CONTINUE
C---------INITIALIZE V-FIELDS
      DO 110 I=1,NI
	DO 110 J=1,NJ
	V(I,J)=VIN
	IF(I.EQ.1)THEN
	V(I,J)=0.0
	END IF 
	IF(J.EQ.1.OR.J.EQ.2.OR.J.EQ.NJ)THEN
	V(I,J)=0.0
	END IF
110   CONTINUE
C---------INITIALIZE C-FIELDS
      DO 120 I=1,NI
	DO 120 J=1,NJ



c     for flow case
	TE(I,J)=TEIN
	ED(I,J)=EDIN
	CPT(I,J)=(0.167*300.0+1173.0)*300.0
	T(I,J)=300.0
	YCH4(I,J)=0.00
	YO2(I,J)=0.23*1.0
	YN2(I,J)=0.77*1.0
	YCO2(I,J)=0.0
	YH2O(I,J)=0.0
	C(I,J)=0.0
	IF(I.EQ.1.AND.J.LE.10)THEN
	CPT(I,J)=(0.167*300+1173.0)*300.0
	T(I,J)=300.0
	YCH4(I,J)=0.0
	YO2(I,J)=0.23*1.0
	YN2(I,J)=0.77*1.0
	YCO2(I,J)=0.0
	YH2O(I,J)=0.0
	C(I,J)=0.0
	END IF
	IF(I.EQ.1.AND.J.GE.11)THEN
	TE(I,J)=0.0
	ED(I,J)=0.0
	CPT(I,J)=(0.167*300+1173.0)*300.0
	T(I,J)=300.0
	YCH4(I,J)=0.0
	YO2(I,J)=0.23*1.0
	YN2(I,J)=0.77*1.0
	YCO2(I,J)=0.0
	YH2O(I,J)=0.0
	C(I,J)=0.0
	END IF
	IF(J.EQ.NJ)THEN
	TE(I,J)=0.0
	ED(I,J)=0.0
	CPT(I,J)=(0.167*1000+1173.0)*1000.0
	T(I,J)=1000.0
	YCH4(I,J)=0.0
	YO2(I,J)=0.23*1.0
	YN2(I,J)=0.77*1.0
	YCO2(I,J)=0.0
	YH2O(I,J)=0.0
	C(I,J)=0.0
      END IF
	if(i.eq.1.and.j.eq.4)then
	c(i,j)=1.0
	end if
c     for flow case

c     for combustion case
	TE(I,J)=TEIN
	ED(I,J)=EDIN
	CPT(I,J)=(0.167*1700.0+1173.0)*1700.0
	T(I,J)=1700.0
	YCH4(I,J)=0.03
	YO2(I,J)=0.23*0.97
	YN2(I,J)=0.77*0.97
	YCO2(I,J)=0.0
	YH2O(I,J)=0.0
	C(I,J)=0.0
	IF(I.EQ.1.AND.J.LE.10)THEN
	CPT(I,J)=(0.167*500+1173.0)*500.0
	T(I,J)=500.0
	YCH4(I,J)=0.03
	YO2(I,J)=0.23*0.97
	YN2(I,J)=0.77*0.97
	YCO2(I,J)=0.0
	YH2O(I,J)=0.0
	C(I,J)=0.0
	END IF
	IF(I.EQ.1.AND.J.GE.11)THEN
	TE(I,J)=0.0
	ED(I,J)=0.0
	CPT(I,J)=(0.167*500+1173.0)*500.0
	T(I,J)=500.0
	YCH4(I,J)=0.0
	YO2(I,J)=0.23*1.0
	YN2(I,J)=0.77*1.0
	YCO2(I,J)=0.0
	YH2O(I,J)=0.0
	C(I,J)=0.0
	END IF
	IF(J.EQ.NJ)THEN
	TE(I,J)=0.0
	ED(I,J)=0.0
	CPT(I,J)=(0.167*500+1173.0)*500.0
	T(I,J)=500.0
	YCH4(I,J)=0.0
	YO2(I,J)=0.23*1.0
	YN2(I,J)=0.77*1.0
	YCO2(I,J)=0.0
	YH2O(I,J)=0.0
	C(I,J)=0.0
      END IF
	if(i.eq.1.and.j.eq.4)then
	c(i,j)=1.0
	end if
c     for combustin case
120   CONTINUE
228   CONTINUE
      do 15 i=1,ni
	do 15 j=1,nj
	XM=YCH4(I,J)/XMCH4+YO2(I,J)/XMO2+YCO2(I,J)/XMCO2+
     1   YH2O(I,J)/XMH2O+YN2(I,J)/XMN2
	DEN(I,J)=1.01325E5/(8314.6*T(I,J)*XM)       
15    CONTINUE
C---------CALCULATING INLET QUANTITIES
      FLOWIN=0.0                                                        
      ARDEN=0.0                                                         
      ARDENT=0.                                                         
      XMONIN=0.
	TMONIN=0.
      DO 20 J=2,NJM1
	AREAW=SSN(J)*RCV(J)
	FLOWIN=FLOWIN+DEN(1,J)*U(2,J)*AREAW
	ARDENT=ARDENT+DEN(1,J)*AREAW
	XMONIN=XMONIN+DEN(1,J)*U(2,J)*U(2,J)*AREAW
	TMONIN=TMONIN+DEN(1,J)*U(2,J)*CPT(1,J)*AREAW
20    CONTINUE
      UMEAN=FLOWIN/ARDENT 
      IF(INPROV) CALL PROPS
	RETURN
	END
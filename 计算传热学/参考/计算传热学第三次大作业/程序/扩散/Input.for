      SUBROUTINE INPUT
      INCLUDE 'HEAD'
C---------FLUID PROPERITIES                                             
      DENSIT=1.29
      VISCOS=1.86E-5
C---------TURBULENCE CONSTANTS
      CD=1.0
      CMU=0.09
      C1=1.44
 	C2=1.92
	CF=0.3
	CK=CF*SQRT(CMU)/4.
C---------PRANT NUMBER	
	PRTE=1.3
      PRED=1.0
	PRC=1.0
	PRT=1.0
	PRCH4=1.0
	PRO2=1.0
	PRCO2=1.0
	PRH2O=1.0
C---------INLET BOUNDARY VALUES
	DIN=0.095
      UIN=1.
	VIN=0.0
      TEIN=0.002*(UIN**2+VIN**2)
	EDIN=2*TEIN**1.5*CMU**0.75/DIN/0.42
C---------PRESSURE CALCULATION                                          
      IPREF=NIM1
      JPREF=2
C---------SWEEP NUMBERS
      NSWPU=1
      NSWPV=1
      NSWPP=10
      NSWPK=1                                                          
      NSWPD=1 
	NSWPT=1
	NSCH4=1
	NSWO2=1
	NSCO2=1
	NSH2O=1
	NSWC=1
C---------UNDERRELAXATION FACTORS
      URFU=0.5
      URFV=0.5
      URFP=1.0
      URFK=0.5
      URFE=0.5
	URFC=0.5
      URFVIS=0.5
	URFT=0.25
	UCH4=0.25
	URO2=0.25
	UCO2=0.25
	UH2O=0.25
C---------MOLECULAR WEIGHT
      XMCH4=16.0
	XMO2=32.0
	XMCO2=44.0
	XMH2O=18.0
	XMN2=28.0	
      RETURN
      END
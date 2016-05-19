      SUBROUTINE PRINT(NI,NJ,IT,JT,X,Y,PHI,HEAD,nfile)
	CHARACTER*4 HEAD
CA**********************************************************************
C                                                                       
      DIMENSION PHI(IT,JT),X(IT),Y(JT)
	WRITE(nfile,*)
      WRITE(nfile,200)HEAD
	ISTEP=5
	IS=1
	IE=ISTEP
10    CONTINUE
	IF(NI.LE.IE)THEN
	IE=NI
	END IF
      WRITE(nfile,210)(I,I=IS,IE)
      WRITE(nfile,220)(X(I),I=IS,IE)
      DO 30 J=NJ,1,-1
      WRITE(nfile,240)J,Y(J),(PHI(I,J),I=IS,IE)
30    CONTINUE
      WRITE(nfile,230)
      IF(IE.EQ.NI)GOTO 20
	IS=IE
	IE=IE+ISTEP
	IS=IS+1
	GOTO 10
20    CONTINUE	
      RETURN
200   FORMAT(1X,20(1H*),7X,A4,7X,20(1H*))
210   FORMAT(1X,6X,2HI=,6X,15(4X,I3,5X))
220   FORMAT(1X,6X,2HX=,3X,15(2x,f9.4,1X))
230   FORMAT(1X,2X,1HJ,2X,3X,1HY)
240   FORMAT(1X,I4,1X,F7.5,1X,15(E11.3,1X))
      END
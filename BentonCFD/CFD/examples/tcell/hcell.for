      PROGRAM MAIN
C***********************************************************************
C
C  HCELL SOLVES THE NAVIER-STOKES EQUATIONS IN AN H-CELL REGION.
C
C  DISCUSSION:
C
C    THE FLUID FLOW PROBLEM IS FORMULATED IN TERMS OF
C    PRIMITIVE VARIABLES-U,V,AND P.
C
C    U_T-LAPLACIAN U+(U.GRAD)U+GRAD P=F
C                                DIV U=0
C
C    BOUNDARY CONDITIONS:  (U,V)=(0,0)ON TOP
C                          (U,V)=0 ON LEFT,RIGHT AND BOTTOM
C
C    THIS VERSION USES FINITE ELEMENT TECHNIQUES
C    WITH PIECEWISE LINEAR FUNCTIONS ON TRIANGLES TO APPROXIMATE
C    THE PRESSURE AND QUADRATICS ON TRIANGLES FOR THE VELOCITY
C    (TAYLOR-HOOD ELEMENT),ISOPARAMETRIC ELEMENT
C
C  INPUT FILES:
C
C    UP000.TXT CONTAINS THE INITIAL VALUES OF THE SOLUTION COEFFICIENTS.
C
C  LOCAL PARAMETERS:
C
C    DOUBLE PRECISION AREA(NELEMN),THE AREA OF EACH ELEMENT.
C
C    INTEGER BC_TYPE SELECTS THE BOUNDARY CONDITIONS,BY CONTROLLING THE VALUE OF ALPHA.
C    LEGAL VALUES ARE 1,FOR A STEP FUNCTION,2 FOR A "HAT" FUNCTION,3 FOR A SINUSOID.
C
C    INTEGER INDX(MAXND,NUK),LISTS THE INDICES OF THE U,V,AND P
C    VARIABLES ASSOCIATED WITH THE NODE.  IT WOULD BE USEFUL IF THE
C    INDEX ASSOCIATED WITH PRESSURE ALSO INDICATED WHETHER THERE WAS
C    NO PRESSURE VARIABLE ASSOCIATED WITH THE NODE,OR THAT IT WAS
C    PRESCRIBED.  THIS COULD BE DONE BY ASSIGNING INDX(NODE,3)=0
C    FOR THE MIDSIDE NODES OF THE 6 NODE QUADRATIC ELEMENTS.
C
C    INTEGER MAXEL IS AN OVERESTIMATE OF THE NUMBER OF ELEMENTS.
C    INTEGER MAXND IS AN OVERESTIMATE OF THE NUMBER OF NODES.
C    INTEGER MAXUN IS AN OVERESTIMATE OF THE NUMBER OF UNKNOWNS.
C    INTEGER MINUN IS AN ESTIMATE OF THE NECESSARY BANDWIDTH OF THE SYSTEM MATRIX.
C    INTEGER MX COUNTS THE NUBER OF COLUMNS OF NODES.
C    INTEGER MY COUNTS THE NUMBER OF ROWS OF NODES.
C    INTEGER NBAND,THE BANDWIDTH FOR THE FINITE ELEMENT MATRIX.
C    INTEGER NELEMN,THE NUMBER OF ELEMENTS.
C    INTEGER NEQN,THE TOTAL NUMBER OF UNKNOWNS.
C    INTEGER NLBAND,NUBAND,THE LOWWER AND UPPER HALF BANDWIDTHS
C    FOR THE FINITE ELEMENT MATRIX.
C    INTEGER NNODES,THE NUMBER OF NODES PER ELEMENT.
C    INTEGER NODE(MAXEL,NNODES),THE NODES THAT MAKE UP EACH ELEMENT.
C    INTEGER NP,THE NUMBER OF NODES.
C    INTEGER NQUAD,THE NUMBER OF QUADRATURE POINTS USED IN ASSEMBLY.
C    (THIS IS 3)
C    INTEGER NUK,THE MAXIMUM NUMBER OF UNKNOWNS ASSOCIATED WITH ONE NODE.
C    (THIS IS 3)
C    INTEGER NX COUNTS,NOT QUITE ALL THE ELEMENTS IN THE X DIRECTION,BUT THE NUMBER
C    OF ELEMENTS PLUS 1.
C    INTEGER NX1,NX2,NX3,COUNT THE ELEMENTS IN THE X DIRECTION IN THE THREE SUBREGIONS,
C    INTEGER NY COUNTS,NOT QUITE ALL THE ELEMENTS IN THE Y DIRECTION,BUT THE NUMBER
C    OF ELEMENTS PLUS 1.
C    NY1,NY2,NY3,COUNT THE ELEMENTS IN THE Y DIRECTION IN THE THREE SUBREGIONS.
C    IREG_DENSITY_X(3),IREG_DENSITY_Y(3),SPECIFIES THE
C    DENSITY OF ELEMENTS IN THE TWO COORDINATE DIRECTIONS.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      PARAMETER(NX1=10)
      PARAMETER(NX2=10)
      PARAMETER(NX3=10)
      PARAMETER(NY1=10)
      PARAMETER(NY2=10)
      PARAMETER(NY3=10)
      PARAMETER(NX=NX1+NX2+NX3+1)
      PARAMETER(NY=NY1+NY2+NY3+1)
      PARAMETER(MX=2*NX-1)
      PARAMETER(MY=2*NY-1)
      PARAMETER(MAXEL= 2*(NX-1)*(NY-1))
      PARAMETER(MAXND= MX*MY)
      PARAMETER(MAXUN= 2*MX*MY+NX*NY)
      PARAMETER(MINUN= 27*NY)
      PARAMETER(NNODES=6)
      PARAMETER(NUK=3)
      PARAMETER(NQUAD=3)
      PARAMETER(N_TIME=10)
      DIMENSION A(MINUN,MAXUN)
      DIMENSION AREA(MAXEL)
      DIMENSION F(MAXUN)
      DIMENSION G(MAXUN)
      DIMENSION REGION_X(4)
      DIMENSION REGION_Y(4)
      DIMENSION UOLD(MAXUN)
      DIMENSION XC(MAXND)
      DIMENSION XM(MAXEL,NQUAD)
      DIMENSION YC(MAXND)
      DIMENSION YM(MAXEL,NQUAD)
      DIMENSION INDX(MAXND,NUK)
      DIMENSION IPIVOT(MAXUN)
      DIMENSION NODE(MAXEL,NNODES)
      DIMENSION IREG_DENSITY_X(3)
      DIMENSION IREG_DENSITY_Y(3)
      LOGICAL NODE_MASK(MAXND)
C
      WRITE(*,'(A)')'HREGION VERSION 4:'
      WRITE(*,'(A)')'SOLVE THE NAVIER STOKES FLUID FLOW'
      WRITE(*,'(A)')'EQUATIONS IN AN H-SHAPED REGION,'
      WRITE(*,'(A)')'USING FINITE ELEMENTS.'
      IBC_TYPE=1
C
C  THE DENSITY OF ELEMENTS IN EACH REGION IS ASSUMED TO BE SOME MULTIPLE OF THE BASE.
C  FOR NOW,WE WILL JUST USE THE BASE DENSITIES.
C
      NREFINE=1
      IREG_DENSITY_X(1)=NX1*NREFINE
      IREG_DENSITY_X(2)=NX2*NREFINE
      IREG_DENSITY_X(3)=NX3*NREFINE
      IREG_DENSITY_Y(1)=NY1*NREFINE
      IREG_DENSITY_Y(2)=NY2*NREFINE
      IREG_DENSITY_Y(3)=NY3*NREFINE
      WRITE(*,'(A)')'THE X DIRECTION IS DIVIDED INTO THREE'
      WRITE(*,'(A)')'REGIONS,WITH ELEMENT DENSITIES:'
      WRITE(*,'(4X,3I6)')( IREG_DENSITY_X(I),I=1,3)
      WRITE(*,'(A,I6)')'CORRESPONDING NX=',NX
      WRITE(*,'(A)')'THE Y DIRECTION IS DIVIDED INTO THREE'
      WRITE(*,'(A)')'REGIONS,WITH ELEMENT DENSITIES:'
      WRITE(*,'(4X,3I6)')( IREG_DENSITY_Y(I),I=1,3)
      WRITE(*,'(A,I6)')'CORRESPONDING NY=',NY
C
C  DEFINE THE BREAKPOINTS THAT DIVIDE THE REGION INTO 9 LOGICAL BLOCKS.
C
      REGION_X(1)=0D0
      REGION_X(2)=3D0
      REGION_X(3)=6D0
      REGION_X(4)=9D0
      REGION_Y(1)=0D0
      REGION_Y(2)=3D0
      REGION_Y(3)=6D0
      REGION_Y(4)=9D0
      WRITE(*,'(A)')'THE X SUBREGIONS ARE DEMARKED BY 4 VALUES:'
      WRITE(*,'(4X,4F10.4)')( REGION_X(I),I=1,4)
      WRITE(*,'(A)')'THE Y SUBREGIONS ARE DEMARKED BY 4 VALUES:'
      WRITE(*,'(4X,4F10.4)')( REGION_Y(I),I=1,4)
      REYNLD=1.0D0
      MROW1=MINUN
      NSIM=3
      NSTEPS=15
      TOLNS=1D-6
      TOLOPT=1D-6
      PI=4.D0*DATAN(1.0D0)
      WRITE(*,'(A,I6)')'MAXIMUM NUMBER OF NODES=   ',MAXND
      WRITE(*,'(A,I6)')'MAXIMUM NUMBER OF ELEMENTS=',MAXEL
      WRITE(*,'(A,I6)')'MAXIMUM NUMBER OF UNKNOWNS=',MAXUN
      WRITE(*,'(A,I6)')'MAXIMUM MATRIX DIMENSION 1=',MINUN
C
C  SETGRD CONSTRUCTS GRID,NUMBERS UNKNOWNS,CALCULATES AREAS,
C  AND POINTS FOR MIDPOINT QUADRATURE RULE,BANDWIDTH AND NEQN
C
      CALL SETGRD(XC,YC,AREA,XM,YM,IREG_DENSITY_X,
     &   IREG_DENSITY_Y,REGION_X,REGION_Y,
     &   NODE,INDX,NLBAND,NUBAND,NBAND,NELEMN,NP,
     &   NNODES,NUK,NQUAD,NEQN,MAXND,MAXEL)
      WRITE(*,'(A,I6)')'NUMBER OF NODES=   ',NP
      WRITE(*,'(A,I6)')'NUMBER OF ELEMENTS=',NELEMN
      WRITE(*,'(A,I6)')'NUMBER OF UNKNOWNS=',NEQN
C
C  WRITE THE ELEMENT NODE MATRIX TO A FILE.
C
      WRITE(*,*)'writing: hcell.2dv'
      OPEN(1,FILE='hcell.2dv',STATUS='UNKNOWN')
      WRITE(1,'(I4)')NP
      DO I=1,NP
        WRITE(1,'(2E13.5)')XC(I),YC(I)
      ENDDO
      WRITE(1,'(I4)')NELEMN
      DO I=1,NELEMN
        WRITE(1,*)(NODE(I,J),J=1,3)
      ENDDO
      CLOSE(1)
C
C  MAKE A PLOT OF THE NODES.
C
      NROW1=NLBAND+NLBAND+NUBAND+1
      NCOL1=NEQN
      DO I=1,NEQN
        F(I)=0.D0
      ENDDO
C
C  READ THE INITIAL VALUE OF THE SOLUTION FROM A FILE.
C
      DO I=1,NEQN
        UOLD(I)=0.5D0
      ENDDO
      DELTAT=0.0002D0
      RDEL=1.D0/DELTAT
      WRITE(*,*)'DELTA T=',DELTAT
C
C  CARRY OUT THE TIME ITERATION.
C
      DO ITER=1,N_TIME
         WRITE(*,*)'TIME STEP ',ITER
        IF(IBC_TYPE.EQ.1)THEN
          IF(ITER.LE.250)THEN
            ALPHA=5.D0
          ELSE
            ALPHA=1.D0
          ENDIF
        ELSE IF(IBC_TYPE.EQ.2)THEN
          IF(ITER.LE.250)THEN
            ALPHA=80.D0*DBLE(ITER)*DELTAT+1.D0
          ELSE
            ALPHA=-80.D0*DBLE(ITER)*DELTAT+9.D0
          ENDIF
        ELSE IF(IBC_TYPE.EQ.3)THEN
          ALPHA=2.D0*SIN(DBLE(ITER)*0.01*PI)
        ENDIF
        DO I=1,NEQN
          G(I)=F(I)
        ENDDO
        DO I=1,NEQN
          F(I)=0.D0
        ENDDO
        CALL NSTOKE(XC,YC,AREA,XM,YM,
     &     A,F,G,UOLD,REYNLD,TOLNS,XLNGTH,YLNGTH,
     &     NODE,INDX,IPIVOT,MROW1,
     &     NLBAND,NUBAND,NBAND,NROW1,NCOL1,
     &     NELEMN,NP,NNODES,NUK,NQUAD,NEQN,
     &     NSTEPS,NSIM,MAXND,MAXEL,RDEL,ALPHA)
C
C  SAVE U=(GX,GY)TO 'UE.DAT' FOR 'T.F'
C
        DO I=1,NEQN
          UOLD(I)=F(I)
        ENDDO
      ENDDO
      WRITE(*,*)'writing: hcell.v2d'
      OPEN(1,FILE='hcell.v2d',STATUS='UNKNOWN')
      WRITE(*,*)'writing: hcell.plt'
      OPEN(2,FILE='hcell.plt',STATUS='UNKNOWN')
      WRITE(2,1000)
 1000 FORMAT('TITLE="H-CELL"')
      WRITE(2,1001)
 1001 FORMAT('VARIABLES="X" "Y" "P" "U" "V"')
      WRITE(2,1002)NP,NELEMN
 1002 FORMAT('ZONE T="RESULTS" N=',I4,', E=',I4,
     &', F=FEPOINT, ET=TRIANGLE'/
     &'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)')
      DO N=1,NP
        I=INDX(N,1)
        IF(0.LT.I)THEN
          U=F(I)
        ELSE IF(I.LT.0)THEN
          U=UBDRY(1,N,I,XC,YC)
        ELSE
          U=0.0D0
        ENDIF
        I=INDX(N,2)
        IF(0.LT.I)THEN
          V=F(I)
        ELSE IF(I.LT.0)THEN
          V=UBDRY(2,N,I,XC,YC)
        ELSE
          V=0.0D0
        ENDIF
        I=INDX(N,3)
        IF(0.LT.I)THEN
          P=F(I)
        ELSE
          P=0.0D0
        ENDIF
        WRITE(1,'(4G15.6)')XC(N),YC(N),U,V
        WRITE(2,'(5G15.6)')XC(N),YC(N),P,U,V
      ENDDO
      CLOSE(1)
      DO I=1,NELEMN
        WRITE(2,*)(NODE(I,J),J=1,3)
      ENDDO
      CLOSE(2)
      STOP
      END
      SUBROUTINE SETGRD(XC,YC,AREA,XM,YM,IREG_DENSITY_X,
     &  IREG_DENSITY_Y,REGION_X,REGION_Y,
     &  NODE,INDX,NLBAND,NUBAND,NBAND,NELEMN,NP,
     &  NNODES,NUK,NQUAD,NEQN,MAXND,MAXEL)
C
C***********************************************************************
C
C  SETGRD SETS UP THE GRID FOR THE PROBLEM.
C
C  PARAMETERS:
C
C    OUTPUT,DOUBLE PRECISION XC(NP),YC(NP),THE COORDINATES OF THE NODES.
C    OUTPUT,DOUBLE PRECISION AREA(NELEMN),THE AREA OF EACH ELEMENT.
C    OUTPUT,DOUBLE PRECISION XM(MAXEL,NQUAD),YM(MAXEL,NQUAD),THE COORDINATES
C    OF QUADRATURE POINTS IN EACH ELEMENT.
C    INPUT, IREG_DENSITY_X(3),IREG_DENSITY_Y(3),SPECIFIES THE
C    DENSITY OF ELEMENTS IN THE TWO COORDINATE DIRECTIONS.
C    INPUT,DOUBLE PRECISION REGION_X(4),REGION_Y(4),THE COORDINATES OF
C    BREAKPOINTS THAT DEFINE 9 LOGICAL SUBRECTANGLES.
C    OUTPUT, NODE(MAXEL,NNODES),THE NODES THAT MAKE UP EACH ELEMENT.
C    OUTPUT, INDX(MAXND,NUK),LISTS THE INDICES OF THE U,V,AND P
C    VARIABLES ASSOCIATED WITH THE NODE.
C    OUTPUT, NLBAND,NUBAND,THE LOWWER AND UPPER HALF BANDWIDTHS
C    FOR THE FINITE ELEMENT MATRIX.
C    OUTPUT, NBAND,THE BANDWIDTH FOR THE FINITE ELEMENT MATRIX.
C    OUTPUT, NELEMN,THE NUMBER OF ELEMENTS.
C    OUTPUT, NP,THE NUMBER OF NODES.
C    INPUT, NNODES,THE NUMBER OF NODES PER ELEMENT.
C    INPUT, NUK,THE MAXIMUM NUMBER OF UNKNOWNS ASSOCIATED WITH ONE NODE.
C    INPUT, NQUAD,THE NUMBER OF QUADRATURE POINTS.
C    OUTPUT, NEQN,THE TOTAL NUMBER OF UNKNOWNS.
C    INPUT, MAXND,THE MAXIMUM NUMBER OF NODES.
C    INPUT, MAXEL,THE MAXIMUM NUMBER OF ELEMENTS.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION AREA(MAXEL)
      DIMENSION INDX(MAXND,NUK)
      DIMENSION IREG_DENSITY_X(3)
      DIMENSION IREG_DENSITY_Y(3)
      DIMENSION NODE(MAXEL,NNODES)
      DIMENSION REGION_X(4)
      DIMENSION REGION_Y(4)
      DIMENSION XC(MAXND)
      DIMENSION XM(MAXEL,NQUAD)
      DIMENSION YC(MAXND)
      DIMENSION YM(MAXEL,NQUAD)
      INTEGER DOF_COUNT
C
C  DETERMINE NP,THE NUMBER OF NODES.
C
      CALL HCELL_NODE_COUNT(IREG_DENSITY_X,IREG_DENSITY_Y,
     &  NP)
C
C  ASSIGN COORDINATES TO THE NODES,XC AND YC.
C
      CALL HCELL_NODE_XY(IREG_DENSITY_X,IREG_DENSITY_Y,
     &  NP,REGION_X,REGION_Y,XC,YC)
C
C  DETERMINE NELEMN,THE NUMBER OF ELEMENTS.
C
      CALL HCELL_ELEMENT_COUNT(IREG_DENSITY_X,
     &  IREG_DENSITY_Y,NELEMN)
C
C  ASSIGN NODES TO ELEMENTS IN NODE.
C
      CALL HCELL_ELEMENT_NODE(IREG_DENSITY_X,IREG_DENSITY_Y,
     &  MAXEL,NELEMN,NNODES,NODE)
C
C  DETERMINE NEQN,THE NUMBER OF DEGREES OF FREEDOM.
C  ASSIGN DEGREES OF FREEDOM IN INDX.
C
      CALL HCELL_DOF_SET(IREG_DENSITY_X,
     &  IREG_DENSITY_Y,MAXND,NP,INDX,NEQN)
C
C  FOR THE ASSEMBLY ROUTINE,DETERMINE THE QUADRATURE DATA
C  NQUAD,AREA,XM AND YM.
C
      CALL QUAD_A_SET(MAXEL,NELEMN,NQUAD,AREA,XM,YM)
C
C  COMPUTE THE BANDWIDTHS NLBAND,NUBAND AND NBAND.
C
      CALL ELEMENT_NODE_BANDWIDTH(MAXND,MAXEL,NNODES,
     &   NELEMN,NODE,NEQN,NP,INDX,NLBAND,NUBAND,NBAND)
      RETURN
      END
      SUBROUTINE HCELL_ELEMENT_COUNT(IREG_DENSITY_X,
     &  IREG_DENSITY_Y,ELEMENT_NUM)
C
C*******************************************************************************
C
C  HCELL_ELEMENT_COUNT DETERMINES THE NUMBER OF ELEMENTS IN THE REGION.
C
C  DIAGRAM:
C
C          +----------------------------+
C          |              :     :       |
C    ROW 3 |   (3,1)     :(3,2): (3,3)|
C          |              :     :       |
C          +--------------+.....+-------+
C                         |     |
C    ROW 2     EMPTY      |(2,2)|  EMPTY
C                         |     |
C          +--------------+.....+-------+
C          |              :     :       |
C    ROW 1 |   (1,1)     :(1,2): (1,3)|
C          |              :     :       |
C          +----------------------------+
C
C              COL 1       COL 2  COL 3
C
C  DISCUSSION:
C
C    THE REGION IS DIVIDED INTO A 3 BY 3 GRID.  SUBREGION(I,J)
C    IS DIVIDED INTO ELEMENT_DENSITY_X(J)*ELEMENT_DENSITY_Y(I)SQUARES.
C    THEN EACH SQUARE IS SPLIT INTO TWO TRIANGLES,WITH THE DIAGONAL
C    GOING FROM THE UPPER LEFT TO LOWER RIGHT.
C
C  AUTHOR:
C
C    JOHN BURKARDT
C
C  PARAMETERS:
C
C    INPUT, ELEMENT_DENSITY_X(3),THE DENSITY OF ELEMENTS
C    IN THE THREE COLUMNS.
C    INPUT, ELEMENT_DENSITY_Y(3),THE DENSITY OF ELEMENTS
C    IN THE THREE ROWS.
C    OUTPUT, ELEMENT_NUM,THE NUMBER OF ELEMENTS.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER ELEMENT_NUM
      DIMENSION IREG_DENSITY_X(3)
      DIMENSION IREG_DENSITY_Y(3)
C
      ELEMENT_NUM=
     &    2*IREG_DENSITY_X(1)*IREG_DENSITY_Y(1)
     & +2*IREG_DENSITY_X(1)*IREG_DENSITY_Y(3)
     & +2*IREG_DENSITY_X(2)*IREG_DENSITY_Y(1)
     & +2*IREG_DENSITY_X(2)*IREG_DENSITY_Y(2)
     & +2*IREG_DENSITY_X(2)*IREG_DENSITY_Y(3)
     & +2*IREG_DENSITY_X(3)*IREG_DENSITY_Y(1)
     & +2*IREG_DENSITY_X(3)*IREG_DENSITY_Y(3)
      WRITE(*,'(A)')'HCELL_ELEMENT_COUNT:'
      WRITE(*,'(A,I6)')'NUMBER OF ELEMENTS=',ELEMENT_NUM
      RETURN
      END
      SUBROUTINE HCELL_ELEMENT_NODE(IREG_DENSITY_X,
     &  IREG_DENSITY_Y,MAXEL,ELEMENT_NUM,NNODES,ELEMENT_NODE)
C
C*******************************************************************************
C
C  HCELL_ELEMENT_NODE DETERMINES THE NODES THAT MAKE UP EACH ELEMENT.
C
C  DIAGRAM:
C
C          +----------------------------+
C          |              :     :       |
C    ROW 3 |   (3,1)     :(3,2): (3,3)|
C          |              :     :       |
C          +--------------+.....+-------+
C                         |     |
C    ROW 2     EMPTY      |(2,2)|  EMPTY
C                         |     |
C          +--------------+.....+-------+
C          |              :     :       |
C    ROW 1 |   (1,1)     :(1,2): (1,3)|
C          |              :     :       |
C          +----------------------------+
C
C              COL 1       COL 2  COL 3
C
C  DISCUSSION:
C
C    THE REGION IS DIVIDED INTO A 3 BY 3 GRID.  SUBREGION(I,J)
C    IS DIVIDED INTO ELEMENT_DENSITY_X(J)*ELEMENT_DENSITY_Y(I)SQUARES.
C    THEN EACH SQUARE IS SPLIT INTO TWO TRIANGLES,WITH THE DIAGONAL
C    GOING FROM THE UPPER LEFT TO LOWER RIGHT.
C
C  AUTHOR:
C
C    JOHN BURKARDT
C
C  PARAMETERS:
C
C    INPUT, ELEMENT_DENSITY_X(3),THE DENSITY OF ELEMENTS
C    IN THE THREE COLUMNS.
C    INPUT, ELEMENT_DENSITY_Y(3),THE DENSITY OF ELEMENTS
C    IN THE THREE ROWS.
C    INPUT, MAXEL,THE MAXIMUM NUMBER OF ELEMENTS.
C    INPUT, ELEMENT_NUM,THE NUMBER OF ELEMENTS.
C    INPUT, NNODES,THE NUMBER OF NODES PER ELEMENT.
C    OUTPUT, ELEMENT_NODE(MAXEL,NNODES),THE NODES THAT MAKE UP
C    EACH ELEMENT.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION IREG_DENSITY_X(3)
      DIMENSION IREG_DENSITY_Y(3)
      INTEGER COL
      INTEGER COL2
      INTEGER ELEMENT
      INTEGER ELEMENT_NODE(MAXEL,NNODES)
      INTEGER ELEMENT_NUM
      INTEGER ROW
      INTEGER ROW2
C
      ELEMENT=0
      DO COL=1,3
        DO COL2=1,IREG_DENSITY_X(COL)
          DO ROW=1,3
            IF(ROW.NE.2.OR.COL.EQ.2)THEN
              IF(COL.EQ.1)THEN
                IF(COL2.LT.IREG_DENSITY_X(1))THEN
                  IF(ROW.EQ.1)THEN
                    IF(COL2.EQ.1)THEN
                      NODE_SW=1
                    ELSE
                      NODE_SW=NODE_SW+INC1+1
                    ENDIF
                  ELSE
                    NODE_SW=NODE_SW+1
                  ENDIF
                  INC1=( 2*IREG_DENSITY_Y(1)+1)
     &                 +( 2*IREG_DENSITY_Y(3)+1)
                  INC2=( 2*IREG_DENSITY_Y(1)+1)
     &                 +( 2*IREG_DENSITY_Y(3)+1)
                ELSE IF(ROW.EQ.1)THEN
                  NODE_SW=NODE_SW+INC1+1
                   INC1=( 2*IREG_DENSITY_Y(1)+1)
     &                 +( 2*IREG_DENSITY_Y(3)+1)
                   INC2=( 2*IREG_DENSITY_Y(1)+1)
     &                 +( 2*IREG_DENSITY_Y(3)+1)
                ELSE IF(ROW.EQ.3)THEN
                  NODE_SW=NODE_SW+1
                  INC1=( 2*IREG_DENSITY_Y(1)+1)
     &                 +( 2*IREG_DENSITY_Y(3)+1)
                  INC2=( 2*IREG_DENSITY_Y(1)+1)
     &                 +( 2*IREG_DENSITY_Y(2)-1)
     &                 +( 2*IREG_DENSITY_Y(3)+1)
                ENDIF
              ELSE IF(COL.EQ.2)THEN
                IF(ROW.EQ.1)THEN
                  NODE_SW=NODE_SW+INC1+1
                ENDIF
                INC1=( 2*IREG_DENSITY_Y(1)+1)
     &                    +( 2*IREG_DENSITY_Y(2)-1)
     &               +( 2*IREG_DENSITY_Y(3)+1)
                INC2=( 2*IREG_DENSITY_Y(1)+1)
     &                +( 2*IREG_DENSITY_Y(2)-1)
     &                +( 2*IREG_DENSITY_Y(3)+1)
              ELSE IF(COL.EQ.3)THEN
                IF(1.LT.COL2)THEN
                  IF(ROW.EQ.1)THEN
                    NODE_SW=NODE_SW+INC1+1
                  ELSE
                    NODE_SW=NODE_SW+1
                  ENDIF
                  INC1=( 2*IREG_DENSITY_Y(1)+1)
     &                  +( 2*IREG_DENSITY_Y(3)+1)
                  INC2=( 2*IREG_DENSITY_Y(1)+1)
     &                  +( 2*IREG_DENSITY_Y(3)+1)
                ELSE IF(ROW.EQ.1)THEN
                  NODE_SW=NODE_SW+INC1+1
                  INC1=( 2*IREG_DENSITY_Y(1)+1)
     &                  +( 2*IREG_DENSITY_Y(2)-1)
     &                  +( 2*IREG_DENSITY_Y(3)+1)
                  INC2=( 2*IREG_DENSITY_Y(1)+1)
     &                  +( 2*IREG_DENSITY_Y(3)+1)
                ELSE IF(ROW.EQ.3)THEN
                  NODE_SW=NODE_SW
     &              +( 2*IREG_DENSITY_Y(2)-1)+1
                  INC1=( 2*IREG_DENSITY_Y(1)+1)
     &                 +( 2*IREG_DENSITY_Y(3)+1)
                  INC2=( 2*IREG_DENSITY_Y(1)+1)
     &                 +( 2*IREG_DENSITY_Y(3)+1)
                ENDIF
              ENDIF
              DO ROW2=1,IREG_DENSITY_Y(ROW)
                ELEMENT=ELEMENT+1
                ELEMENT_NODE(ELEMENT,1)=NODE_SW
                ELEMENT_NODE(ELEMENT,2)=NODE_SW+INC1+INC2
                ELEMENT_NODE(ELEMENT,3)=NODE_SW              +2
                ELEMENT_NODE(ELEMENT,4)=NODE_SW+INC1
                ELEMENT_NODE(ELEMENT,5)=NODE_SW+INC1       +1
                ELEMENT_NODE(ELEMENT,6)=NODE_SW              +1
                ELEMENT=ELEMENT+1
                ELEMENT_NODE(ELEMENT,1)=NODE_SW+INC1+INC2+2
                ELEMENT_NODE(ELEMENT,2)=NODE_SW              +2
                ELEMENT_NODE(ELEMENT,3)=NODE_SW+INC1+INC2
                ELEMENT_NODE(ELEMENT,4)=NODE_SW+INC1       +2
                ELEMENT_NODE(ELEMENT,5)=NODE_SW+INC1       +1
                ELEMENT_NODE(ELEMENT,6)=NODE_SW+INC1+INC2+1
                NODE_SW=NODE_SW+2
              ENDDO
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE HCELL_NODE_COUNT(IREG_DENSITY_X,IREG_DENSITY_Y,
     &  NODE_NUM)
C
C*******************************************************************************
C
C  HCELL_NODE_NUM DETERMINES THE NUMBER OF NODES IN THE REGION.
C
C  DIAGRAM:
C
C          +----------------------------+
C          |              :     :       |
C    ROW 3 |   (3,1)     :(3,2): (3,3)|
C          |              :     :       |
C          +--------------+.....+-------+
C                         |     |
C    ROW 2     EMPTY      |(2,2)|  EMPTY
C                         |     |
C          +--------------+.....+-------+
C          |              :     :       |
C    ROW 1 |   (1,1)     :(1,2): (1,3)|
C          |              :     :       |
C          +----------------------------+
C
C              COL 1       COL 2  COL 3
C
C  DISCUSSION:
C
C    THE REGION IS DIVIDED INTO A 3 BY 3 GRID.  SUBREGION(I,J)
C    IS DIVIDED INTO ELEMENT_DENSITY_X(J)*ELEMENT_DENSITY_Y(I)SQUARES.
C    THEN EACH SQUARE IS SPLIT INTO TWO TRIANGLES,WITH THE DIAGONAL
C    GOING FROM THE UPPER LEFT TO LOWER RIGHT.
C
C  AUTHOR:
C
C    JOHN BURKARDT
C
C  PARAMETERS:
C
C    INPUT, ELEMENT_DENSITY_X(3),THE DENSITY OF ELEMENTS
C    IN THE THREE COLUMNS.
C    INPUT, ELEMENT_DENSITY_Y(3),THE DENSITY OF ELEMENTS
C    IN THE THREE ROWS.
C    OUTPUT, NODE_NUM,THE NUMBER OF NODES.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION IREG_DENSITY_X(3)
      DIMENSION IREG_DENSITY_Y(3)
      NODE_NUM=
     &   (2*IREG_DENSITY_X(1)+1)
     & *( 2*IREG_DENSITY_Y(1)+1)
     & +( 2*IREG_DENSITY_X(1)+1)
     & *( 2*IREG_DENSITY_Y(3)+1)
     & +( 2*IREG_DENSITY_X(2)-1)
     & *( 2*IREG_DENSITY_Y(1)+1)
     & +( 2*IREG_DENSITY_X(2)+1)
     & *( 2*IREG_DENSITY_Y(2)-1)
     & +( 2*IREG_DENSITY_X(2)-1)
     & *( 2*IREG_DENSITY_Y(3)+1)
     & +( 2*IREG_DENSITY_X(3)+1)
     & *( 2*IREG_DENSITY_Y(1)+1)
     & +( 2*IREG_DENSITY_X(3)+1)
     & *( 2*IREG_DENSITY_Y(3)+1)
      WRITE(*,'(A)')'HCELL_NODE_COUNT:'
      WRITE(*,'(A,I6)')'NUMBER OF NODES=',NODE_NUM
      RETURN
      END
      SUBROUTINE HCELL_NODE_XY(IREG_DENSITY_X,IREG_DENSITY_Y,
     &  NODE_NUM,REGION_X,REGION_Y,X_NODE,Y_NODE)
C
C*******************************************************************************
C
C  HCELL_NODE_XY ASSIGNS COORDINATES TO EACH NODE.
C
C  DIAGRAM:
C
C          +----------------------------+
C          |              :     :       |
C    ROW 3 |   (3,1)      :(3,2): (3,3) |
C          |              :     :       |
C          +--------------+.....+-------+
C                         |     |
C    ROW 2     EMPTY      |(2,2)|  EMPTY
C                         |     |
C          +--------------+.....+-------+
C          |              :     :       |
C    ROW 1 |   (1,1)      :(1,2): (1,3) |
C          |              :     :       |
C          +----------------------------+
C
C              COL 1       COL 2  COL 3
C
C  DISCUSSION:
C
C    THE REGION IS DIVIDED INTO A 3 BY 3 GRID.  SUBREGION(I,J)
C    IS DIVIDED INTO ELEMENT_DENSITY_X(J)*ELEMENT_DENSITY_Y(I)SQUARES.
C    THEN EACH SQUARE IS SPLIT INTO TWO TRIANGLES,WITH THE DIAGONAL
C    GOING FROM THE UPPER LEFT TO LOWER RIGHT.
C
C  AUTHOR:
C
C    JOHN BURKARDT
C
C  PARAMETERS:
C
C    INPUT, IREG_DENSITY_X(3),THE DENSITY OF ELEMENTS
C    IN THE THREE COLUMNS.
C    INPUT, IREG_DENSITY_Y(3),THE DENSITY OF ELEMENTS
C    IN THE THREE ROWS.
C    INPUT, NODE_NUM,THE NUMBER OF NODES.
C    INPUT,DOUBLE PRECISION REGION_X(4),REGION_Y(4),THE COORDINATES OF
C    BREAKPOINTS THAT DEFINE 9 LOGICAL SUBRECTANGLES.
C    OUTPUT,DOUBLE PRECISION X_NODE(NODE_NUM),Y_NODE(NODE_NUM),
C    THE X AND Y COORDINATES OF THE NODES.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION IREG_DENSITY_X(3)
      DIMENSION IREG_DENSITY_Y(3)
      DIMENSION REGION_X(4)
      DIMENSION REGION_Y(4)
      DIMENSION X_NODE(NODE_NUM)
      DIMENSION Y_NODE(NODE_NUM)
      INTEGER COL
      INTEGER ROW
C
      NODE=0
      J=0
C
C  WORKING IN COLUMN 1,EXCEPT FOR THE EXTREME RIGHT.
C
      DO COL=1,2*IREG_DENSITY_X(1)
        J=J+1
        I=0
C
C  WORKING IN ROW 1.
C
C  +--+--+--+
C  |        |
C  +--+  +--+
C     |  |
C  +--+  +--+
C  |11      |
C  +--+--+--+
C
        DO ROW=1,2*IREG_DENSITY_Y(1)+1
          I=I+1
          NODE=NODE+1
          X_NODE(NODE)=
     &     (DBLE(2*IREG_DENSITY_X(1)+1-COL)*REGION_X(1)
     &     +DBLE(                       -1+COL)*REGION_X(2))
     &     /DBLE(2*IREG_DENSITY_X(1)         )
          Y_NODE(NODE)=
     &     (DBLE(2*IREG_DENSITY_Y(1)+1-ROW)*REGION_Y(1)
     &     +DBLE(                       -1+ROW)*REGION_Y(2))
     &     /DBLE(2*IREG_DENSITY_Y(1)         )
        ENDDO
C
C  WORKING IN ROW 3.
C
C  +--+--+--+
C  |31      |
C  +--+  +--+
C     |  |
C  +--+  +--+
C  |        |
C  +--+--+--+
C
        I=I+2*IREG_DENSITY_Y(2)-1
        DO ROW=1,2*IREG_DENSITY_Y(3)+1
          I=I+1
          NODE=NODE+1
          X_NODE(NODE)=
     &     (DBLE(2*IREG_DENSITY_X(1)+1-COL)*REGION_X(1)
     &     +DBLE(                       -1+COL)*REGION_X(2))
     &     /DBLE(2*IREG_DENSITY_X(1)         )
          Y_NODE(NODE)=
     &     (DBLE(2*IREG_DENSITY_Y(3)+1-ROW)*REGION_Y(3)
     &     +DBLE(                       -1+ROW)*REGION_Y(4))
     &     /DBLE(2*IREG_DENSITY_Y(3)         )
        ENDDO
      ENDDO
C
C  WORKING IN COLUMN 2,INCLUDING EXTREME LEFT AND RIGHT.
C
      DO COL=1,2*IREG_DENSITY_X(2)+1
C
C  WORKING IN ROW 1.
C
C  +--+--+--+
C  |        |
C  +--+  +--+
C     |  |
C  +--+  +--+
C  |   12   |
C  +--+--+--+
C
        J=J+1
        I=0
        DO ROW=1,2*IREG_DENSITY_Y(1)+1
          I=I+1
          NODE=NODE+1
          X_NODE(NODE)=
     &     (DBLE(2*IREG_DENSITY_X(2)+1-COL)*REGION_X(2)
     &     +DBLE(                       -1+COL)*REGION_X(3))
     &     /DBLE(2*IREG_DENSITY_X(2)         )
          Y_NODE(NODE)=
     &     (DBLE(2*IREG_DENSITY_Y(1)+1-ROW)*REGION_Y(1)
     &     +DBLE(                       -1+ROW)*REGION_Y(2))
     &     /DBLE(2*IREG_DENSITY_Y(1)         )
        ENDDO
C
C  WORKING IN ROW 2.
C
C  +--+--+--+
C  |        |
C  +--+  +--+
C     |22|
C  +--+  +--+
C  |        |
C  +--+--+--+
C
        DO ROW=2,2*IREG_DENSITY_Y(2)
          I=I+1
          NODE=NODE+1
          X_NODE(NODE)=
     &     (DBLE(2*IREG_DENSITY_X(2)+1-COL)*REGION_X(2)
     &     +DBLE(                       -1+COL)*REGION_X(3))
     &     /DBLE(2*IREG_DENSITY_X(2)         )
          Y_NODE(NODE)=
     &     (DBLE(2*IREG_DENSITY_Y(2)+1-ROW)*REGION_Y(2)
     &     +DBLE(                       -1+ROW)*REGION_Y(3))
     &     /DBLE(2*IREG_DENSITY_Y(2)         )
        ENDDO
C
C  WORKING IN ROW 3.
C
C  +--+--+--+
C  |   32   |
C  +--+  +--+
C     |  |
C  +--+  +--+
C  |        |
C  +--+--+--+
C
        DO ROW=1,2*IREG_DENSITY_Y(3)+1
          I=I+1
          NODE=NODE+1
          X_NODE(NODE)=
     &     (DBLE(2*IREG_DENSITY_X(2)+1-COL)*REGION_X(2)
     &     +DBLE(                       -1+COL)*REGION_X(3))
     &     /DBLE(2*IREG_DENSITY_X(2)         )
          Y_NODE(NODE)=
     &     (DBLE(2*IREG_DENSITY_Y(3)+1-ROW)*REGION_Y(3)
     &     +DBLE(                       -1+ROW)*REGION_Y(4))
     &     /DBLE(2*IREG_DENSITY_Y(3)         )
        ENDDO
      ENDDO
C
C  WORKING IN COLUMN 3,EXCEPT FOR EXTREME LEFT.
C
      DO COL=2,2*IREG_DENSITY_X(3)+1
C
C  WORKING IN ROW 1.
C
C  +--+--+--+
C  |        |
C  +--+  +--+
C     |  |
C  +--+  +--+
C  |      13|
C  +--+--+--+
C
        J=J+1
        I=0
        DO ROW=1,2*IREG_DENSITY_Y(1)+1
          I=I+1
          NODE=NODE+1
          X_NODE(NODE)=
     &     (DBLE(2*IREG_DENSITY_X(3)+1-COL)*REGION_X(3)
     &     +DBLE(                       -1+COL)*REGION_X(4))
     &     /DBLE(2*IREG_DENSITY_X(3)         )
          Y_NODE(NODE)=
     &     (DBLE(2*IREG_DENSITY_Y(1)+1-ROW)*REGION_Y(1)
     &     +DBLE(                       -1+ROW)*REGION_Y(2))
     &     /DBLE(2*IREG_DENSITY_Y(1)         )
        ENDDO
        I=I+2*IREG_DENSITY_Y(2)-1
C
C  WORKING IN ROW 3.
C
C  +--+--+--+
C  |      33|
C  +--+  +--+
C     |  |
C  +--+  +--+
C  |        |
C  +--+--+--+
C
        DO ROW=1,2*IREG_DENSITY_Y(3)+1
          I=I+1
          NODE=NODE+1
          X_NODE(NODE)=
     &     (DBLE(2*IREG_DENSITY_X(3)+1-COL)*REGION_X(3)
     &     +DBLE(                       -1+COL)*REGION_X(4))
     &     /DBLE(2*IREG_DENSITY_X(3)         )
          Y_NODE(NODE)=
     &     (DBLE(2*IREG_DENSITY_Y(3)+1-ROW)*REGION_Y(3)
     &     +DBLE(                       -1+ROW)*REGION_Y(4))
     &     /DBLE(2*IREG_DENSITY_Y(3)         )
        ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE HCELL_DOF_COUNT(IREG_DENSITY_X,IREG_DENSITY_Y,
     &  DOF_NUM)
C
C*******************************************************************************
C
C  HCELL_DOF_COUNT DETERMINES THE NUMBER OF DEGREES OF FREEDOM IN THE REGION.
C
C  DISCUSSION:
C
C    THE COUNT PRODUCED BY THIS ROUTINE DOES NOT CORRESPOND TO CURRENT PRACTICE,
C    BECAUSE DEGREES OF FREEDOM ASSOCIATED WITH BOUNDARY CONDITIONS AND OTHE
C    SPECIFICATIONS ARE NOW NOT BEING COUNTED.
C
C  DIAGRAM:
C
C          +----------------------------+
C          |              :     :       |
C    ROW 3 |   (3,1)     :(3,2): (3,3)|
C          |              :     :       |
C          +--------------+.....+-------+
C                         |     |
C    ROW 2     EMPTY      |(2,2)|  EMPTY
C                         |     |
C          +--------------+.....+-------+
C          |              :     :       |
C    ROW 1 |   (1,1)     :(1,2): (1,3)|
C          |              :     :       |
C          +----------------------------+
C
C              COL 1       COL 2  COL 3
C
C  DISCUSSION:
C
C    THE REGION IS DIVIDED INTO A 3 BY 3 GRID.  SUBREGION(I,J)
C    IS DIVIDED INTO ELEMENT_DENSITY_X(J)*ELEMENT_DENSITY_Y(I)SQUARES.
C    THEN EACH SQUARE IS SPLIT INTO TWO TRIANGLES,WITH THE DIAGONAL
C    GOING FROM THE UPPER LEFT TO LOWER RIGHT.
C
C    EACH ELEMENT,IN TURN,IS MADE UP OF 6 NODES.  THE CORNER
C    NODES HAVE 3 DEGREES OF FREEDOM(U,V,AND P),WHILE THE
C    SIDE NODES HAVE 2 DEGREES OF FREEDOM (U AND V ONLY).
C
C  AUTHOR:
C
C    JOHN BURKARDT
C
C  PARAMETERS:
C
C    INPUT, ELEMENT_DENSITY_X(3),THE DENSITY OF ELEMENTS
C    IN THE THREE COLUMNS.
C    INPUT, ELEMENT_DENSITY_Y(3),THE DENSITY OF ELEMENTS
C    IN THE THREE ROWS.
C    OUTPUT, DOF_NUM,THE NUMBER OF DEGREES OF FREEDOM.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION IREG_DENSITY_X(3)
      DIMENSION IREG_DENSITY_Y(3)
      INTEGER DOF_NUM
C
      NODE_NUM_U=
     &   (2*IREG_DENSITY_X(1)+1)
     & *( 2*IREG_DENSITY_Y(1)+1)
     & +( 2*IREG_DENSITY_X(1)+1)
     & *( 2*IREG_DENSITY_Y(3)+1)
     & +( 2*IREG_DENSITY_X(2)-1)
     & *( 2*IREG_DENSITY_Y(1)+1)
     & +( 2*IREG_DENSITY_X(2)+1)
     & *( 2*IREG_DENSITY_Y(2)-1)
     & +( 2*IREG_DENSITY_X(2)-1)
     & *( 2*IREG_DENSITY_Y(3)+1)
     & +( 2*IREG_DENSITY_X(3)+1)
     & *( 2*IREG_DENSITY_Y(1)+1)
     & +( 2*IREG_DENSITY_X(3)+1)
     & *( 2*IREG_DENSITY_Y(3)+1)
      NODE_NUM_P=
     &   (IREG_DENSITY_X(1)+1)
     & *( IREG_DENSITY_Y(1)+1)
     & +( IREG_DENSITY_X(1)+1)
     & *( IREG_DENSITY_Y(3)+1)
     & +( IREG_DENSITY_X(2)-1)
     & *( IREG_DENSITY_Y(1)+1)
     & +( IREG_DENSITY_X(2)+1)
     & *( IREG_DENSITY_Y(2)-1)
     & +( IREG_DENSITY_X(2)-1)
     & *( IREG_DENSITY_Y(3)+1)
     & +( IREG_DENSITY_X(3)+1)
     & *( IREG_DENSITY_Y(1)+1)
     & +( IREG_DENSITY_X(3)+1)
     & *( IREG_DENSITY_Y(3)+1)
      DOF_NUM=2*NODE_NUM_U+NODE_NUM_P
      WRITE(*,'(A)')'HCELL_DOF_COUNT:'
      WRITE(*,'(A,I6)')'NUMBER OF DEGREES OF FREEDOM=',DOF_NUM
      RETURN
      END
      SUBROUTINE HCELL_DOF_SET(IREG_DENSITY_X,
     &  IREG_DENSITY_Y,MAXND,NODE_NUM,NODE_DOF_INDEX,DOF)
C
C*******************************************************************************
C
C  HCELL_DOF_SET ASSIGNS DEGREES OF FREEDOM TO EACH NODE.
C
C  DIAGRAM:
C
C          +----------------------------+
C          |              :     :       |
C    ROW 3 |   (3,1)     :(3,2): (3,3)|
C          |              :     :       |
C          +--------------+.....+-------+
C                         |     |
C    ROW 2     EMPTY      |(2,2)|  EMPTY
C                         |     |
C          +--------------+.....+-------+
C          |              :     :       |
C    ROW 1 |   (1,1)     :(1,2): (1,3)|
C          |              :     :       |
C          +----------------------------+
C
C              COL 1       COL 2  COL 3
C
C  DISCUSSION:
C
C    THE REGION IS DIVIDED INTO A 3 BY 3 GRID.  SUBREGION(I,J)
C    IS DIVIDED INTO ELEMENT_DENSITY_X(J)*ELEMENT_DENSITY_Y(I)RECTANGLES.
C    THEN EACH RECTANGLE IS SPLIT INTO TWO TRIANGLES,WITH THE DIAGONAL
C    GOING FROM THE UPPER LEFT TO LOWER RIGHT.
C
C  AUTHOR:
C
C    JOHN BURKARDT
C
C  PARAMETERS:
C
C    INPUT, IREG_DENSITY_X(3),THE DENSITY OF ELEMENTS
C    IN THE THREE COLUMNS.
C    INPUT, IREG_DENSITY_Y(3),THE DENSITY OF ELEMENTS
C    IN THE THREE ROWS.
C    INPUT, MAXND,THE MAXIMUM NUMBER OF NODES.
C    INPUT, NODE_NUM,THE NUMBER OF NODES.
C    OUTPUT, NODE_DOF_INDEX(NODE_NUM,3),THE NODAL DEGREES OF FREEDOM.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION IREG_DENSITY_X(3)
      DIMENSION IREG_DENSITY_Y(3)
      DIMENSION NODE_DOF_INDEX(MAXND,3)
      INTEGER COL
      INTEGER DOF
      INTEGER ROW
      LOGICAL FOUND
C
      NODE=0
      DOF=0
      J=0
C
C  WORKING IN COLUMN 1,AND NOT HANDLING NODES IN EXTREME RIGHT.
C
      DO COL=1,2*IREG_DENSITY_X(1)
C
C  WORKING IN ROW 1.
C
C  +--+--+--+
C  |        |
C  +--+  +--+
C     |  |
C  +--+  +--+
C  |11      |
C  +--+--+--+
C
        J=J+1
        I=0
        DO ROW=1,2*IREG_DENSITY_Y(1)+1
          I=I+1
          NODE=NODE+1
       IF(COL.EQ.1)THEN
         NODE_DOF_INDEX(NODE,1)=0
         NODE_DOF_INDEX(NODE,2)=0
         IF(ROW.GE.(2*2*IREG_DENSITY_Y(1))/5+1
     &        .AND.
     &        ROW.LE.(3*2*IREG_DENSITY_Y(1))/5+1)THEN
            NODE_DOF_INDEX(NODE,1)=-3
         ENDIF
       ELSE
         IF(ROW.EQ.1)THEN
            NODE_DOF_INDEX(NODE,1)=0
            NODE_DOF_INDEX(NODE,2)=0
         ELSE
           IF(ROW.EQ.2*IREG_DENSITY_Y(1)+1)THEN
             NODE_DOF_INDEX(NODE,1)=0
             NODE_DOF_INDEX(NODE,2)=0
             IF(COL.GE.(10*2*IREG_DENSITY_X(1))/45+1
     &            .AND.
     &            COL.LE.(11*2*IREG_DENSITY_X(1))/45+1)THEN
                NODE_DOF_INDEX(NODE,2)=-3
             ENDIF
             IF(COL.GE.(22*2*IREG_DENSITY_X(1))/45+1
     &            .AND.
     &            COL.LE.(23*2*IREG_DENSITY_X(1))/45+1)THEN
                NODE_DOF_INDEX(NODE,2)=-3
             ENDIF
             IF(COL.GE.(34*2*IREG_DENSITY_X(1))/45+1
     &            .AND.
     &            COL.LE.(35*2*IREG_DENSITY_X(1))/45+1)THEN
                NODE_DOF_INDEX(NODE,2)=-3
             ENDIF
           ELSE
             DOF=DOF+1
             NODE_DOF_INDEX(NODE,1)=DOF
             DOF=DOF+1
             NODE_DOF_INDEX(NODE,2)=DOF
           ENDIF
         ENDIF
       ENDIF
          IF(MOD(J,2).EQ.1.AND.MOD(I,2).EQ.1)THEN
            DOF=DOF+1
            NODE_DOF_INDEX(NODE,3)=DOF
          ELSE
            NODE_DOF_INDEX(NODE,3)=0
          ENDIF
        ENDDO
C
C  WORKING IN ROW 3.
C
C  +--+--+--+
C  |31     |
C  +--+  +--+
C     |  |
C  +--+  +--+
C  |        |
C  +--+--+--+
C
        I=I+2*IREG_DENSITY_Y(2)-1
        DO ROW=1,2*IREG_DENSITY_Y(3)+1
          I=I+1
          NODE=NODE+1
       IF(COL.EQ.1)THEN
         NODE_DOF_INDEX(NODE,1)=0
         NODE_DOF_INDEX(NODE,2)=0
         IF(ROW.GE.(2*2*IREG_DENSITY_Y(3))/5+1
     &        .AND.
     &        ROW.LE.(3*2*IREG_DENSITY_Y(3))/5+1)THEN
            NODE_DOF_INDEX(NODE,1)=-1
         ENDIF
       ELSE
         IF(ROW.EQ.2*IREG_DENSITY_Y(3)+1)THEN
           NODE_DOF_INDEX(NODE,1)=0
           NODE_DOF_INDEX(NODE,2)=0
           IF(COL.GE.(10*2*IREG_DENSITY_X(1))/45+1
     &          .AND.
     &          COL.LE.(11*2*IREG_DENSITY_X(1))/45+1)THEN
              NODE_DOF_INDEX(NODE,2)=-5
           ENDIF
           IF(COL.GE.(22*2*IREG_DENSITY_X(1))/45+1
     &          .AND.
     &          COL.LE.(23*2*IREG_DENSITY_X(1))/45+1)THEN
              NODE_DOF_INDEX(NODE,2)=-5
           ENDIF
           IF(COL.GE.(34*2*IREG_DENSITY_X(1))/45+1
     &          .AND.
     &          COL.LE.(35*2*IREG_DENSITY_X(1))/45+1)THEN
              NODE_DOF_INDEX(NODE,2)=-5
           ENDIF
         ELSE
           IF(ROW.EQ.1)THEN
             NODE_DOF_INDEX(NODE,1)=0
             NODE_DOF_INDEX(NODE,2)=0
             IF(COL.GE.(10*2*IREG_DENSITY_X(1))/45+1
     &            .AND.
     &            COL.LE.(11*2*IREG_DENSITY_X(1))/45+1)THEN
                NODE_DOF_INDEX(NODE,2)=-1
             ENDIF
             IF(COL.GE.(22*2*IREG_DENSITY_X(1))/45+1
     &            .AND.
     &            COL.LE.(23*2*IREG_DENSITY_X(1))/45+1)THEN
                NODE_DOF_INDEX(NODE,2)=-1
             ENDIF
             IF(COL.GE.(34*2*IREG_DENSITY_X(1))/45+1
     &            .AND.
     &            COL.LE.(35*2*IREG_DENSITY_X(1))/45+1)THEN
                NODE_DOF_INDEX(NODE,2)=-1
             ENDIF
           ELSE
             DOF=DOF+1
             NODE_DOF_INDEX(NODE,1)=DOF
             DOF=DOF+1
             NODE_DOF_INDEX(NODE,2)=DOF
           ENDIF
         ENDIF
       ENDIF
          IF(MOD(J,2).EQ.1.AND.MOD(I,2).EQ.1)THEN
            DOF=DOF+1
            NODE_DOF_INDEX(NODE,3)=DOF
          ELSE
            NODE_DOF_INDEX(NODE,3)=0
          ENDIF
        ENDDO
      ENDDO
C
C  WORKING IN COLUMN 2,INCLUDING NODES ON EXTREME RIGHT AND LEFT.
C
      DO COL=1,2*IREG_DENSITY_X(2)+1
C
C  WORKING IN ROW 1.
C
C  +--+--+--+
C  |        |
C  +--+  +--+
C     |  |
C  +--+  +--+
C  |   12   |
C  +--+--+--+
C
        J=J+1
        I=0
        DO ROW=1,2*IREG_DENSITY_Y(1)+1
          I=I+1
          NODE=NODE+1
       IF(ROW.EQ.1)THEN
         NODE_DOF_INDEX(NODE,1)=0
         NODE_DOF_INDEX(NODE,2)=0
       ELSE
         IF(ROW.EQ.2*IREG_DENSITY_Y(1)+1
     &       .AND.
     &       COL.EQ.1)THEN
               NODE_DOF_INDEX(NODE,1)=0
               NODE_DOF_INDEX(NODE,2)=0
         ELSE
           IF(ROW.EQ.2*IREG_DENSITY_Y(1)+1
     &         .AND.
     &         COL.EQ.2*IREG_DENSITY_X(2)+1)THEN
                 NODE_DOF_INDEX(NODE,1)=0
                 NODE_DOF_INDEX(NODE,2)=0
           ELSE
             DOF=DOF+1
             NODE_DOF_INDEX(NODE,1)=DOF
             DOF=DOF+1
             NODE_DOF_INDEX(NODE,2)=DOF
           ENDIF
         ENDIF
       ENDIF
          IF(MOD(J,2).EQ.1.AND.MOD(I,2).EQ.1)THEN
            DOF=DOF+1
            NODE_DOF_INDEX(NODE,3)=DOF
          ELSE
            NODE_DOF_INDEX(NODE,3)=0
          ENDIF
        ENDDO
C
C  WORKING IN ROW 2.
C
C  +--+--+--+
C  |        |
C  +--+  +--+
C     |22|
C  +--+  +--+
C  |        |
C  +--+--+--+
C
        DO ROW=2,2*IREG_DENSITY_Y(2)
          I=I+1
          NODE=NODE+1
       IF(COL.EQ.1)THEN
          NODE_DOF_INDEX(NODE,1)=0
          NODE_DOF_INDEX(NODE,2)=0
       ELSE
         IF(COL.EQ.2*IREG_DENSITY_X(2)+1)THEN
            NODE_DOF_INDEX(NODE,1)=0
            NODE_DOF_INDEX(NODE,2)=0
         ELSE
             DOF=DOF+1
             NODE_DOF_INDEX(NODE,1)=DOF
             DOF=DOF+1
             NODE_DOF_INDEX(NODE,2)=DOF
         ENDIF
       ENDIF
          IF(MOD(J,2).EQ.1.AND.MOD(I,2).EQ.1)THEN
            DOF=DOF+1
            NODE_DOF_INDEX(NODE,3)=DOF
          ELSE
            NODE_DOF_INDEX(NODE,3)=0
          ENDIF
        ENDDO
C
C  WORKING IN ROW 3.
C
C  +--+--+--+
C  |   32   |
C  +--+  +--+
C     |  |
C  +--+  +--+
C  |        |
C  +--+--+--+
C
        DO ROW=1,2*IREG_DENSITY_Y(3)+1
          I=I+1
          NODE=NODE+1
       IF(ROW.EQ.1
     &     .AND.
     &     COL.EQ.1)THEN
          NODE_DOF_INDEX(NODE,1)=0
          NODE_DOF_INDEX(NODE,2)=0
       ELSE
         IF(ROW.EQ.1
     &       .AND.
     &       COL.EQ.2*IREG_DENSITY_X(2)+1)THEN
            NODE_DOF_INDEX(NODE,1)=0
            NODE_DOF_INDEX(NODE,2)=0
         ELSE
           IF(ROW.EQ.2*IREG_DENSITY_Y(3)+1)THEN
             NODE_DOF_INDEX(NODE,1)=0
             NODE_DOF_INDEX(NODE,2)=0
             IF(COL.GE.(5*2*IREG_DENSITY_X(2))/15
     &           .AND.
     &           COL.LE.(10*2*IREG_DENSITY_X(2))/15)THEN
               DOF=DOF+1
C
C  18 MAY 2004,FOLLOWING LINE CHANGED FOR MDG:
C
               NODE_DOF_INDEX(NODE,2)=DOF
             ENDIF
           ELSE
             DOF=DOF+1
             NODE_DOF_INDEX(NODE,1)=DOF
             DOF=DOF+1
             NODE_DOF_INDEX(NODE,2)=DOF
           ENDIF
         ENDIF
       ENDIF
          IF(MOD(J,2).EQ.1.AND.MOD(I,2).EQ.1)THEN
            DOF=DOF+1
            NODE_DOF_INDEX(NODE,3)=DOF
          ELSE
            NODE_DOF_INDEX(NODE,3)=0
          ENDIF
        ENDDO
      ENDDO
C
C  WORKING IN COLUMN 3,AND NOT HANDLING NODES IN EXTREME LEFT.
C
      DO COL=2,2*IREG_DENSITY_X(3)+1
C
C  WORKING IN ROW 1.
C
C  +--+--+--+
C  |        |
C  +--+  +--+
C     |  |
C  +--+  +--+
C  |      13|
C  +--+--+--+
C
        J=J+1
        I=0
        DO ROW=1,2*IREG_DENSITY_Y(1)+1
          I=I+1
          NODE=NODE+1
       IF(COL.EQ.2*IREG_DENSITY_X(3)+1)THEN
         NODE_DOF_INDEX(NODE,1)=0
         NODE_DOF_INDEX(NODE,2)=0
         IF(ROW.GE.(2*2*IREG_DENSITY_Y(1))/5+1
     &        .AND.
     &        ROW.LE.(3*2*IREG_DENSITY_Y(1))/5+1)THEN
            NODE_DOF_INDEX(NODE,1)=-4
         ENDIF
       ELSE
         IF(ROW.EQ.1)THEN
            NODE_DOF_INDEX(NODE,1)=0
            NODE_DOF_INDEX(NODE,2)=0
         ELSE
           IF(ROW.EQ.2*IREG_DENSITY_Y(1)+1)THEN
             NODE_DOF_INDEX(NODE,1)=0
             NODE_DOF_INDEX(NODE,2)=0
             IF(COL.GE.(10*2*IREG_DENSITY_X(3))/45+1
     &            .AND.
     &            COL.LE.(11*2*IREG_DENSITY_X(3))/45+1)THEN
                NODE_DOF_INDEX(NODE,2)=-4
             ENDIF
             IF(COL.GE.(22*2*IREG_DENSITY_X(3))/45+1
     &            .AND.
     &            COL.LE.(23*2*IREG_DENSITY_X(3))/45+1)THEN
                NODE_DOF_INDEX(NODE,2)=-4
             ENDIF
             IF(COL.GE.(34*2*IREG_DENSITY_X(3))/45+1
     &            .AND.
     &            COL.LE.(35*2*IREG_DENSITY_X(3))/45+1)THEN
                NODE_DOF_INDEX(NODE,2)=-4
             ENDIF
           ELSE
             DOF=DOF+1
             NODE_DOF_INDEX(NODE,1)=DOF
             DOF=DOF+1
             NODE_DOF_INDEX(NODE,2)=DOF
           ENDIF
         ENDIF
       ENDIF
          IF(MOD(J,2).EQ.1.AND.MOD(I,2).EQ.1)THEN
            DOF=DOF+1
            NODE_DOF_INDEX(NODE,3)=DOF
          ELSE
            NODE_DOF_INDEX(NODE,3)=0
          ENDIF
        ENDDO
C
C  WORKING IN ROW 3.
C
C  +--+--+--+
C  |      33|
C  +--+  +--+
C     |  |
C  +--+  +--+
C  |        |
C  +--+--+--+
C
        I=I+2*IREG_DENSITY_Y(2)-1
        DO ROW=1,2*IREG_DENSITY_Y(3)+1
          I=I+1
          NODE=NODE+1
       IF(COL.EQ.2*IREG_DENSITY_X(3)+1)THEN
         NODE_DOF_INDEX(NODE,1)=0
         NODE_DOF_INDEX(NODE,2)=0
         IF(ROW.GE.(2*2*IREG_DENSITY_Y(3))/5+1
     &        .AND.
     &        ROW.LE.(3*2*IREG_DENSITY_Y(3))/5+1)THEN
            NODE_DOF_INDEX(NODE,1)=-2
         ENDIF
       ELSE
         IF(ROW.EQ.2*IREG_DENSITY_Y(3)+1)THEN
           NODE_DOF_INDEX(NODE,1)=0
           NODE_DOF_INDEX(NODE,2)=0
           IF(COL.GE.(10*2*IREG_DENSITY_X(3))/45+1
     &          .AND.
     &          COL.LE.(11*2*IREG_DENSITY_X(3))/45+1)THEN
              NODE_DOF_INDEX(NODE,2)=-6
           ENDIF
           IF(COL.GE.(22*2*IREG_DENSITY_X(3))/45+1
     &          .AND.
     &          COL.LE.(23*2*IREG_DENSITY_X(3))/45+1)THEN
              NODE_DOF_INDEX(NODE,2)=-6
           ENDIF
           IF(COL.GE.(34*2*IREG_DENSITY_X(3))/45+1
     &          .AND.
     &          COL.LE.(35*2*IREG_DENSITY_X(3))/45+1)THEN
              NODE_DOF_INDEX(NODE,2)=-6
           ENDIF
         ELSE
           IF(ROW.EQ.1)THEN
             NODE_DOF_INDEX(NODE,1)=0
             NODE_DOF_INDEX(NODE,2)=0
             IF(COL.GE.(10*2*IREG_DENSITY_X(3))/45+1
     &            .AND.
     &            COL.LE.(11*2*IREG_DENSITY_X(3))/45+1)THEN
                NODE_DOF_INDEX(NODE,2)=-2
             ENDIF
             IF(COL.GE.(22*2*IREG_DENSITY_X(3))/45+1
     &            .AND.
     &            COL.LE.(23*2*IREG_DENSITY_X(3))/45+1)THEN
                NODE_DOF_INDEX(NODE,2)=-2
             ENDIF
             IF(COL.GE.(34*2*IREG_DENSITY_X(3))/45+1
     &            .AND.
     &            COL.LE.(35*2*IREG_DENSITY_X(3))/45+1)THEN
                NODE_DOF_INDEX(NODE,2)=-2
             ENDIF
           ELSE
             DOF=DOF+1
             NODE_DOF_INDEX(NODE,1)=DOF
             DOF=DOF+1
             NODE_DOF_INDEX(NODE,2)=DOF
           ENDIF
         ENDIF
       ENDIF
          IF(MOD(J,2).EQ.1.AND.MOD(I,2).EQ.1)THEN
            DOF=DOF+1
            NODE_DOF_INDEX(NODE,3)=DOF
          ELSE
            NODE_DOF_INDEX(NODE,3)=0
          ENDIF
        ENDDO
      ENDDO
      RETURN
      END
      SUBROUTINE QUAD_A_SET(MAXEL,NELEMN,NQUAD,AREA,XM,YM)
C
C***********************************************************************
C
C  QUAD_A_SET SETS QUADRATURE INFORMATION FOR THE ASSEMBLY ROUTINE.
C
C  AUTHOR:
C
C    JOHN BURKARDT
C
C  PARAMETERS:
C
C    INPUT, MAXEL,THE MAXIMUM NUMBER OF ELEMENTS.
C    INPUT, NELEMN,THE NUMBER OF ELEMENTS.
C    INPUT, NQUAD,THE NUMBER OF QUADRATURE POINTS PER ELEMENT.
C    OUTPUT,DOUBLE PRECISION AREA(NELEMN),THE NOMINAL AREA OF EACH ELEMENT.
C    OUTPUT,DOUBLE PRECISION XM(MAXEL,3),YM(MAXEL,3),THE COORDINATES OF
C    THE QUADRATURE POINTS.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION AREA(MAXEL)
      DIMENSION XM(MAXEL,NQUAD)
      DIMENSION YM(MAXEL,NQUAD)
      INTEGER ELEMENT
C
      DO ELEMENT=1,NELEMN
        XM(ELEMENT,1)=0.5D0
        XM(ELEMENT,2)=0.5D0
        XM(ELEMENT,3)=0.0D0
      ENDDO
      DO ELEMENT=1,NELEMN
        YM(ELEMENT,1)=0.0D0
        YM(ELEMENT,2)=0.5D0
        YM(ELEMENT,3)=0.5D0
      ENDDO
      DO ELEMENT=1,NELEMN
        AREA(ELEMENT)=0.5D0
      ENDDO
      RETURN
      END
      SUBROUTINE ELEMENT_NODE_BANDWIDTH(MAXND,MAXEL,NNODES,
     &   NELEMN,NODE,NEQN,NP,INDX,NLBAND,NUBAND,NBAND)
C
C***********************************************************************
C
C  ELEMENT_NODE_BANDWIDTH DETERMINES THE BANDWIDTH ASSOCIATED WITH THE GRID.
C
C  AUTHOR:
C
C    JOHN BURKARDT
C
C  PARAMETERS:
C
C    INPUT, MAXND,THE MAXIMUM NUMBER OF NODES
C    INPUT, MAXEL,THE MAXIMUM NUMBER OF ELEMENTS.
C    INPUT, NNODES,THE ORDER OF THE ELEMENTS.
C    INPUT, NELEMN,THE NUMBER OF ELEMENTS.
C    INPUT, NODE(MAXEL,NNODES),THE NODES IN EACH ELEMENT.
C    INPUT, NEQN,THE NUMBER OF DEGREES OF FREEDOM.
C    INPUT, NP,THE NUMBER OF NODES.
C    INPUT, INDX(MAXND,3),THE NODAL DEGREES OF FREEDOM.
C    OUTPUT, NLBAND,NUBAND,NBAND,THE LOWER,UPPER AND TOTAL BANDWIDTHS.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION INDX(MAXND,3)
      DIMENSION NODE(MAXEL,NNODES)
      INTEGER DOF_MAX(NEQN)
      INTEGER DOF_MIN(NEQN)
      INTEGER ELEMENT
      INTEGER P1
      INTEGER P2
      INTEGER U1
      INTEGER U2
      INTEGER V1
      INTEGER V2
C
      DO IEQN=1,NEQN
        DOF_MIN(IEQN)=IEQN
        DOF_MAX(IEQN)=IEQN
      ENDDO
      DO ELEMENT=1,NELEMN
        DO L1=1,NNODES
          N1=NODE(ELEMENT,L1)
          U1=INDX(N1,1)
          V1=INDX(N1,2)
          P1=INDX(N1,3)
          DO L2=1,NNODES
            N2=NODE(ELEMENT,L2)
            U2=INDX(N2,1)
            V2=INDX(N2,2)
            P2=INDX(N2,3)
            IF(1.LE.U1.AND.1.LE.U2)THEN
              DOF_MIN(U1)=MIN(DOF_MIN(U1),U2)
              DOF_MAX(U1)=MAX(DOF_MAX(U1),U2)
            ENDIF
            IF(1.LE.U1.AND.1.LE.V2)THEN
              DOF_MIN(U1)=MIN(DOF_MIN(U1),V2)
              DOF_MAX(U1)=MAX(DOF_MAX(U1),V2)
            ENDIF
            IF(1.LE.U1.AND.1.LE.P2)THEN
              DOF_MIN(U1)=MIN(DOF_MIN(U1),P2)
              DOF_MAX(U1)=MAX(DOF_MAX(U1),P2)
            ENDIF
            IF(1.LE.V1.AND.1.LE.U2)THEN
              DOF_MIN(V1)=MIN(DOF_MIN(V1),U2)
              DOF_MAX(V1)=MAX(DOF_MAX(V1),U2)
            ENDIF
            IF(1.LE.V1.AND.1.LE.V2)THEN
              DOF_MIN(V1)=MIN(DOF_MIN(V1),V2)
              DOF_MAX(V1)=MAX(DOF_MAX(V1),V2)
            ENDIF
            IF(1.LE.V1.AND.1.LE.P2)THEN
              DOF_MIN(V1)=MIN(DOF_MIN(V1),P2)
              DOF_MAX(V1)=MAX(DOF_MAX(V1),P2)
            ENDIF
            IF(1.LE.P1.AND.1.LE.U2)THEN
              DOF_MIN(P1)=MIN(DOF_MIN(P1),U2)
              DOF_MAX(P1)=MAX(DOF_MAX(P1),U2)
            ENDIF
            IF(1.LE.P1.AND.1.LE.V2)THEN
              DOF_MIN(P1)=MIN(DOF_MIN(P1),V2)
              DOF_MAX(P1)=MAX(DOF_MAX(P1),V2)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      NLBAND=0
      NUBAND=0
      DO IEQN=1,NEQN
        NLBAND=MAX(NLBAND,IEQN-DOF_MIN(IEQN))
        NUBAND=MAX(NUBAND,DOF_MAX(IEQN)-IEQN)
      ENDDO
      NBAND=NLBAND+NUBAND+1
      WRITE(*,'(A)')'ELEMENT_NODE_BANDWIDTH:'
      WRITE(*,'(A,I6)')'LOWER HALF BANDWIDTH=',NLBAND
      WRITE(*,'(A,I6)')'UPPER HALF BANDWIDTH=',NUBAND
      WRITE(*,'(A,I6)')'TOTAL BANDWIDTH=     ',NBAND
      RETURN
      END
      SUBROUTINE NSTOKE(XC,YC,AREA,XM,YM,
     &   A,F,G,UOLD,REYNLD,TOLNS,XLNGTH,YLNGTH,
     &   NODE,INDX,IPIVOT,MROW1,
     &   NLBAND,NUBAND,NBAND,NROW1,NCOL1,
     &   NELEMN,NP,NNODES,NUK,NQUAD,NEQN,
     &   NSTEPS,NSIM,MAXND,MAXEL,RDEL,ALPHA)
C
C**********************************************************************
C
C  NSTOKE SOLVES THE NAVIER-STOKES EQUATIONS USING TAYLOR-HOOD ELEMENTS.
C
C  AUTHOR:
C
C    HYUNG-CHUN LEE,
C    DEPARTMENT OF MATHEMATICS,
C    AJOU UNIVERSITY,KOREA
C
C  PARAMETERS:
C
C    INPUT, NSIM,THE NUMBER OF SIMPLE ITERATIONS PERFORMED.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION A(MROW1,*)
      DIMENSION AREA(*)
      DIMENSION F(*)
      DIMENSION G(*)
      DIMENSION INDX(MAXND,3)
      DIMENSION IPIVOT(*)
      DIMENSION NODE(MAXEL,*)
      DIMENSION UN(2)
      DIMENSION UNX(2)
      DIMENSION UNY(2)
      DIMENSION UOLD(*)
      DIMENSION XC(*)
      DIMENSION XM(MAXEL,*)
      DIMENSION YC(*)
      DIMENSION YM(MAXEL,*)
      INTEGER UNK_U
      INTEGER UNK_V
C
      VISC=1.D0/REYNLD
C
C  MATRIX ASSEMBLY TRIANGLE BY TRIANGLE
C
      DO ITER=1,NSTEPS
        NITER=ITER
        IF(ITER.LT.NSIM)THEN
          CSIM=0.0D0
        ELSE
          CSIM=1.0D0
        ENDIF
        DO I=1,NROW1
          DO J=1,NCOL1
            A(I,J)=0.D0
          ENDDO
        ENDDO
        DO IT=1,NELEMN
          ARR=AREA(IT)/3.D0
          DO IQUAD=1,NQUAD
            Y=YM(IT,IQUAD)
            X=XM(IT,IQUAD)
            CALL TRANS(IT,X,Y,DET,XIX,XIY,ETAX,ETAY,XC,YC,NODE,MAXEL)
            AR=ARR*DET
            DO KK=1,2
              UN(KK)=0.D0
              UNY(KK)=0.D0
              UNX(KK)=0.D0
            ENDDO
            UOLD_QP=0.0D0
            VOLD_QP=0.0D0
            DO IQ=1,NNODES
              CALL REFQBF(X,Y,IQ,BB,TBX,TBY)
              BX=TBX*XIX+TBY*ETAX
              BY=TBX*XIY+TBY*ETAY
              IP=NODE(IT,IQ)
              DO IUK=1,2
                IUN=INDX(IP,IUK)
                IF(0.LT.IUN)THEN
                  UN(IUK)=UN(IUK)+BB*G(IUN)
                  UNX(IUK)=UNX(IUK)+BX*G(IUN)
                  UNY(IUK)=UNY(IUK)+BY*G(IUN)
                  IF(IUK.EQ.1)UOLD_QP=UOLD_QP+BB*UOLD(IUN)
                  IF(IUK.EQ.2)VOLD_QP=VOLD_QP+BB*UOLD(IUN)
                ELSE IF(IUN.LT.0)THEN
                  UBC=ALPHA*UBDRY(1,IUK,IP,XC,YC)
                  UN(IUK)=UN(IUK)+BB*UBC
                  UNX(IUK)=UNX(IUK)+BX*UBC
                  UNY(IUK)=UNY(IUK)+BY*UBC
                  IF(IUK.EQ.1)UOLD_QP=UOLD_QP+BB*UBC
                  IF(IUK.EQ.2)VOLD_QP=VOLD_QP+BB*UBC
                ENDIF
              ENDDO
            ENDDO
            DO IQ=1,NNODES
              IP=NODE(IT,IQ)
              CALL REFQBF(X,Y,IQ,BB,TBX,TBY)
              BX=TBX*XIX+TBY*ETAX
              BY=TBX*XIY+TBY*ETAY
              IF(IQ.LE.3)THEN
                BBL=REFBSP(X,Y,IQ)
              ENDIF
              DO 210 IUK=1,NUK
                I=INDX(IP,IUK)
                IF(I.LE.0)GOTO 210
                IF(IUK.EQ.1)THEN
                  F(I)=F(I)+CSIM *
     &              ((UN(1)*UNX(1)+UN(2)*UNY(1))*BB)*AR
     &               +RDEL*UOLD_QP*BB*AR
                ELSE IF(IUK.EQ.2)THEN
                  F(I)=F(I)+CSIM *
     &              ((UN(1)*UNX(2)+UN(2)*UNY(2))*BB)*AR
     &               +RDEL*VOLD_QP*BB*AR
                ENDIF
                DO IQQ=1,NNODES
                  IPP=NODE(IT,IQQ)
                  CALL REFQBF(X,Y,IQQ,BBB,TBBX,TBBY)
                  BBX=TBBX*XIX+TBBY*ETAX
                  BBY=TBBX*XIY+TBBY*ETAY
                  IF(IQQ.LE.3)THEN
                    BBBL=REFBSP(X,Y,IQQ)
                  ENDIF
                  DO 190 IUKK=1,NUK
                    J=INDX(IPP,IUKK)
                    IF(J.EQ.0)GOTO 190
                    AIJ=0.D0
                    IF(I.EQ.NEQN)GOTO 190
                    IF(IUK.EQ.1)THEN
                      IF(IUKK.EQ.1)THEN
                        AIJ=VISC*(BY*BBY+BX*BBX)
     &                   +(BBB*UNX(1)*BB)*CSIM
     &                   +BB*BBX*UN(1)
     &                   +BB*BBY*UN(2)+RDEL*(BB*BBB)
                      ELSE IF(IUKK.EQ.2)THEN
                        AIJ=CSIM*(BB*BBB*UNY(1))
                      ELSE IF(IUKK.EQ.3)THEN
                        AIJ=-BX*BBBL
                      ENDIF
                    ELSE IF(IUK.EQ.2)THEN
                      IF(IUKK.EQ.1)THEN
                        AIJ= CSIM*(BB*BBB*UNX(2))
                      ELSE IF(IUKK.EQ.2)THEN
                        AIJ=(VISC*(BY*BBY+BX*BBX)
     &                     +(BB*BBB*UNY(2))*CSIM
     &                     +BB*BBY*UN(2)
     &                     +BB*BBX*UN(1))+RDEL*(BB*BBB)
                      ELSE IF(IUKK.EQ.3)THEN
                        AIJ=-BY*BBBL
                      ENDIF
                    ELSE IF(IUK.EQ.3)THEN
                      IF(IUKK.EQ.1)THEN
                        AIJ=BBX*BBL
                      ELSE IF(IUKK.EQ.2)THEN
                        AIJ=BBY*BBL
                      ENDIF
                    ENDIF
C
C  THE COEFFICIENT AIJ IS ADDED TO THE MATRIX ENTRY IF J REPRESENTS AN UNKNOWN,
C  OR IS SUBTRACTED FROM THE RIGHT HAND SIDE IF J CORRESPONDS TO A KNOWN "VARIABLE".
C
                  IF(0.LT.J)THEN
                    IUSE=I-J+NBAND
                    A(IUSE,J)=A(IUSE,J)+AIJ*AR
                  ELSE
                    F(I)=F(I)-AR*ALPHA*UBDRY(IUKK,IPP,J,XC,YC)*AIJ
                  ENDIF
 190              CONTINUE
                ENDDO
 210          CONTINUE
            ENDDO
          ENDDO
        ENDDO
C
C  REPLACE LAST EQUATION TO SET A REFERENCE PRESSURE TO 0.
C
        F(NEQN)=0.D0
        DO J=NEQN-NLBAND,NEQN-1
          I=NEQN-J+NBAND
          A(I,J)=0.D0
        ENDDO
        A(NBAND,NEQN)=1.0D0
C
C  SOLVE THE SYSTEM.
C
        JOB=0
        CALL DGBFA(A,MROW1,NEQN,NLBAND,NUBAND,IPIVOT,INFO)
        CALL DGBSL(A,MROW1,NEQN,NLBAND,NUBAND,IPIVOT,F,JOB)
C
C  CHECK FOR CONVERGENCE
C
        DIFF=0.D0
        DO I=1,NP
          UNK_U=INDX(I,1)
          IF(0.LT.UNK_U)THEN
            DIFF=DIFF+( G(UNK_U)-F(UNK_U))**2
          ENDIF
          UNK_V=INDX(I,2)
          IF(0.LT.UNK_V)THEN
            DIFF=DIFF+( G(UNK_V)-F(UNK_V))**2
          ENDIF
        ENDDO
        DIFF=SQRT(DIFF)
        WRITE(*,*)'ITERATION ',ITER,'  DIFFERENCE IS ',DIFF
        IF(DIFF.LE.TOLNS)THEN
          RETURN
        ENDIF
        DO I=1,NEQN
          G(I)=F(I)
          F(I)=0.D0
        ENDDO
      ENDDO
      RETURN
      END
      FUNCTION UBDRY(IUK,IP,IUN,XC,YC)
C
C***********************************************************************
C
C  UBDRY EVALUATES BOUNDARY CONDITIONS AT A NODE.
C
C  AUTHOR:
C
C    MAX GUNZBURGER
C
C  PARAMETERS:
C
C    INPUT, IUK,INDICATES THE CLASS OF VARIABLE BEING SET.
C   *1=HORIZONTAL VELOCITY;
C   *2=VERTICAL VELOCITY;
C   *3=PRESSURE.
C    INPUT, IP,THE INDEX OF THE NODE.
C    INPUT, IUN,THE VALUE OF INDX(IP,IUK),THAT IS,THE GLOBAL
C    UNKNOWN NUMBER THAT WAS ASSIGNED TO THIS "DEGREE OF FREEDOM".
C    INPUT,DOUBLE PRECISION XC(*),YC(*),THE COORDINATES OF NODES.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XC(*),YC(*)
      UBDRY=0D0
      X=XC(IP)
      Y=YC(IP)
      IF(IUK.EQ.1)THEN
        IF(Y.LE.3D0)THEN
          IF(X.LT.5D0)THEN
            UBDRY=1D0
          ELSE
            UBDRY=-1D0
          ENDIF
        ELSEIF(Y.GE.6D0)THEN
          IF(X.LT.5D0)THEN
            UBDRY=-1D0
          ELSE
            UBDRY=1D0
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END
      FUNCTION REFBSP(X,Y,IQ)
C***********************************************************************
C
C  REFBSP EVALUATES A LINEAR BASIS FUNCTIONS ON THE REFERENCE TRIANGLE.
C
C  AUTHOR:
C
C    HYUNG-CHUN LEE,
C    DEPARTMENT OF MATHEMATICS,
C    AJOU UNIVERSITY,KOREA
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IF(IQ.EQ.1)THEN
        REFBSP=1.0D0-X-Y
      ELSE IF(IQ.EQ.2)THEN
        REFBSP=X
      ELSE IF(IQ.EQ.3)THEN
        REFBSP=Y
      ELSE
        STOP'REFBSP-FATAL ERROR'
      ENDIF
      RETURN
      END
      SUBROUTINE REFQBF(X,Y,IN,BB,BX,BY)
C************************************************************************
C
C  REFQBF EVALUATES QUADRATIC BASIS FUNCTIONS ON THE REFERENCE TRIANGLE.
C
C  AUTHOR:
C
C    HYUNG-CHUN LEE,
C    DEPARTMENT OF MATHEMATICS,
C    AJOU UNIVERSITY,KOREA
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      IF(IN.EQ.1)THEN
        BB=(1.D0-X-Y)*(1.D0-2.D0*X-2.D0*Y)
        BX=-3.D0+4.D0*X+4.D0*Y
        BY=-3.D0+4.D0*X+4.D0*Y
      ELSE IF(IN.EQ.2)THEN
        BB=X*(2.D0*X-1.D0)
        BX=4.D0*X-1.D0
        BY=0.D0
      ELSE IF(IN.EQ.3)THEN
        BB=Y*(2.D0*Y-1.D0)
        BX=0.D0
        BY=4.D0*Y-1.D0
      ELSE IF(IN.EQ.4)THEN
        BB=4.D0*X*(1.D0-X-Y)
        BX=4.D0*(1.D0-2.D0*X-Y)
        BY=-4.D0*X
      ELSE IF(IN.EQ.5)THEN
        BB=4.D0*X*Y
        BX=4.D0*Y
        BY=4.D0*X
      ELSE IF(IN.EQ.6)THEN
        BB=4.D0*Y*(1.D0-X-Y)
        BX=-4.D0*Y
        BY=4.D0*(1.D0-X-2.D0*Y)
      ELSE
        BB=0.0D0
        BX=0.0D0
        BY=0.0D0
      ENDIF
      RETURN
      END
      SUBROUTINE TRANS(IT,XQ,YQ,DET,PJ11,PJ21,PJ12,PJ22,
     &  XC,YC,NODE,MAXEL)
C
C ***********************************************************************
C
C  TRANS TRANSFORMS DATA BETWEEN THE REFERENCE AND PHYSICAL ELEMENTS.
C
C  AUTHOR:
C
C    HYUNG-CHUN LEE,
C    DEPARTMENT OF MATHEMATICS,
C    AJOU UNIVERSITY,KOREA
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XC(*)
      DIMENSION YC(*)
      DIMENSION NODE(MAXEL,*)
C
      I1=NODE(IT,1)
      I2=NODE(IT,2)
      I3=NODE(IT,3)
      I4=NODE(IT,4)
      I5=NODE(IT,5)
      I6=NODE(IT,6)
      X1=XC(I1)
      Y1=YC(I1)
      X2=XC(I2)
      Y2=YC(I2)
      X3=XC(I3)
      Y3=YC(I3)
      X4=XC(I4)
      Y4=YC(I4)
      X5=XC(I5)
      Y5=YC(I5)
      X6=XC(I6)
      Y6=YC(I6)
C
C  COMPUTE PARTIAL DERIVATIVES AT POINT (XQ,YQ)
C
      F1X=X1*(-3.D0+4.D0*XQ+4.D0*YQ)
     &       +X2*(4.D0*XQ-1.D0)
     &       +X4*4.D0*(1.D0-2.D0*XQ-YQ)
     &       +X5*4.D0*YQ+X6*4.D0*(-YQ)
      F1Y=X1*(-3.D0+4.D0*XQ+4.D0*YQ)
     &       +X3*(4.D0*YQ-1.D0)
     &       +X4*4.D0*(-XQ)+X5*4.D0*XQ
     &       +X6*4.D0*(1.D0-XQ-2.D0*YQ)
      F2X=Y1*(-3.D0+4.D0*XQ+4.D0*YQ)
     &       +Y2*(4.D0*XQ-1.D0)
     &       +Y4*4.D0*(1.D0-2.D0*XQ-YQ)
     &       +Y5*4.D0*YQ+Y6*4.D0*(-YQ)
      F2Y=Y1*(-3.D0+4.D0*XQ+4.D0*YQ)
     &       +Y3*(4.D0*YQ-1.D0)
     &       +Y4*4.D0*(-XQ)+Y5*4.D0*XQ
     &       +Y6*4.D0*(1.D0-XQ-2.D0*YQ)
C
C  COMPUTE DETERMINANT OF TRANSFORMATION EVALUATED AT POINT (XQ,YQ)
C
      DET=F1X*F2Y-F1Y*F2X
C
C  COMPUTE J11,J22,J21,J22
C
      PJ11=F2Y/DET
      PJ12=-F2X/DET
      PJ21=-F1Y/DET
      PJ22=F1X/DET
      DET=DABS(DET)
      RETURN
      END
      SUBROUTINE DGBFA(ABD,LDA,N,ML,MU,IPVT,INFO)
C
C***********************************************************************
C
C  DGBFA FACTORS A DOUBLE PRECISION BAND MATRIX BY ELIMINATION.
C
C     DGBFA IS USUALLY CALLED BY DGBCO,BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C
C     ON ENTRY
C
C        ABD     DOUBLE PRECISION(LDA,N)
C                CONTAINS THE MATRIX IN BAND STORAGE.  THE COLUMNS
C                OF THE MATRIX ARE STORED IN THE COLUMNS OF  ABD  AND
C                THE DIAGONALS OF THE MATRIX ARE STORED IN ROWS
C                ML+1 THROUGH 2*ML+MU+1 OF  ABD .
C                SEE THE COMMENTS BELOW FOR DETAILS.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C                LDA MUST BE.GE.2*ML+MU+1 .
C
C        N       THE ORDER OF THE ORIGINAL MATRIX.
C
C        ML      NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
C                0.LE.ML.LT.N .
C
C        MU      NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C                0.LE.MU.LT.N .
C                MORE EFFICIENT IF  ML.LE.MU .
C     ON RETURN
C
C        ABD     AN UPPER TRIANGULAR MATRIX IN BAND STORAGE AND
C                THE MULTIPLIERS WHICH WERE USED TO OBTAIN IT.
C                THE FACTORIZATION CAN BE WRITTEN  A=L*U  WHERE
C                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
C                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
C
C        IPVT    INTEGER(N)
C                AN INTEGER VECTOR OF PIVOT INDICES.
C
C        INFO    INTEGER
C             =0  NORMAL VALUE.
C             =K  IF  U(K,K).EQ. 0.0 .  THIS IS NOT AN ERROR
C                     CONDITION FOR THIS SUBROUTINE,BUT IT DOES
C                     INDICATE THAT DGBSL WILL DIVIDE BY ZERO IF
C                     CALLED.  USE  RCOND  IN DGBCO FOR A RELIABLE
C                     INDICATION OF SINGULARITY.
C
C     BAND STORAGE
C
C           IF  A  IS A BAND MATRIX,THE FOLLOWING PROGRAM SEGMENT
C           WILL SET UP THE INPUT.
C
C                   ML=(BAND WIDTH BELOW THE DIAGONAL)
C                   MU=(BAND WIDTH ABOVE THE DIAGONAL)
C                   M=ML+MU+1
C                   DO J=1,N
C                      I1=MAX(1,J-MU)
C                      I2=MIN(N,J+ML)
C                      DO I=I1,I2
C                         K=I-J+M
C                         ABD(K,J)=A(I,J)
C                      ENDDO
C                   ENDDO
C
C           THIS USES ROWS  ML+1  THROUGH  2*ML+MU+1  OF  ABD .
C           IN ADDITION,THE FIRST  ML  ROWS IN  ABD  ARE USED FOR
C           ELEMENTS GENERATED DURING THE TRIANGULARIZATION.
C           THE TOTAL NUMBER OF ROWS NEEDED IN  ABD  IS  2*ML+MU+1 .
C           THE  ML+MU BY ML+MU  UPPER LEFT TRIANGLE AND THE
C           ML BY ML  LOWER RIGHT TRIANGLE ARE NOT REFERENCED.
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER,UNIVERSITY OF NEW MEXICO,ARGONNE NATIONAL LAB.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ABD(LDA,1),IPVT(1)
C
      M=ML+MU+1
      INFO=0
C
C     ZERO INITIAL FILL-IN COLUMNS
C
      J0=MU+2
      J1=MIN(N,M)-1
      DO JZ=J0,J1
        I0=M+1-JZ
        DO I=I0,ML
          ABD(I,JZ)=0.0D0
        ENDDO
      ENDDO
      JZ=J1
      JU=0
C
C     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
C
      NM1=N-1
      DO 120 K=1,NM1
         KP1=K+1
C
C        ZERO NEXT FILL-IN COLUMN
C
         JZ=JZ+1
         IF(JZ.GT.N)GOTO 50
           DO I=1,ML
             ABD(I,JZ)=0.0D0
           ENDDO
   50    CONTINUE
C
C        FIND L=PIVOT INDEX
C
         LM=MIN(ML,N-K)
         L=IDAMAX(LM+1,ABD(M,K),1)+M-1
         IPVT(K)=L+K-M
C
C        ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
C
         IF(ABD(L,K).EQ.0.0D0)GOTO 100
C
C           INTERCHANGE IF NECESSARY
C
            IF(L.NE.M)THEN
               T=ABD(L,K)
               ABD(L,K)=ABD(M,K)
               ABD(M,K)=T
            ENDIF
C
C  COMPUTE MULTIPLIERS
C
            T=-1.0D0/ABD(M,K)
            CALL DSCAL(LM,T,ABD(M+1,K),1)
C
C  ROW ELIMINATION WITH COLUMN INDEXING
C
            JU=MIN(MAX (JU,MU+IPVT(K)),N)
            MM=M
            DO J=KP1,JU
               L=L-1
               MM=MM-1
               T=ABD(L,J)
               IF(L.NE.MM)THEN
                  ABD(L,J)=ABD(MM,J)
                  ABD(MM,J)=T
               ENDIF
               CALL DAXPY(LM,T,ABD(M+1,K),1,ABD(MM+1,J),1)
            ENDDO
         GO TO 110
  100    CONTINUE
            INFO=K
  110    CONTINUE
  120 CONTINUE
  130 CONTINUE
      IPVT(N)=N
      IF(ABD(M,N).EQ.0.0D0)INFO=N
      RETURN
      END
      SUBROUTINE DGBSL(ABD,LDA,N,ML,MU,IPVT,B,JOB)
C
C***********************************************************************
C
C  DGBSL SOLVES A DOUBLE PRECISION BAND LINEAR SYSTEM.
C
C    THE LINEAR SYSTEM HAS THE FORM
C
C     A*X=B  OR  TRANS(A)*X=B
C
C     THE MATRIX HAS BEEN FACTORED BY DGBCO OR DGBFA.
C
C     ON ENTRY
C
C        ABD     DOUBLE PRECISION(LDA,N)
C                THE OUTPUT FROM DGBCO OR DGBFA.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  ABD .
C
C        N       THE ORDER OF THE ORIGINAL MATRIX.
C
C        ML      NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.
C
C        MU      NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DGBCO OR DGBFA.
C
C        B       DOUBLE PRECISION(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C             =0         TO SOLVE  A*X=B ,
C             =NONZERO   TO SOLVE  TRANS(A)*X=B ,WHERE
C                            TRANS(A) IS THE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF DGBCO HAS SET RCOND.GT.0.0
C        OR DGBFA HAS SET INFO.EQ.0 .
C
C     TO COMPUTE  INVERSE(A)*C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL DGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
C           IF(RCONDISTOOSMALL)GOTO ...
C           DO J=1,P
C              CALL DGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
C           ENDDO
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER,UNIVERSITY OF NEW MEXICO,ARGONNE NATIONAL LAB.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ABD(LDA,*),B(*),IPVT(*)
C
      M=MU+ML+1
      NM1=N-1
      IF(JOB.NE.0)GOTO 50
C
C        JOB=0 ,SOLVE  A*X=B
C        FIRST SOLVE L*Y=B
C
         IF(ML.NE.0)THEN
            DO K=1,NM1
               LM=MIN(ML,N-K)
               L=IPVT(K)
               T=B(L)
               IF(L.NE.K)THEN
                  B(L)=B(K)
                  B(K)=T
               ENDIF
               CALL DAXPY(LM,T,ABD(M+1,K),1,B(K+1),1)
            ENDDO
         ENDIF
C
C        NOW SOLVE  U*X=Y
C
         DO KB=1,N
            K=N+1-KB
            B(K)=B(K)/ABD(M,K)
            LM=MIN(K,M)-1
            LA=M-LM
            LB=K-LM
            T=-B(K)
            CALL DAXPY(LM,T,ABD(LA,K),1,B(LB),1)
         ENDDO
      GO TO 100
   50 CONTINUE
C
C        JOB=NONZERO,SOLVE  TRANS(A)*X=B
C        FIRST SOLVE  TRANS(U)*Y=B
C
         DO K=1,N
            LM=MIN(K,M)-1
            LA=M-LM
            LB=K-LM
            T=DDOT(LM,ABD(LA,K),1,B(LB),1)
            B(K)=(B(K)-T)/ABD(M,K)
         ENDDO
C
C        NOW SOLVE TRANS(L)*X=Y
C
         IF(ML.EQ.0)GOTO 90
            DO KB=1,NM1
               K=N-KB
               LM=MIN (ML,N-K)
               B(K)=B(K)+DDOT(LM,ABD(M+1,K),1,B(K+1),1)
               L=IPVT(K)
               IF(L.NE.K)THEN
                  T=B(L)
                  B(L)=B(K)
                  B(K)=T
               ENDIF
            ENDDO
   90    CONTINUE
  100 CONTINUE
      RETURN
      END
      SUBROUTINE DSCAL(N,DA,DX,INCX)
C
C***********************************************************************
C
C  DSCAL SCALES A VECTOR BY A CONSTANT.
C
C     USES UNROLLED LOOPS FOR INCREMENT EQUAL TO ONE.
C     JACK DONGARRA,LINPACK,3/11/78.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DX(*)
C
      IF(N.LE.0)THEN
        RETURN
      ENDIF
      IF(INCX.EQ.1)GOTO 20
C
C   CODE FOR INCREMENT NOT EQUAL TO 1
C
      NINCX=N*INCX
      DO I=1,NINCX,INCX
        DX(I)=DA*DX(I)
      ENDDO
      RETURN
C
C   CODE FOR INCREMENT EQUAL TO 1
C
C   CLEAN-UP LOOP
C
   20 M=MOD(N,5)
      DO I=1,M
        DX(I)=DA*DX(I)
      ENDDO
      DO I=M+1,N,5
        DX(I)=DA*DX(I)
        DX(I+1)=DA*DX(I+1)
        DX(I+2)=DA*DX(I+2)
        DX(I+3)=DA*DX(I+3)
        DX(I+4)=DA*DX(I+4)
      ENDDO
      RETURN
      END
      FUNCTION IDAMAX(N,DX,INCX)
C
C***********************************************************************
C
C  IDAMAX FINDS THE VECTOR ELEMENT OF LARGEST MAGNITUDE.
C
C     JACK DONGARRA,LINPACK,3/11/78.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DX(*)
C
      IDAMAX=0
      IF(N.LT.1)THEN
        RETURN
      ENDIF
      IDAMAX=1
      IF(N.EQ.1)THEN
        RETURN
      ENDIF
      IF(INCX.EQ.1)GOTO 20
C
C   CODE FOR INCREMENT NOT EQUAL TO 1
C
      IX=1
      DMAX=DABS(DX(1))
      IX=IX+INCX
      DO I=2,N
        IF(DMAX.LT.DABS(DX(IX)))THEN
          IDAMAX=I
          DMAX=DABS(DX(IX))
          IX=IX+INCX
        ENDIF
      ENDDO
      RETURN
C
C   CODE FOR INCREMENT EQUAL TO 1
C
   20 DMAX=DABS(DX(1))
      DO I=2,N
        IF(DMAX.LT.DABS(DX(I)))THEN
          IDAMAX=I
          DMAX=DABS(DX(I))
        ENDIF
      ENDDO
      RETURN
      END
      FUNCTION DDOT(N,DX,INCX,DY,INCY)
C
C***********************************************************************
C
C  DDOT FORMS THE DOT PRODUCT OF TWO VECTORS.
C
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA,LINPACK,3/11/78.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DX(*)
      DIMENSION DY(*)
C
      DDOT=0.0D0
      DTEMP=0.0D0
      IF(N.LE.0)THEN
        RETURN
      ENDIF
      IF(INCX.EQ.1.AND.INCY.EQ.1)GOTO 20
C
C    CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C    NOT EQUAL TO 1
C
      IX=1
      IY=1
      IF(INCX.LT.0)IX=(-N+1)*INCX+1
      IF(INCY.LT.0)IY=(-N+1)*INCY+1
      DO I=1,N
        DTEMP=DTEMP+DX(IX)*DY(IY)
        IX=IX+INCX
        IY=IY+INCY
      ENDDO
      DDOT=DTEMP
      RETURN
C
C    CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C    CLEAN-UP LOOP
C
   20 M=MOD(N,5)
      DO I=1,M
        DTEMP=DTEMP+DX(I)*DY(I)
      ENDDO
      IF(N.LT.5)GOTO 60
      DO I=M+1,N,5
      DTEMP=DTEMP+DX(I)*DY(I)+DX(I+1)*DY(I+1)+
     &  DX(I+2)*DY(I+2)+DX(I+3)*DY(I+3)+DX(I+4)*DY(I+4)
      ENDDO
   60 DDOT=DTEMP
      RETURN
      END
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
C
C***********************************************************************
C
C  DAXPY ADDS A MULTIPLE OF ONE VECTOR TO ANOTHER.
C
C  AUTHOR:
C
C    JACK DONGARRA
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DX(*)
      DIMENSION DY(*)
C
      IF(N.LE.0)THEN
        RETURN
      ENDIF
      IF(DA.EQ.0.0D0)THEN
        RETURN
      ENDIF
      IF(INCX.EQ.1.AND.INCY.EQ.1)GOTO 20
C
C    CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C    NOT EQUAL TO 1
C
      IX=1
      IY=1
      IF(INCX.LT.0)IX=(-N+1)*INCX+1
      IF(INCY.LT.0)IY=(-N+1)*INCY+1
      DO I=1,N
        DY(IY)=DY(IY)+DA*DX(IX)
        IX=IX+INCX
        IY=IY+INCY
      ENDDO
      RETURN
C
C   CODE FOR BOTH INCREMENTS EQUAL TO 1
C
   20 M=MOD(N,4)
      DO I=1,M
        DY(I)=DY(I)+DA*DX(I)
      ENDDO
      DO I=M+1,N,4
        DY(I)=DY(I) +DA*DX(I)
        DY(I+1)=DY(I+1)+DA*DX(I+1)
        DY(I+2)=DY(I+2)+DA*DX(I+2)
        DY(I+3)=DY(I+3)+DA*DX(I+3)
      ENDDO
      RETURN
      END
C     OPEN(1,FILE='hcell.plt',STATUS='UNKNOWN')
C     WRITE(1,*)'TITLE="HCELL DATA"'
C     WRITE(1,*)'VARIABLES="X","Y","U","V"'
C     WRITE(1,*)'ZONE N=',NODE_NUM,',E=',4*ELEMENT_NUM,
C    &  ',F=FEPOINT,ET=TRIANGLE'
C     DO NODE=1,NODE_NUM
C       I=INDX(NODE,1)
C       IF(0.LT.I)THEN
C         U=C(I)
C       ELSE IF(I.LT.0)THEN
C         U=UBDRY(1,NODE,I,X_NODE,Y_NODE)
C       ELSE
C         U=0.0D0
C       ENDIF
C       I=INDX(NODE,2)
C       IF(0.LT.I)THEN
C         V=C(I)
C       ELSE IF(I.LT.0)THEN
C         V=UBDRY(2,NODE,I,X_NODE,Y_NODE)
C       ELSE
C         V=0.0D0
C       ENDIF
C       WRITE(1,'(4G15.6)')X_NODE(NODE),Y_NODE(NODE),U,V
C     ENDDO
C     DO ELEMENT=1,ELEMENT_NUM
C       WRITE(1,'(3I6)')(ELEMENT_NODE(ELEMENT,J),J=1,3)
C     ENDDO
C     CLOSE(1)

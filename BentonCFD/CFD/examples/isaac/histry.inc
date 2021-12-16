C     Common that holds the residual history for this run.
C     This common restricts the run to have 10000 iterations.
C
C     R2ONE  : The residual on the very first iteration of the case
C
      PARAMETER (MXHIST = 40000)
      COMMON /IHSTRY/ ITTOT
      COMMON /HSTRY / R2(MXHIST), R2ONE

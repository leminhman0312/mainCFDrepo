SET BIN=C:\DFORT61\BIN
SET LIB=.;C:\DFORT61\LIB
SET PATH=C:\DFORT61\BIN;C:\UTILS;%windir%\SYSTEM32;%windir%
SET INCLUDE=.;C:\DFORT61\INCLUDE

f90 /nodebug nsc2ke.for
if errorlevel 1 goto error

:error

SET BIN=C:\DFORT61\BIN
SET LIB=C:\DFORT61\LIB
SET PATH=C:\DFORT61\BIN;C:\;C:\UTILS;%windir%\SYSTEM32;%windir%

f90 /nodebug unstruct2d.f90
if errorlevel 1 goto error

erase *.mod

:error

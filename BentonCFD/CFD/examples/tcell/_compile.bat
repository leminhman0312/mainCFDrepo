SET BIN=C:\VC32\BIN
SET LIB=.;C:\VC32\LIB
SET PATH=C:\VC32\BIN;C:\UTILS;%windir%\SYSTEM32;%windir%
SET INCLUDE=.;C:\VC32\INCLUDE

cl /Ox /W3 /WX tcell.c
if errorlevel 1 goto error

cl /Ox /W3 /WX hcell.c
if errorlevel 1 goto error

cl /Ox /W3 /WX fdm.c
if errorlevel 1 goto error

erase *.obj

:error

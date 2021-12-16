rename nast2d.c nast2d.c

SET BIN=C:\VC32\BIN
SET LIB=.;C:\VC32\LIB
SET PATH=C:\VC32\BIN;C:\UTILS;%windir%\SYSTEM32;%windir%
SET INCLUDE=.;C:\VC32\INCLUDE

cl /Ox /W3 /WX nast2d.c
if errorlevel 1 goto error

erase *.obj
erase 0?????.BMP

:error

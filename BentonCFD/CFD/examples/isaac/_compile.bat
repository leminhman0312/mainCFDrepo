rename isaac.for isaac.for
rename common.inc common.inc
rename histry.inc histry.inc

SET BIN=C:\DFORT61\BIN
SET LIB=C:\DFORT61\LIB
SET PATH=C:\DFORT61\BIN;C:\;C:\UTILS;%windir%\SYSTEM32;%windir%

f90 /nodebug isaac.for

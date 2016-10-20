@echo off
call setEnv
call %BIN_PATH%ExpIA.exe 1>cout.txt 2>cerr.txt
pause
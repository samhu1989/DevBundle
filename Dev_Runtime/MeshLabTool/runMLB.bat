@echo off
call setEnv
set ipath=%DATA_PATH%03\input\
set opath=%DATA_PATH%03\input_binary\
set spath=%BASE_PATH%mlx\fn2vn.mlx
@echo %ipath%
@echo %opath%
@echo %spath%
for /f "delims=\" %%i in ('dir /b "%ipath%*.ply"') do (
  @echo %%i
  call D:\MeshLab\meshlabserver.exe -i %ipath%%%i -s %spath% -o %opath%%%i -om vc vn
)
pause
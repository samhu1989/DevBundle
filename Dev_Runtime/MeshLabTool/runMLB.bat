@echo off
call setEnv
set ipath=%DATA_PATH%05\input_extracted\
set opath=%DATA_PATH%03\input_down\
set spath=%BASE_PATH%MeshLabTool\mlx\down.mlx
@echo %ipath%
@echo %opath%
@echo %spath%
for /f "delims=\" %%i in ('dir /b "%ipath%*.ply"') do (
  @echo %%i
  call D:\MeshLab\meshlabserver.exe -i %ipath%%%i -s %spath% -o %opath%%%i -om vc vn
)
pause
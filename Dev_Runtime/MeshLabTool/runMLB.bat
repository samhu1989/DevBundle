@echo off
call setEnv
set ipath=%DATA_PATH%09\aligned\
set opath=%DATA_PATH%09\down\
set spath=%BASE_PATH%MeshLabTool\mlx\sampling.mlx
@echo %ipath%
@echo %opath%
@echo %spath%
for /f "delims=\" %%i in ('dir /b "%ipath%*.ply"') do (
  @echo %%i
  call D:\MeshLab\meshlabserver.exe -i %ipath%%%i -s %spath% -o %opath%%%i -om vc vn
)
pause
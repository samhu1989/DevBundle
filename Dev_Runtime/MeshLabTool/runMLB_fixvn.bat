@echo off
call setEnv
set ipath1=%DATA_PATH%05\input_extracted_corrupted_vn\
set ipath2=%DATA_PATH%05\input_extracted_plane_mannual\
set opath=%DATA_PATH%05\input_extracted\
set spath=%BASE_PATH%mlx\fixvn.mlx
@echo %ipath1%
@echo %ipath2%
@echo %opath%
@echo %spath%
for /f "delims=\" %%i in ('dir /b "%ipath1%*.ply"') do (
  @echo %%i
  call D:\MeshLab\meshlabserver.exe -i %ipath1%%%i -i %ipath2%%%i -s %spath% -o %opath%%%i -om vc vn
)
pause	
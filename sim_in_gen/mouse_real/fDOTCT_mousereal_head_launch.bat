D:
cd D:\documents\istar\mcGPU_AlexWorkingCopy_VRT\build\Release
for /L %%A in (1,1,360) do mcGPU_vrt.exe in_mousereal\CTTB_study_%%A.in ..\vol\seg_digimouse.txt 100000000 out_mousereal\0%%A C:\Users\aortega\Documents\GitHub\sim_in_gen\mouse_real
C:
cd C:\Users\aortega\Documents\GitHub\sim_in_gen\mouse_real
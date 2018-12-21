function katsevich_bckMex(prjFile,outFile,nPrjs,prjSize,pxSize,volSize,volOr,vxSize,angSpan,startAng,P,R0,D)

%  
% katsevich_bckMex(prjFile,outFile,nPrjs,prjSize,pxSize,volSize,volOr,vxSize,angSpan,startAng,P,R0,D)
%
% Function to perform the backprojection stage of the Katsevich's reconstruction
% algorithm, according to the approaches by Noo, Wunderlich and Tan, and modified by Fontaine.
% 
% Refs: 
%       Noo et al., Phys Med Biol, 48, 3787, 2003. 
%       Wunderlich, Dept of Math Oregon Sta Univ., MS Thesis, 2006.
%       Tan et al., Med Phys 39, 4695, 2012.
%		Fontaine and Lee, ICPADS, 2007.
%
% Input:
%         prjFile:  File storing projection data - Format: float32
%         outFile:  Name of the file to store output data - Format: float32
%         nPrjs:    Number of projections in the set
%         prjSize:  Size of projections. Vector [N_rows, N_cols] (px)
%         pxSize:   Pixel size [S_row S_col] (mm)
%         volSize:  Size of volume. Vector [N_rows, N_cols, N_slices] (vx)
%         volOr:    Origin coordinate for the volume. Vector [or_x or_y or_z] (vx)
%         vxSize:   Voxel size [S_row S_col S_slice] (mm)
%         angSpan:  Angular span of the acquisition (deg)
%         startAng: First angular position (deg)
%         P:        Helix pitch (mm)
%         R0:       Helix radius (mm)
%         D:        Source-Detector-Distance (mm)
%
% Output:
%         The result of the calculation is stored in the file indicated
%         by "outFile"
%
%
% LIM - BiiG - UC3M
% Author: ASC
% Version 0 - Nov 2013
% Revision May 2014 - AOG
%

end
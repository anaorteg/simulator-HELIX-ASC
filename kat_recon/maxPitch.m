function P = maxPitch(Nrows,R0,D,dw,um)

%  
% maxPitch(Nrows,R0,D,dw,um)
%
% Function to estimate the maximum allowable pitch as a function of 
% system configuration.
% 
% Refs: 
%       Noo et al., Phys Med Biol, 48, 3787, 2003. 
%
% Input:
%         Nrows: Number of rows in the detector
%         R0:    Helix radius (mm)
%         D:     Source-Detector-Distance (mm)
%         dw:    Vertical pixel pitch (mm)
%         um:    Limit lateral coordinate of the detector (mm)
%
% Output:
%         P:     Maximum pitch (mm)
%
%
% LIM - BiiG - UC3M
% Author: ASC
% Version 0 - Nov 2013
%


  P = ((Nrows - 1)*(pi*R0*D*dw))/((um^2 + D^2)*(pi/2 + atan(um/D)));

end
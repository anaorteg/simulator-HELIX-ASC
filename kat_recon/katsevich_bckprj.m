function katsevich_bckprj(prjFile,outFile,nPrjs,prjSize,pxSize,volSize,volOr,vxSize,angSpan,startAng,P,R0,D)

%  
% katsevich_bckprj(prjFile,)
%
% Function to perform the last stages of the Katsevich's reconstruction
% algorithm, according to the approaches by Noo, Wunderlich and also Tan.
% 
% Refs: 
%       Noo et al., Phys Med Biol, 48, 3787, 2003. 
%       Wunderlich, Dept of Math Oregon Sta Univ., MS Thesis, 2006.
%       Tan et al., Med Phys 39, 4695, 2012.
%
% Input: - ToDo!!!
%         prjFile:  File storing processed projection data - Format: float32
%         nPrjs:    Number of projections in the set
%         prjSize:  Size of projections. Vector [N_rows, N_cols] (px)
%         pxSize:   Pixel size [S_row S_col] (mm)
%         angSpan:  Angular span of the acquisition (deg)
%         startAng: Starting angular point (deg)
%         P:        Helix pitch (mm)
%         R0:       Helux radius (mm)
%         D:        Source-Detector-Distance (mm)
%
% Output:
%
%
% LIM - BiiG - UC3M
% Author: ASC
% Version 0 - Nov 2013
%

% Projection sizes
uSize = prjSize(1);
wSize = prjSize(2);
uStep = pxSize(1);
wStep = pxSize(2);

% Volume sizes
xSize = volSize(1);
ySize = volSize(2);
zSize = volSize(3);
xStep = vxSize(1);
yStep = vxSize(2);
zStep = vxSize(3);
x_0   = volOr(1); % For the time being , we assumed this to be zero, only offset in z
y_0   = volOr(2); % For the time being , we assumed this to be zero, only offset in z
z_0   = volOr(3); % Offset in z

% Angular params - Be careful with deg-rad mistakes
angSpan_rd  = (angSpan*pi)/180.0;
lambdaStep  = angSpan_rd/nPrjs;
lambdaStart = (startAng*pi)/180.0;
lambdaEnd   = lambdaStart + angSpan_rd - lambdaStep;
% Generate vector for angular coords
angCoords   = lambdaStart:lambdaStep:lambdaEnd;

% Constant quantities
TD_cte = P/(2*pi*R0*D);
Dsq    = D^2;
pi2    = pi/2.0;
mag    = (D/R0);


% Storage
prjStack = zeros(uSize,wSize); % We need at least three projections to allow the backprojection according to Noo et al.
volStack = zeros(volSize);

% Open files
fdPrj = fopen(prjFile,'rb');

% Loop through projections, then, for each projection through the entire
% volume
% This is not an efficient approach in Matlab, but it is indeed the closest
% approach for implementing it in C/C++ and in CUDA, replacing the loops by
% CUDA kernels. Since the code is aimed at be implemented in CUDA, we
% prefer to get a non-optimized Matlab implementation easy to transcript to
% CUDA.
revStr = '';

% DEBUG -> show results
figure;
colormap gray;

for l = 1:nPrjs,
    
    msg = sprintf('Processing projection %d of %d',l,nPrjs);
    fprintf([revStr, msg]);
    revStr = repmat(sprintf('\b'), 1, length(msg));
    % Read data and update index
    tmp = fread(fdPrj,uSize*wSize,'float32');
    prjStack(:,:) = reshape(tmp,[uSize wSize]);
    
    % Calculations dependent only on angular position
    lambdaCur = angCoords(l);
    lambdaCos = cos(lambdaCur);
    lambdaSin = sin(lambdaCur);
    
    % Previous and next angles
    lambdaPre    = angCoords(l) - lambdaStep;
    lambdaCosPre = cos(lambdaPre);
    lambdaSinPre = sin(lambdaPre);
    lambdaPos    = angCoords(l) + lambdaStep;
    lambdaCosPos = cos(lambdaPos);
    lambdaSinPos = sin(lambdaPos);
    
    % Volume loop - first x, then y, last z -> See Noo etal for the reason
    for ix = 1:xSize,
        % Calculations not influenced by y neither z
        x = (ix - xSize/2 - 1) * xStep; % x coord
        
        % Y loop
        for iy = 1:ySize,
            % Calculations not influenced by z
            y = (iy - ySize/2 - 1) * yStep; % y coord
            % Detector indexes            
            v_st = R0 - x*lambdaCos - y*lambdaSin; % st stands for star, keep notation coherent with refs
            u_st = (D/v_st)*(-x*lambdaSin + y*lambdaCos);
            % Also for previous and post-projection
            v_st_pre = R0 - x*lambdaCosPre - y*lambdaSinPre;
            u_st_pre = (D/v_st_pre)*(-x*lambdaSinPre + y*lambdaCosPre);
            v_st_pos = R0 - x*lambdaCosPos - y*lambdaSinPos;
            u_st_pos = (D/v_st_pos)*(-x*lambdaSinPos + y*lambdaCosPos);
            
            % Tam-Danielsson window limits
            w_bot = -TD_cte * (u_st^2 + Dsq) * (pi2 + atan(u_st/D));
            w_top =  TD_cte * (u_st^2 + Dsq) * (pi2 - atan(u_st/D));
            z_bot = (v_st/D)*w_bot + (P/(2*pi))*lambdaCur;
            z_top = (v_st/D)*w_top + (P/(2*pi))*lambdaCur;
            
            % Compute limits for previous and post-projections and convert
            % them to z
            %w_bot_pre = -TD_cte * (u_st_pre^2 + Dsq) * (pi2 + atan(u_st_pre/D));
            w_bot_pos = -TD_cte * (u_st_pos^2 + Dsq) * (pi2 + atan(u_st_pos/D));
            w_top_pre =  TD_cte * (u_st_pre^2 + Dsq) * (pi2 - atan(u_st_pre/D));
            %w_top_pos =  TD_cte * (u_st_pos^2 + Dsq) * (pi2 - atan(u_st_pos/D));
            
            % Convert to z
            %z_bot_pre = (v_st/D)*w_bot_pre + (P/2*pi)*lambdaPre;
            z_bot_pos = (v_st_pos/D)*w_bot_pos + (P/(2*pi))*lambdaPos;
            z_top_pre = (v_st_pre/D)*w_top_pre + (P/(2*pi))*lambdaPre;
            %z_top_pos = (v_st/D)*w_top_pos + (P/2*pi)*lambdaPos;
            
            % Compute constant scale factor
            sca = lambdaStep/(2*pi*v_st);
                        
            % Z loop - For this we use the cone beam cover method described
            % by Fontaine
            % Not sure it is so accurate as the PI-line segment approach
            % but it is much more efficient
            z        = (wSize/(2*mag))*zStep - z_0;
            w_st     = (D/v_st)*(z - (P/(2*pi))*lambdaCur);
            w_st_inc = (D/v_st)*zStep;
            
            % Speedup the loop, loop only through the valid z
            
            
            for iz = 1:zSize,
                z = z - zStep; % z coord
                % Detector index
                w_st = w_st - w_st_inc;
                
                % Is the pixel within the Tam-Danielsson window?
                % or within the PI-line "extended" limits
                if (w_st > w_bot) && (w_st < w_top),
                    % Get weight
                    if (z < z_bot_pos),
                        weigh = abs(z_bot - z)/abs(z_bot - z_bot_pos) + 0.5;
                    else if (z > z_top_pre),
                            weigh = abs(z_top - z)/abs(z_top - z_top_pre) + 0.5;
                        else
                            weigh = 1;
                        end
                    end
                    
                    % Get the interpolated value
                    w_tmp = w_st/wStep + wSize/2;
                    u_tmp = u_st/uStep + uSize/2;
                    w_ind = floor(w_tmp);
                    u_ind = floor(u_tmp);
                    w_res = w_tmp - w_ind;
                    u_res = u_tmp - u_ind;
                    if (u_ind < 1) || (u_ind >= uSize) || (w_ind < 1) || (w_ind >= wSize),
                        continue;
                    end
                    
                    g_int = prjStack(u_ind,w_ind)*(1-u_res)*(1-w_res) + prjStack(u_ind+1,w_ind)*u_res*(1-w_res) + prjStack(u_ind,w_ind+1)*(1-u_res)*w_res + prjStack(u_ind+1,w_ind+1)*u_res*w_res;
                    
                    % Update voxel
                    volStack(ix,iy,iz) = volStack(ix,iy,iz) + sca*weigh*g_int;
                    
                end % End if
                
            end % End for z
        end % End for y
    end % End for x
    % Move to next prj
    imagesc(squeeze(volStack(:,:,zSize/2)));
    drawnow();
end % End for l

% Close prj file
fclose(fdPrj);

% Save volume file
fdVol = fopen(outFile,'wb');
fwrite(fdVol,volStack,'float32');
fclose(fdVol);

% Get a new line in the prompt
fprintf('\n');
end
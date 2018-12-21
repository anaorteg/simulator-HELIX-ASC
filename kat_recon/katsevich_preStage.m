function katsevich_preStage(prjFile,outFile,nPrjs,prjSize,pxSize,angSpan,P,R0,D)

%  
% katsevich_preStage(prjFile,nPrjs,angSpan,P,R0,prjSize)
%
% Function to perform the first stages of the Katsevich's reconstruction
% algorithm, according to the approaches by Noo, Wunderlich and Tan.
% 
% Refs: 
%       Noo et al., Phys Med Biol, 48, 3787, 2003. 
%       Wunderlich, Dept of Math Oregon Sta Univ., MS Thesis, 2006.
%       Tan et al., Med Phys 39, 4695, 2012.
%
% Input:
%         prjFile: File storing projection data - Format: float32
%         outFile: Name of the file to store output data - Format: float32
%         nPrjs:   Number of projections in the set
%         prjSize: Size of projections. Vector [N_rows, N_cols] (px)
%         pxSize:  Pixel size [S_row S_col] (mm)
%         angSpan: Angular span of the acquisition (deg)
%         P:       Helix pitch (mm)
%         R0:      Helix radius (mm)
%         D:       Source-Detector-Distance (mm)
%
% Output:
%         The result of the calculation is stored in the file indicated
%         by "outFile"
%
%
% LIM - BiiG - UC3M
% Author: ASC
% Version 0 - Nov 2013
% Author: AOG
% Version 0.1 - Junio 2014
%



% Notation:
%           Following Noo et al. we use the following notation:
%           
%           v      - Vector joining piercing point and source focus
%           u      - x detector coord, perpendicukar tu v
%           w      - z detector coord, perpendicular to u and v
%           lambda - Angular rotation for the projection
%           alpha  - Fan angle
%           D      - Source-Detector-Distance (SDD)
%           R0     - Helix radius (SAD) 
%           r      - FOV radius (Determined by size of detector here)
%           gx     - Data at step number x

% This first part of the method produces filtered images and corresponds to
% the tasks marked as stage 1 in Noo et al.
% The stage is subdivided into 5 substeps

% Get necessary data - For the time being dummy data
% Projection sizes
uSize = prjSize(1);
wSize = prjSize(2);
uStep = pxSize(1);
wStep = pxSize(2);
% Angular params - Be careful with deg-rad mistakes

% AO- angSpan must be an array:
%       one angular position per proyection
%       lambdaStep is also an array, used later for determining the mean
%           step as a function of the previous and post positions.
angSpan_rd = (angSpan.*pi)./180.0;
max_rd = 360*pi/180; %rad in 360º(0º)
 for l=1:nPrjs
        pos_l = l + 1;
        pre_l = l-1;
        if(pos_l>nPrjs||pre_l<=0)
        if (pos_l>nPrjs) 
            pos_l = (l+1) - nPrjs;
            lambdaStep(l)=abs(((max_rd+angSpan_rd(pos_l)) - angSpan_rd(l))+(angSpan_rd(l) - angSpan_rd(pre_l)))/2;
        end
        if (pre_l<=0)
            pre_l =(l-1) + nPrjs;
            lambdaStep(l)=abs((angSpan_rd(pos_l) - angSpan_rd(l))+(angSpan_rd(l) - (max_rd-angSpan_rd(pre_l))))/2;
        end
        else
        lambdaStep(l)=abs((angSpan_rd(pos_l) - angSpan_rd(l))+(angSpan_rd(l) - angSpan_rd(pre_l)))/2;
        end
 end %End of l

% Constant quantities
mag     = D./R0; % AO- R0 must be an array. Variable magnification due to R0 variability, D assumed constant
FOV_r   = ((uSize/2)*uStep)./mag;
alpha_m = mean(asin(FOV_r./R0),1); % AO - since the maximum difference in the vector is 10e-16, alpha_m can be considered as ctt
nKcurve = 2*wSize;
incPhi  = (pi + alpha_m)/nKcurve;
phiVal = zeros(nPrjs,nKcurve+1);
%for k = 1:nPrjs
phiVal  = (-nKcurve/2:1:nKcurve/2)*incPhi;
%end %end for k
 
% Compute remapping coords 
% for each projection - AO
WKcte   = (D*P)./((2*pi).*R0);
wk      = zeros(uSize,nKcurve+1,nPrjs);
indNaN       = find(phiVal == 0,1,'first');
% AO - there's a wk for each projection since theres a R0 for each angular
% position
for fr = 1:nPrjs 
    fr
for n = 1:uSize,
    u = (n - uSize/2)*uStep;
    for k = 1:(nKcurve+1)
        wk(n,k,fr) = WKcte(fr)*(phiVal(k) + (phiVal(k)/tan(phiVal(k)))*(u/D)); %AO - WKcte(k) is a vector 
        %disp([tan(phiVal(k)) wk(n,k,fr)]);
    end % End for k
end % End for n
wk(:,indNaN,fr) = (wk(:,indNaN-1,fr) + wk(:,indNaN+1,fr))/2.0;
%% 04/08 ALGO ESTA MAL, NO PUEDE HABER INDICES NEGATIVOS o > 512 Y WKIND LOS TIENE. a la derecha esta la version de alex 
% 18/08 Está relacionado con la geometría de nuestro sistema (puede que R>R0), si pongo los valores fijos de Alex el código funciona con las posiciones angulares  
wkInd(:,:,fr) = fix(wk(:,:,fr)/wStep) + wSize/2;
end % End for fr
% Filtering kernel - We gotta take into account the half pixel shift caused
% by the differentiation operation.
% ToDo the half-pixel shift
h_H = zeros(2*uSize,1);
for n = 1:2*uSize,
    ind_n = n - uSize;
    %ind_n = n;
    %if mod(n,2) == 0,
        %h_H(n) = 2/(pi*(n-1)*uStep);
        h_H(n) = 2/(pi*(ind_n-1)*uStep);
    %end
end
h_H(uSize+1) = 0;
% Perform the FFT
FFT_h_H = fft(h_H);
% % Half pixel shift
% FFT_h_H = FFT_h_H.*((exp((-1i*pi*(0:length(FFT_h_H)-1))*(1/length(FFT_h_H))))');
% Window the filter FFT
% Create the Hanning window
hanWin = zeros(size(FFT_h_H));
for n = 1:2*uSize,
    ind_n = n - uSize;
    hanWin(n) = 0.5 + 0.5*cos(pi*((abs(ind_n)-uSize)/(uSize)));
end
FFT_h_H = FFT_h_H.*hanWin;
% For posterior calculus
tmpRow = zeros(length(h_H),1);

% Find the limits of the Tam-Danielsson window
w_bot = zeros(uSize,nPrjs);
w_top = zeros(uSize,nPrjs);
for n = 1:uSize,
    u = (n - uSize/2)*uStep;
    w_bot(n,:) = -(P./(2*pi.*R0*D)) .* (u^2 + D^2) * ((pi/2)+atan(u/D)); %size (uSize,nFrames)
    w_top(n,:) = P./(2*pi.*R0*D) .* (u^2 + D^2) * ((pi/2)-atan(u/D)); %size (uSize,nFrames)
end

% Further storage
prjSet = zeros(uSize,wSize,2); % To store the two prjs to differentiate
% Storage intermediate steps
g1 = zeros(uSize,wSize);
g2 = zeros(uSize,wSize);
g3 = zeros(uSize,nKcurve+1);
g4 = zeros(uSize,nKcurve+1);
g5 = zeros(uSize,wSize);
% storage for interpolated weights
wkInd_aux = zeros(uSize,nKcurve+1);

% Open file for reading
fd    = fopen(prjFile,'rb');
 if fd ==-1
              disp(['Unable to open file: ' prjFile])
 end
fdOut = fopen(outFile,'wb');
 if fdOut ==-1
              disp(['Unable to open file: ' outFile])
 end


% The processing keeps looping till we have processed all the projection
% data.
% The loop approach is suboptimal in Matlab, but allows direct translation
% to C/C++ or CUDA and reduces memory burden. Since this code is just a
% proof of concept to serve as a basis for an optimized implementation, the
% loop approach provides a workflow similar to the one to implement in
% C/C++. 

% Read first projection
tmp           = fread(fd,uSize*wSize,'float32');
%%SOLO PARA QUE PUEDA COGER LAS IMAGENES SIMULADAS POR INES 31/10
tmp = tmp./200;
%end

prjSet(:,:,2) = reshape(tmp,[wSize uSize 1]);
% Perform I0 correction - Not implemented for the time being

% Perform log correction
prjSet(:,:,2) = -log(prjSet(:,:,2));

for nP = 2:nPrjs,
    disp(nP);
    drawnow();
    % Read next projection and move the previous one...
    % This would be much more efficient if we just modify the pointer, but
    % let's keep it like this for the sake of clarity -> take into account
    % for the final implementation
    tmp           = fread(fd,uSize*wSize,'float32');
    %%SOLO PARA QUE PUEDA COGER LAS IMAGENES SIMULADAS POR INES 31/10
tmp = log(tmp);
%end

    prjSet(:,:,1) = prjSet(:,:,2);
    prjSet(:,:,2) = reshape(tmp,[wSize uSize 1]);

    % Perform I0 correction - Not implemented for the time being
    
    % Perform log correction
    prjSet(:,:,2) = -log(prjSet(:,:,2));

    % The next section performs the following tasks
    %
    % FF1 - Angular differentiation
    %       Implemented using the chain rule and loops for the same reason
    %
    % FF2 - Length correction to compensate that the detector is flat
    
    for n = 1:uSize-1,
        % Compute coords
        u = (n - uSize/2)*uStep;
        % Inner loop - columns
        for j = 1:wSize-1,
            % Compute coords
            w = (j - wSize/2)*wStep;
            
            % FF1
            dLambda = ((prjSet(n,j,2)   - prjSet(n,j,1)) + (prjSet(n,j+1,2)   - prjSet(n,j+1,1)) + (prjSet(n+1,j,2) - prjSet(n+1,j,1)) + (prjSet(n+1,j+1,2) - prjSet(n+1,j+1,1))) / (4*lambdaStep(nP));%lambdaStep(nP) o lambdaStep(nP-1) o el promedio?
            du      = ((prjSet(n+1,j,1) - prjSet(n,j,1)) + (prjSet(n+1,j+1,1) - prjSet(n,j+1,1)) + (prjSet(n+1,j,2) - prjSet(n,j,2))   + (prjSet(n+1,j+1,2) - prjSet(n,j+1,2)))   / (4*uStep);
            dw      = ((prjSet(n,j+1,1) - prjSet(n,j,1)) + (prjSet(n+1,j+1,1) - prjSet(n+1,j,1)) + (prjSet(n,j+1,2) - prjSet(n,j,2))   + (prjSet(n+1,j+1,2) - prjSet(n+1,j,2)))   / (4*wStep);
            weighU  = ((u + uStep/2)^2 + D^2)/D;
            weighW  = ((u + uStep/2) * (w + wStep/2))/D;
            
            g1(n,j) = dLambda + weighU*du + weighW*dw;
            
            % FF2
            weighL  = D/sqrt(u^2 + D^2 + w^2);
            g2(n,j) = g1(n,j) * weighL;
            
        end % End for j
    end % End for n
      

    % This new section performs the filtering, on the rows remapped to
    % K-lines
    % FF3 - Remap w values to K lines.

    wkInd_aux = wkInd(:,:,nP)+wkInd(:,:,nP-1);
    wkInd_aux = fix(wkInd_aux/2.0);
    for n = 1:uSize,
        % we have several k-lines per row in prj...
        for k = 1:nKcurve+1,
            % Ok, wkInd is angular-position-dependent therefore should we
            % used the interpolation of both weights or just the one of the
            % first/second projection?
            % for the time being, I chose the interpolation option
            g3(n,k) = g2(n,wkInd_aux(n,k));
        end % End for k
    end % End for n

    % FF4 - Proper filtering but now along the k-lines
    %for k = 1:(2*nKcurve+1),
    for k = 1:nKcurve+1,
        tmpRow(1:length(g3(:,k)))                = 0;
        tmpRow(length(g3(:,k))+1:length(tmpRow)) = g3(:,k);
        % Perform FFT
        tmpFFT  = ifft(FFT_h_H.*fft(tmpRow));
        g4(:,k) = tmpFFT(1:length(g3(:,k)));
    end
    
    % FF5 - undo remapping
    %       We have to get back the results to the original grid, but in
    %       this case we first have to find the k-curve with the lowest phi
    for n = 2:uSize-1,
        % Compute coords
        u = (n - uSize/2)*uStep;
        % Inner loop - columns
        for j = 2:wSize-1,
            % Compute coords
            w = (j - wSize/2)*wStep;
            % Is w within the Tam-Danielsson window?
            if (w_bot(n) > w) || (w_top(n) < w),
                continue; % Outside of the window
            end
            
            % Find the proper phi value - the one with the lowest abs value
            % Two options depending on whether the u coord is negative or
            % positive
            if u < 0,
                %k = 2*nKcurve + 1;
                k = nKcurve + 1;
                while ((wk(n,k-1) < wk(n,k)) && (k > 2)),
                    if (wk(n,k-1) < w) && (wk(n,k) >= w),
                        % Compute interpolation values
                        c = (w - wk(n,k-1))/(wk(n,k) - wk(n,k-1));
                        % Perform interpolation
                        w_int = (1 - c)*g4(n,k-1) + c*g4(n,k);
                    end % End if
                    k = k - 1;
                end % End while
            else
                k = 1;
                while ((wk(n,k) < wk(n,k+1)) && (k < nKcurve)),
                    if (wk(n,k) < w) && (wk(n,k+1) >= w),
                        % Compute interpolation values
                        c = (w - wk(n,k))/(wk(n,k+1) - wk(n,k));
                        % Perform interpolation
                        w_int = (1 - c)*g4(n,k) + c*g4(n,k+1);
                    end % End if
                    k = k + 1;
                end % End while
            end % End if
            g5(n,j) = w_int;
        end % End for j
    end % End for n
    
    % Save result to output file
    fwrite(fdOut,g5,'float32');
end % End for nP

% Close files
fclose(fd);
fclose(fdOut);

end % The end


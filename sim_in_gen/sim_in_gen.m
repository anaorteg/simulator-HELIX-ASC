clc
close all
clear all
%% Calibration data
%Data
calibFilePath='C:/Users/aortega/Documents/GitHub/sim_in_gen/digimouse_calib.txt';
nProy=360;


proyAngle=zeros(nProy,1);
DSO=zeros(nProy,1);
DDO=zeros(nProy,1);
offsetX=zeros(nProy,1);
offsetY=zeros(nProy,1);
etaCalc=zeros(nProy,1);
thetaCalc=zeros(nProy,1);



%% Load Calibration File

    
    fileID = fopen(calibFilePath,'r');
    if fileID ~=-1
%        fgetl(fileID);
%         fgetl(fileID);  %Para cuando hay hueco
        
        for i=1:nProy
            buff = fgets(fileID);
            buff = sscanf(buff,'%f');
              % 30/07 cuidado, las , de Sedecal no se leen como floats, hay que cambiarlas a .   
            proyAngle(i)=buff(1);
            DSO(i)=buff(2);
            DDO(i)=buff(3);
            % DSO(i) = 1;
            offsetX(i)=buff(4);
            offsetY(i)=buff(5);
            etaCalc(i)=buff(6);
            thetaCalc(i)=buff(7);
            
%              fgetl(fileID);    %Para cuando hay hueco
        end
        
        fclose(fileID);
        res=0;
    else
        disp(['Unable to open file: ' calibFilePath])
    end


% System settings manually modify
% Set Spectrum_File = ..\spectra\spec_mouse_68kV02mmCu.txt
% Set ..\spectra\HamamatsuC7940.txt
% Voxel_Img_Filename = ..\vol\digimouse_flt_380x208x992.img  

% Base name convention: System_study_#projection.in
baseName = './mouse_real/in_mousereal/CTTB_study_'; 
z_origin = -14.5; % Insert float number for the origin in z axis
%px_size = 0.145;
px_size = 0.3;
zstep = 130/(nProy)*px_size; % bed displacement in px
z = z_origin;
fin = fopen('template.in');

%% Loop for every projection 
for it= 1:nProy
 % Return to the beginnning of the template
frewind(fin);    
% Create the new .in file
outName = strcat(baseName, num2str(it),'.in');
fout = fopen(outName,'a+');
% Writting the new file
while ~feof(fin)
   % Get line 
   s = fgetl(fin);
   % Modify initial angle (degrees)
   s = strrep(s, 'AngStart_Prj = ', strcat('AngStart_Prj = ', num2str(proyAngle(it))));
   % Modify geometric configuration
   s = strrep(s, 'SAD = ', strcat('SAD =  ',num2str(DSO(it)))); % mm
   s = strrep(s, 'u0 = ',  strcat('u0 = ', num2str(256.0+offsetX(it),'%.2f'))); % px
   s = strrep(s, 'v0 = ', strcat('v0 = ', num2str(256.0+offsetY(it),'%.2f'))); % px
   % Modify z: Volume_origin_z = float
   s = strrep(s, 'Volume_origin_z = ', strcat('Volume_origin_z = ', num2str(z,'%.2f'))); % px
   % Write down the resulting string in the new .in
   fprintf(fout,'%s',s);
   fprintf(fout,'\n');
   %disp(s)
end % modifying the .in file
z = z - zstep; % + or - depending on the sense of motion
fclose(fout);
end % projection iteration
fclose(fin);
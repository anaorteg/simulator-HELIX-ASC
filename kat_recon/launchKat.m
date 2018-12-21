clc
close all
clear all
%% Calibration data
%Data
calibFilePath='D:/MATLAB/kat_recon/kat_recon agosto/mangoose_sim/digimouse_calib.txt';
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
%            DSO(i) = 1;
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
     DSD = mean(DSO+DDO,1);
%% Call recon functions
%  maxPitch(1944,113,201,0.145,150) --> ans = 144.49
% katsevich_preStage('./p512/helix_digimouse_p512.img','./p512/helix_digimouse_p512_filt.img',360,[512 512],[0.3 0.3],proyAngle,130,DSO,DSD);
proyAngle_vec = proyAngle';
DSO_vec = DSO';
katsevich_preStage('./mangoose_sim/digiHel.ctf','./p512/digiHel_filt.img',360,[512 512],[0.3 0.3],proyAngle,130,DSO,DSD);
%katsevich_bckMex('./p512/helix_digimouse_p512_filt.img','./p512/helixvol_p512_130.img',359,[512 512],[0.3 0.3],[512 512 1024],[0 0 0],[0.18 0.18 0.18],proyAngle_vec,proyAngle(1),130,DSO_vec,DSD);
%katsevich_bckMex('./p512/digiHel_filt.img','./p512/helixvol_digiHel_130.img',359,[512 512],[0.3 0.3],[512 512 1024],[0 0 0],[0.18 0.18 0.18],proyAngle_vec,proyAngle(1),130,DSO_vec,DSD);





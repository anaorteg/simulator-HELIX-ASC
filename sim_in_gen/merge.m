%% Initialization
% Image size
sizeRow = 512;
sizeCol = 512;
% Number of proyections to merge
nProy = 360;
% Create the new output file
outName = strcat('./mouse_real/out_mousereal/digiHel.ctf');
fout = fopen(outName,'a+');
%% Loop for every projection 
for it= 1:nProy
% Opening each file
inName = strcat('./mouse_real/out_mousereal/0', num2str(it),'_primary.img')
fin = fopen(inName,'rb');
d=fread(fin,sizeRow*sizeCol,'float32');
% Copying data consecutively into the output file 
d=d./max(d);
fwrite(fout, d, 'float32');
%d=reshape(d,[sizeRow sizeCol sizeSlices]);
fclose(fin);
end % projection iteration
fclose(fout);
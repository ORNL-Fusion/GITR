file_re = 'acd50_w.dat';


fileID = fopen(file_re,'r');
tline = fgetl(fileID)
C = strsplit(tline, '/');
[IZMAX, IDMAXD, ITMAXD, IZ1MIN, IZ1MAX] = strread(C{1}, '%f %f %f %f %f');
tline = fgetl(fileID)
DDENSD = fscanf(fileID, '%f',IDMAXD)
DTEVD = fscanf(fileID, '%f',ITMAXD)
tline = fgetl(fileID);

RecombRateCoeff = zeros(ITMAXD,IDMAXD,ISEL);
for i=1:ISEL

tline = fgetl(fileID);


C = strsplit(tline, '/');
    
    D = strsplit(C{2}, '=');
    IPRT = str2num(D{2});
    
    D = strsplit(C{3}, '=');
    IGRD = str2num(D{2});
    
    D = strsplit(C{5}, '=');
    Z1 = str2num(D{2});
    
       if i==1
       RecombState = [IPRT IGRD Z1];
   
   else
       RecombState = [RecombState; [IPRT IGRD Z1]];
       
   end
    
    E = fscanf(fileID, '%f',IDMAXD*ITMAXD);
    E=transpose(reshape(E,[IDMAXD,ITMAXD]));
    
       
   RecombRateCoeff(:,:,i) = E;

    tline = fgetl(fileID);
    
end

fclose(fileID);
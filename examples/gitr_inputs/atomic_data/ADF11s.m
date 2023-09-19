function [Te, dens, RateCoeff, State] = ADF11s(file)
%file = 'ADAS/scd93_c.dat';


fileID = fopen(file,'r');
if fileID == -1
    print('Could not find ionization file')
end
tline = fgetl(fileID);
C = strsplit(tline, ' ');
IZMAX = str2num(C{2});
IDMAXD = str2num(C{3});
ITMAXD= str2num(C{4});
IZ1MIN= str2num(C{5});
IZ1MAX= str2num(C{6});

%NSEL = strread(tline, '%f %*s %*s %*s %*s %*s');




RateCoeff = zeros(IDMAXD,ITMAXD,IZ1MAX);
State = zeros(IZ1MAX,2);

tline = fgetl(fileID);

dens = fscanf(fileID,'%f',IDMAXD);
Te = fscanf(fileID, '%f',ITMAXD);
tline = fgetl(fileID);
for i=1:IZ1MAX
    tline = fgetl(fileID);
    C = strsplit(tline, '/');
    D = strsplit(C{2}, '=');
   
    IZ1 = str2num(D{2});
    IZ = IZ1 - 1;
   
 
RateCoeff(:,:,i) = fscanf(fileID,'%f',[IDMAXD,ITMAXD]);


    
    State(i,:) =[IZ IZ1];

 tline = fgetl(fileID);
end

fclose(fileID);



%plot(Te,RateCoeff(14,:,1))


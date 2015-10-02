clear variables
%User defined variables
file = 'ca09_w.dat';


fileID = fopen(file,'r');
tline = fgetl(fileID);
NSEL = strread(tline, '%f %*s %*s %*s %*s %*s');



for i=1:NSEL
    tline = fgetl(fileID);
    C = strsplit(tline, '/');
    D = strsplit(C{1}, '+');

    SYMB = D{1};
    IZ = str2num(D{2});
    
    D = strsplit(C{2}, '+');
    IZ1 = str2num(D{2});
    NTE = str2num(C{3});
    BWNO = strread(C{4}, '%*s %*s %f');
   METI = strread(C{5}, '%*s %*s %f');
   METF = strread(C{6}, '%*s %*s %f');
   ISEL = strread(C{7}, '%*s %*s %f');

   tline = fgetl(fileID);
   a = strread(tline, '%f');
   tline = fgetl(fileID);
   b = strread(tline, '%f');
   
   if i==1
       State = [IZ IZ1];
   Te = [a; b];
   else
       State = [State; [IZ IZ1]];
       Te = [Te [a;b]];
   end
   
      tline = fgetl(fileID);
   a = strread(tline, '%f');
   tline = fgetl(fileID);
   b = strread(tline, '%f');
   if i==1
   RateCoeff = [a;b];
   else
       RateCoeff = [RateCoeff [a;b]];
   end
end

fclose(fileID);







file = 'ADAS/ionelec_szd#c.dat';


fileID = fopen(file,'r');
if fileID == -1
    print('Could not find ionization file')
end
tline = fgetl(fileID);
NSEL = strread(tline, '%f %*s %*s %*s %*s %*s');

Te = [];
RateCoeff = [];
State = [];

for i=1:NSEL
    tline = fgetl(fileID);
    C = strsplit(tline, '/');
    D = strsplit(C{1}, '+');
    
    SYMB = D{1};
    IZ = str2num(D{2});
    
    D = strsplit(C{2}, '+');
    IZ1 = str2num(D{2});
    D = strsplit(C{5}, '=');
    ICODE = str2num(D{2});
    D = strsplit(C{6}, '=');
    FCODE = str2num(D{2});

    NTE = str2num(C{3});
    BWNO = strread(C{4}, '%*s %*s %f');
    METI = strread(C{5}, '%*s %*s %f');
    METF = strread(C{6}, '%*s %*s %f');
    ISEL = strread(C{7}, '%*s %*s %f');
    
    State = [State; [IZ ICODE IZ1 FCODE]];
    Te_in = [];
    RC_in = [];
    j=0;
    while j<NTE
        
        tline = fgetl(fileID);
        a = strread(tline, '%f');
       


            Te_in = [Te_in;a];

        j = j+length(a);
    end
    k=0;
    while k<NTE
        tline = fgetl(fileID);
        b = strread(tline, '%f');
        

            RC_in = [RC_in; b];
 
        k=  k+ length(b);
    end
  Te = [Te Te_in];
  RateCoeff = [RateCoeff RC_in];
end

fclose(fileID);






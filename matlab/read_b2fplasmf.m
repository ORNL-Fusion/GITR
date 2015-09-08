nx = 98;
ny = 38;
ns=9;

[filename,pathname] = uigetfile('b2fplasmf*');
fid = fopen(strcat(pathname,filename));

% Read in the first line, plus the labels on the second
% Then read in the first variable: bb
junk = fscanf(fid,'%s',3);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'bb')
    disp(['Error!: Expected bb, found ' label]);
    break
end
bb = fscanf(fid,'%e',nelem);

% Next variable: crx
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'crx')
    disp(['Error!: Expected crx, found ' label]);
    break
end
crx = fscanf(fid,'%e',nelem);

% Next variable: cry
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'cry')
    disp(['Error!: Expected cry, found ' label]);
    break
end
cry = fscanf(fid,'%e',nelem);

% Next variable: ffbz
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'ffbz')
    disp(['Error!: Expected ffbz, found ' label]);
    break
end
ffbz = fscanf(fid,'%e',nelem);

% Next variable: fpsi
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fpsi')
    disp(['Error!: Expected fpsi, found ' label]);
    break
end
fpsi = fscanf(fid,'%e',nelem);

% Next variable: gs
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'gs')
    disp(['Error!: Expected gs, found ' label]);
    break
end
gs = fscanf(fid,'%e',nelem);

% Next variable: hx
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'hx')
    disp(['Error!: Expected hx, found ' label]);
    break
end
hx = fscanf(fid,'%e',nelem);

% Next variable: hy
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'hy')
    disp(['Error!: Expected hy, found ' label]);
    break
end
hy = fscanf(fid,'%e',nelem);

% Next variable: qz
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'qz')
    disp(['Error!: Expected qz, found ' label]);
    break
end
qz = fscanf(fid,'%e',nelem);

% Next variable: qc
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'qc')
    disp(['Error!: Expected qc, found ' label]);
    break
end
qc = fscanf(fid,'%e',nelem);

% Next variable: vol
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'vol')
    disp(['Error!: Expected vol, found ' label]);
    break
end
vol = fscanf(fid,'%e',nelem);

% Next variable: pbs
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'pbs')
    disp(['Error!: Expected pbs, found ' label]);
    break
end
pbs = fscanf(fid,'%e',nelem);


% Next variable: fch
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fch')
    disp(['Error!: Expected fch, found ' label]);
    break
end
fch = fscanf(fid,'%e',nelem);

% Next variable: fch0
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fch0')
    disp(['Error!: Expected fch0, found ' label]);
    break
end
fch0 = fscanf(fid,'%e',nelem);

% Next variable: fchp
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fchp')
    disp(['Error!: Expected fchp, found ' label]);
    break
end
fchp = fscanf(fid,'%e',nelem);

% Next variable: fhe
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fhe')
    disp(['Error!: Expected fhe, found ' label]);
    break
end
fhe = fscanf(fid,'%e',nelem);

% Next variable: fhe0
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fhe0')
    disp(['Error!: Expected fhe0, found ' label]);
    break
end
fhe0 = fscanf(fid,'%e',nelem);

% Next variable: fhep
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fhep')
    disp(['Error!: Expected fhep, found ' label]);
    break
end
fhep = fscanf(fid,'%e',nelem);

% Next variable: fhet
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fhet')
    disp(['Error!: Expected fhet, found ' label]);
    break
end
fhet = fscanf(fid,'%e',nelem);

% Next variable: fhi
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fhi')
    disp(['Error!: Expected fhi, found ' label]);
    break
end
fhi = fscanf(fid,'%e',nelem);

% Next variable: fhi0
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fhi0')
    disp(['Error!: Expected fhi0, found ' label]);
    break
end
fhi0 = fscanf(fid,'%e',nelem);

% Next variable: fhip
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fhip')
    disp(['Error!: Expected fhip, found ' label]);
    break
end
fhip = fscanf(fid,'%e',nelem);

% Next variable: fhit
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fhit')
    disp(['Error!: Expected fhit, found ' label]);
    break
end
fhit = fscanf(fid,'%e',nelem);

% Next variable: fna
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fna')
    disp(['Error!: Expected fna, found ' label]);
    break
end
fna = fscanf(fid,'%e',nelem);

% Next variable: fna0
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fna0')
    disp(['Error!: Expected fna0, found ' label]);
    break
end
fna0 = fscanf(fid,'%e',nelem);

% Next variable: fnap
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fnap')
    disp(['Error!: Expected fnap, found ' label]);
    break
end
fnap = fscanf(fid,'%e',nelem);

% Next variable: fne
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fne')
    disp(['Error!: Expected fne, found ' label]);
    break
end
fne = fscanf(fid,'%e',nelem);

% Next variable: fni
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fni')
    disp(['Error!: Expected fni, found ' label]);
    break
end
fni = fscanf(fid,'%e',nelem);

% Next variable: na
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'na')
    disp(['Error!: Expected na, found ' label]);
    break
end
na = fscanf(fid,'%e',nelem);

% Next variable: na0
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'na0')
    disp(['Error!: Expected na0, found ' label]);
    break
end
na0 = fscanf(fid,'%e',nelem);

% Next variable: nap
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'nap')
    disp(['Error!: Expected nap, found ' label]);
    break
end
nap = fscanf(fid,'%e',nelem);

% Next variable: ne
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'ne')
    disp(['Error!: Expected ne, found ' label]);
    break
end
ne = fscanf(fid,'%e',nelem);

% Next variable: ne0
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'ne0')
    disp(['Error!: Expected ne0, found ' label]);
    break
end
ne0 = fscanf(fid,'%e',nelem);

% Next variable: ne2
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'ne2')
    disp(['Error!: Expected ne2, found ' label]);
    break
end
ne2 = fscanf(fid,'%e',nelem);

% Next variable: nep
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'nep')
    disp(['Error!: Expected nep, found ' label]);
    break
end
nep = fscanf(fid,'%e',nelem);

% Next variable: ni
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'ni')
    disp(['Error!: Expected ni, found ' label]);
    break
end
ni = fscanf(fid,'%e',nelem);

% Next variable: ni0
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'ni0')
    disp(['Error!: Expected ni0, found ' label]);
    break
end
ni0 = fscanf(fid,'%e',nelem);

% Next variable: pb
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'pb')
    disp(['Error!: Expected pb, found ' label]);
    break
end
pb = fscanf(fid,'%e',nelem);

% Next variable: po
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'po')
    disp(['Error!: Expected po, found ' label]);
    break
end
po = fscanf(fid,'%e',nelem);

% Next variable: po0
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'po0')
    disp(['Error!: Expected po0, found ' label]);
    break
end
po0 = fscanf(fid,'%e',nelem);

% Next variable: pop
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'pop')
    disp(['Error!: Expected pop, found ' label]);
    break
end
pop = fscanf(fid,'%e',nelem);

% Next variable: te
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'te')
    disp(['Error!: Expected te, found ' label]);
    break
end
te = fscanf(fid,'%e',nelem);

% Next variable: te0
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'te0')
    disp(['Error!: Expected te0, found ' label]);
    break
end
te0 = fscanf(fid,'%e',nelem);

% Next variable: tep
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'tep')
    disp(['Error!: Expected tep, found ' label]);
    break
end
tep = fscanf(fid,'%e',nelem);

% Next variable: ti
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'ti')
    disp(['Error!: Expected ti, found ' label]);
    break
end
ti = fscanf(fid,'%e',nelem);

% Next variable: ti0
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'ti0')
    disp(['Error!: Expected ti0, found ' label]);
    break
end
ti0 = fscanf(fid,'%e',nelem);

% Next variable: tip
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'tip')
    disp(['Error!: Expected tip, found ' label]);
    break
end
tip = fscanf(fid,'%e',nelem);

% Next variable: ua
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'ua')
    disp(['Error!: Expected ua, found ' label]);
    break
end
ua = fscanf(fid,'%e',nelem);

% Next variable: ua0
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'ua0')
    disp(['Error!: Expected ua0, found ' label]);
    break
end
ua0 = fscanf(fid,'%e',nelem);

% Next variable: uap
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'uap')
    disp(['Error!: Expected uap, found ' label]);
    break
end
uap = fscanf(fid,'%e',nelem);

% Next variable: uadia
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'uadia')
    disp(['Error!: Expected uadia, found ' label]);
    break
end
uadia = fscanf(fid,'%e',nelem);

% Next variable: fchdia
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fchdia')
    disp(['Error!: Expected fchdia, found ' label]);
    break
end
fchdia = fscanf(fid,'%e',nelem);

% Next variable: fmo
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fmo')
    disp(['Error!: Expected fmo, found ' label]);
    break
end
fmo = fscanf(fid,'%e',nelem);

% Next variable: fna_32
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fna_32')
    disp(['Error!: Expected fna_32, found ' label]);
    break
end
fna_32 = fscanf(fid,'%e',nelem);

% Next variable: fna_52
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fna_52')
    disp(['Error!: Expected fna_52, found ' label]);
    break
end
fna_52 = fscanf(fid,'%e',nelem);

% Next variable: fni_32
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fni_32')
    disp(['Error!: Expected fni_32, found ' label]);
    break
end
fni_32 = fscanf(fid,'%e',nelem);

% Next variable: fni_52
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fni_52')
    disp(['Error!: Expected fni_52, found ' label]);
    break
end
fni_52 = fscanf(fid,'%e',nelem);

% Next variable: fne_32
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fne_32')
    disp(['Error!: Expected fne_32, found ' label]);
    break
end
fne_32 = fscanf(fid,'%e',nelem);

% Next variable: fne_52
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fne_52')
    disp(['Error!: Expected fne_52, found ' label]);
    break
end
fne_52 = fscanf(fid,'%e',nelem);

% Next variable: wadia
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'wadia')
    disp(['Error!: Expected wadia, found ' label]);
    break
end
wadia = fscanf(fid,'%e',nelem);

% Next variable: vaecrb
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'vaecrb')
    disp(['Error!: Expected vaecrb, found ' label]);
    break
end
vaecrb = fscanf(fid,'%e',nelem);

% Next variable: facdrift
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'facdrift')
    disp(['Error!: Expected facdrift, found ' label]);
    break
end
facdrift = fscanf(fid,'%e',nelem);

% Next variable: fac_ExB
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fac_ExB')
    disp(['Error!: Expected fac_ExB, found ' label]);
    break
end
fac_ExB = fscanf(fid,'%e',nelem);

% Next variable: fchvispar
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fchvispar')
    disp(['Error!: Expected fchvispar, found ' label]);
    break
end
fchvispar = fscanf(fid,'%e',nelem);

% Next variable: fchvisper
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fchvisper')
    disp(['Error!: Expected fchvisper, found ' label]);
    break
end
fchvisper = fscanf(fid,'%e',nelem);

% Next variable: fchin
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fchin')
    disp(['Error!: Expected fchin, found ' label]);
    break
end
fchin = fscanf(fid,'%e',nelem);

% Next variable: fna_nodrift
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fna_nodrift')
    disp(['Error!: Expected fna_nodrift, found ' label]);
    break
end
fna_nodrift = fscanf(fid,'%e',nelem);

% Next variable: fac_vis
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fac_vis')
    disp(['Error!: Expected fac_vis, found ' label]);
    break
end
fac_vis = fscanf(fid,'%e',nelem);

% Next variable: resco
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'resco')
    disp(['Error!: Expected resco, found ' label]);
    break
end
resco = fscanf(fid,'%e',nelem);

% Next variable: reshe
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'reshe')
    disp(['Error!: Expected reshe, found ' label]);
    break
end
reshe = fscanf(fid,'%e',nelem);

% Next variable: reshi
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'reshi')
    disp(['Error!: Expected reshi, found ' label]);
    break
end
reshi = fscanf(fid,'%e',nelem);

% Next variable: resmo
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'resmo')
    disp(['Error!: Expected resmo, found ' label]);
    break
end
resmo = fscanf(fid,'%e',nelem);

% Next variable: resmt
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'resmt')
    disp(['Error!: Expected resmt, found ' label]);
    break
end
resmt = fscanf(fid,'%e',nelem);

% Next variable: respo
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'respo')
    disp(['Error!: Expected respo, found ' label]);
    break
end
respo = fscanf(fid,'%e',nelem);

% Next variable: sch
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'sch')
    disp(['Error!: Expected sch, found ' label]);
    break
end
sch = fscanf(fid,'%e',nelem);

% Next variable: she
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'she')
    disp(['Error!: Expected she, found ' label]);
    break
end
she = fscanf(fid,'%e',nelem);

% Next variable: shi
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'shi')
    disp(['Error!: Expected shi, found ' label]);
    break
end
shi = fscanf(fid,'%e',nelem);

% Next variable: smo
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'smo')
    disp(['Error!: Expected smo, found ' label]);
    break
end
smo = fscanf(fid,'%e',nelem);

% Next variable: smq
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'smq')
    disp(['Error!: Expected smq, found ' label]);
    break
end
smq = fscanf(fid,'%e',nelem);

% Next variable: sna
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'sna')
    disp(['Error!: Expected sna, found ' label]);
    break
end
sna = fscanf(fid,'%e',nelem);

% Next variable: sne
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'sne')
    disp(['Error!: Expected sne, found ' label]);
    break
end
sne = fscanf(fid,'%e',nelem);

% Next variable: rsana
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'rsana')
    disp(['Error!: Expected rsana, found ' label]);
    break
end
rsana = fscanf(fid,'%e',nelem);

% Next variable: rsahi
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'rsahi')
    disp(['Error!: Expected rsahi, found ' label]);
    break
end
rsahi = fscanf(fid,'%e',nelem);

% Next variable: rsamo
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'rsamo')
    disp(['Error!: Expected rsamo, found ' label]);
    break
end
rsamo = fscanf(fid,'%e',nelem);

% Next variable: rrana
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'rrana')
    disp(['Error!: Expected rrana, found ' label]);
    break
end
rrana = fscanf(fid,'%e',nelem);

% Next variable: rrahi
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'rrahi')
    disp(['Error!: Expected rrahi, found ' label]);
    break
end
rrahi = fscanf(fid,'%e',nelem);

% Next variable: rramo
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'rramo')
    disp(['Error!: Expected rramo, found ' label]);
    break
end
rramo = fscanf(fid,'%e',nelem);

% Next variable: rqahe
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'rqahe')
    disp(['Error!: Expected rqahe, found ' label]);
    break
end
rqahe = fscanf(fid,'%e',nelem);

% Next variable: rqrad
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'rqrad')
    disp(['Error!: Expected rqrad, found ' label]);
    break
end
rqrad = fscanf(fid,'%e',nelem);

% Next variable: rqbrm
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'rqbrm')
    disp(['Error!: Expected rqbrm, found ' label]);
    break
end
rqbrm = fscanf(fid,'%e',nelem);

% Next variable: rcxna
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'rcxna')
    disp(['Error!: Expected rcxna, found ' label]);
    break
end
rcxna = fscanf(fid,'%e',nelem);

% Next variable: rcxhi
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'rcxhi')
    disp(['Error!: Expected rcxhi, found ' label]);
    break
end
rcxhi = fscanf(fid,'%e',nelem);

% Next variable: rcxmo
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'rcxmo')
    disp(['Error!: Expected rcxmo, found ' label]);
    break
end
rcxmo = fscanf(fid,'%e',nelem);

% Next variable: b2stbr_sna
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2stbr_sna')
    disp(['Error!: Expected b2stbr_sna, found ' label]);
    break
end
b2stbr_sna = fscanf(fid,'%e',nelem);

% Next variable: b2stbr_smo
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2stbr_smo')
    disp(['Error!: Expected b2stbr_smo, found ' label]);
    break
end
b2stbr_smo = fscanf(fid,'%e',nelem);

% Next variable: b2stbr_she
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2stbr_she')
    disp(['Error!: Expected b2stbr_she, found ' label]);
    break
end
b2stbr_she = fscanf(fid,'%e',nelem);

% Next variable: b2stbr_shi
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2stbr_shi')
    disp(['Error!: Expected b2stbr_shi, found ' label]);
    break
end
b2stbr_shi = fscanf(fid,'%e',nelem);

% Next variable: b2stbr_sch
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2stbr_sch')
    disp(['Error!: Expected b2stbr_sch, found ' label]);
    break
end
b2stbr_sch = fscanf(fid,'%e',nelem);

% Next variable: b2stbr_sne
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2stbr_sne')
    disp(['Error!: Expected b2stbr_sne, found ' label]);
    break
end
b2stbr_sne = fscanf(fid,'%e',nelem);

% Next variable: b2stbc_sna
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2stbc_sna')
    disp(['Error!: Expected b2stbc_sna, found ' label]);
    break
end
b2stbc_sna = fscanf(fid,'%e',nelem);

% Next variable: b2stbc_smo
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2stbc_smo')
    disp(['Error!: Expected b2stbc_smo, found ' label]);
    break
end
b2stbc_smo = fscanf(fid,'%e',nelem);

% Next variable: b2stbc_she
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2stbc_she')
    disp(['Error!: Expected b2stbc_she, found ' label]);
    break
end
b2stbc_she = fscanf(fid,'%e',nelem);

% Next variable: b2stbc_shi
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2stbc_shi')
    disp(['Error!: Expected b2stbc_shi, found ' label]);
    break
end
b2stbc_shi = fscanf(fid,'%e',nelem);

% Next variable: b2stbc_sch
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2stbc_sch')
    disp(['Error!: Expected b2stbc_sch, found ' label]);
    break
end
b2stbc_sch = fscanf(fid,'%e',nelem);

% Next variable: b2stbc_sne
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2stbc_sne')
    disp(['Error!: Expected b2stbc_sne, found ' label]);
    break
end
b2stbc_sne = fscanf(fid,'%e',nelem);

% Next variable: b2stbm_sna
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2stbm_sna')
    disp(['Error!: Expected b2stbm_sna, found ' label]);
    break
end
b2stbm_sna = fscanf(fid,'%e',nelem);

% Next variable: b2stbm_smo
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2stbm_smo')
    disp(['Error!: Expected b2stbm_smo, found ' label]);
    break
end
b2stbm_smo = fscanf(fid,'%e',nelem);

% Next variable: b2stbm_she
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2stbm_she')
    disp(['Error!: Expected b2stbm_she, found ' label]);
    break
end
b2stbm_she = fscanf(fid,'%e',nelem);

% Next variable: b2stbm_shi
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2stbm_shi')
    disp(['Error!: Expected b2stbm_shi, found ' label]);
    break
end
b2stbm_shi = fscanf(fid,'%e',nelem);

% Next variable: b2stbm_sch
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2stbm_sch')
    disp(['Error!: Expected b2stbm_sch, found ' label]);
    break
end
b2stbm_sch = fscanf(fid,'%e',nelem);

% Next variable: b2stbm_sne
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2stbm_sne')
    disp(['Error!: Expected b2stbm_sne, found ' label]);
    break
end
b2stbm_sne = fscanf(fid,'%e',nelem);

% Next variable: b2sihs_divue
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2sihs_divue')
    disp(['Error!: Expected b2sihs_divue, found ' label]);
    break
end
b2sihs_divue = fscanf(fid,'%e',nelem);

% Next variable: b2sihs_divua
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2sihs_divua')
    disp(['Error!: Expected b2sihs_divua, found ' label]);
    break
end
b2sihs_divua = fscanf(fid,'%e',nelem);

% Next variable: b2sihs_exbe
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2sihs_exbe')
    disp(['Error!: Expected b2sihs_exbe, found ' label]);
    break
end
b2sihs_exbe = fscanf(fid,'%e',nelem);

% Next variable: b2sihs_exba
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2sihs_exba')
    disp(['Error!: Expected b2sihs_exba, found ' label]);
    break
end
b2sihs_exba = fscanf(fid,'%e',nelem);

% Next variable: b2sihs_visa
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2sihs_visa')
    disp(['Error!: Expected b2sihs_visa, found ' label]);
    break
end
b2sihs_visa = fscanf(fid,'%e',nelem);

% Next variable: b2sihs_joule
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2sihs_joule')
    disp(['Error!: Expected b2sihs_joule, found ' label]);
    break
end
b2sihs_joule = fscanf(fid,'%e',nelem);

% Next variable: b2sihs_fraa
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2sihs_fraa')
    disp(['Error!: Expected b2sihs_fraa, found ' label]);
    break
end
b2sihs_fraa = fscanf(fid,'%e',nelem);

% Next variable: b2sihs_str
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2sihs_str')
    disp(['Error!: Expected b2sihs_str, found ' label]);
    break
end
b2sihs_str = fscanf(fid,'%e',nelem);

% Next variable: b2npmo_smaf
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2npmo_smaf')
    disp(['Error!: Expected b2npmo_smaf, found ' label]);
    break
end
b2npmo_smaf = fscanf(fid,'%e',nelem);


% % Next variable: b2npht_qie
% junk = fscanf(fid,'%s',2);
% nelem = fscanf(fid,'%i',1);
% label = fscanf(fid,'%s',1);
% if ~strcmp(label,'b2npht_qie')
%     disp(['Error!: Expected b2npht_qie, found ' label]);
%     break
% end
% b2npht_qie = fscanf(fid,'%e',nelem);


% Next variable: b2npmo_smag
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2npmo_smag')
    disp(['Error!: Expected b2npmo_smag, found ' label]);
    break
end
b2npmo_smag = fscanf(fid,'%e',nelem);

% Next variable: b2npmo_smav
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'b2npmo_smav')
    disp(['Error!: Expected b2npmo_smav, found ' label]);
    break
end
b2npmo_smav = fscanf(fid,'%e',nelem);

% Next variable: ext_sna
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'ext_sna')
    disp(['Error!: Expected ext_sna, found ' label]);
    break
end
ext_sna = fscanf(fid,'%e',nelem);

% Next variable: ext_smo
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'ext_smo')
    disp(['Error!: Expected ext_smo, found ' label]);
    break
end
ext_smo = fscanf(fid,'%e',nelem);

% Next variable: ext_she
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'ext_she')
    disp(['Error!: Expected ext_she, found ' label]);
    break
end
ext_she = fscanf(fid,'%e',nelem);

% Next variable: ext_shi
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'ext_shi')
    disp(['Error!: Expected ext_shi, found ' label]);
    break
end
ext_shi = fscanf(fid,'%e',nelem);

% Next variable: ext_sch
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'ext_sch')
    disp(['Error!: Expected ext_sch, found ' label]);
    break
end
ext_sch = fscanf(fid,'%e',nelem);

% Next variable: ext_sne
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'ext_sne')
    disp(['Error!: Expected ext_sne, found ' label]);
    break
end
ext_sne = fscanf(fid,'%e',nelem);

% Next variable: calf
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'calf')
    disp(['Error!: Expected calf, found ' label]);
    break
end
calf = fscanf(fid,'%e',nelem);

% Next variable: cdna
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'cdna')
    disp(['Error!: Expected cdna, found ' label]);
    break
end
cdna = fscanf(fid,'%e',nelem);

% Next variable: cdpa
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'cdpa')
    disp(['Error!: Expected cdpa, found ' label]);
    break
end
cdpa = fscanf(fid,'%e',nelem);

% Next variable: ceqp
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'ceqp')
    disp(['Error!: Expected ceqp, found ' label]);
    break
end
ceqp = fscanf(fid,'%e',nelem);

% Next variable: chce
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'chce')
    disp(['Error!: Expected chce, found ' label]);
    break
end
chce = fscanf(fid,'%e',nelem);

% Next variable: chci
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'chci')
    disp(['Error!: Expected chci, found ' label]);
    break
end
chci = fscanf(fid,'%e',nelem);

% Next variable: chve
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'chve')
    disp(['Error!: Expected chve, found ' label]);
    break
end
chve = fscanf(fid,'%e',nelem);

% Next variable: chvemx
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'chvemx')
    disp(['Error!: Expected chvemx, found ' label]);
    break
end
chvemx = fscanf(fid,'%e',nelem);

% Next variable: chvi
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'chvi')
    disp(['Error!: Expected chvi, found ' label]);
    break
end
chvi = fscanf(fid,'%e',nelem);

% Next variable: chvimx
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'chvimx')
    disp(['Error!: Expected chvimx, found ' label]);
    break
end
chvimx = fscanf(fid,'%e',nelem);

% Next variable: csig
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'csig')
    disp(['Error!: Expected csig, found ' label]);
    break
end
csig = fscanf(fid,'%e',nelem);

% Next variable: cvla
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'cvla')
    disp(['Error!: Expected cvla, found ' label]);
    break
end
cvla = fscanf(fid,'%e',nelem);

% Next variable: cvsa
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'cvsa')
    disp(['Error!: Expected cvsa, found ' label]);
    break
end
cvsa = fscanf(fid,'%e',nelem);

% Next variable: cthe
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'cthe')
    disp(['Error!: Expected cthe, found ' label]);
    break
end
cthe = fscanf(fid,'%e',nelem);

% Next variable: cthi
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'cthi')
    disp(['Error!: Expected cthi, found ' label]);
    break
end
cthi = fscanf(fid,'%e',nelem);

% Next variable: csigin
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'csigin')
    disp(['Error!: Expected csigin, found ' label]);
    break
end
csigin = fscanf(fid,'%e',nelem);

% Next variable: cvsa_cl
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'cvsa_cl')
    disp(['Error!: Expected cvsa_cl, found ' label]);
    break
end
cvsa_cl = fscanf(fid,'%e',nelem);

% Next variable: fllime
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fllime')
    disp(['Error!: Expected fllime, found ' label]);
    break
end
fllime = fscanf(fid,'%e',nelem);

% Next variable: fllimi
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fllimi')
    disp(['Error!: Expected fllimi, found ' label]);
    break
end
fllimi = fscanf(fid,'%e',nelem);

% Next variable: fllim0fna
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fllim0fna')
    disp(['Error!: Expected fllim0fna, found ' label]);
    break
end
fllim0fna = fscanf(fid,'%e',nelem);

% Next variable: fllim0fhi
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fllim0fhi')
    disp(['Error!: Expected fllim0fhi, found ' label]);
    break
end
fllim0fhi = fscanf(fid,'%e',nelem);

% Next variable: fllimvisc
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'fllimvisc')
    disp(['Error!: Expected fllimvisc, found ' label]);
    break
end
fllimvisc = fscanf(fid,'%e',nelem);

% Next variable: sig0
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'sig0')
    disp(['Error!: Expected sig0, found ' label]);
    break
end
sig0 = fscanf(fid,'%e',nelem);

% Next variable: hce0
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'hce0')
    disp(['Error!: Expected hce0, found ' label]);
    break
end
hce0 = fscanf(fid,'%e',nelem);

% Next variable: alf0
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'alf0')
    disp(['Error!: Expected alf0, found ' label]);
    break
end
alf0 = fscanf(fid,'%e',nelem);

% Next variable: hci0
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'hci0')
    disp(['Error!: Expected hci0, found ' label]);
    break
end
hci0 = fscanf(fid,'%e',nelem);

% Next variable: hcib
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'hcib')
    disp(['Error!: Expected hcib, found ' label]);
    break
end
hcib = fscanf(fid,'%e',nelem);

% Next variable: dpa0
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'dpa0')
    disp(['Error!: Expected dpa0, found ' label]);
    break
end
dpa0 = fscanf(fid,'%e',nelem);

% Next variable: dna0
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'dna0')
    disp(['Error!: Expected dna0, found ' label]);
    break
end
dna0 = fscanf(fid,'%e',nelem);

% Next variable: vsa0
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'vsa0')
    disp(['Error!: Expected vsa0, found ' label]);
    break
end
vsa0 = fscanf(fid,'%e',nelem);

% Next variable: vla0
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'vla0')
    disp(['Error!: Expected vla0, found ' label]);
    break
end
vla0 = fscanf(fid,'%e',nelem);

% Next variable: csig_an
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'csig_an')
    disp(['Error!: Expected csig_an, found ' label]);
    break
end
csig_an = fscanf(fid,'%e',nelem);

% Next variable: calf_an
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'calf_an')
    disp(['Error!: Expected calf_an, found ' label]);
    break
end
calf_an = fscanf(fid,'%e',nelem);

% Next variable: nstra
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'nstra')
    disp(['Error!: Expected nstra, found ' label]);
    break
end
nstra = fscanf(fid,'%e',nelem);

% Next variable: sclstra
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'sclstra')
    disp(['Error!: Expected sclstra, found ' label]);
    break
end
sclstra = fscanf(fid,'%e',nelem);

% Next variable: sclrtio
junk = fscanf(fid,'%s',2);
nelem = fscanf(fid,'%i',1);
label = fscanf(fid,'%s',1);
if ~strcmp(label,'sclrtio')
    disp(['Error!: Expected sclrtio, found ' label]);
    break
end
sclrtio = fscanf(fid,'%e',nelem);

% Reshape the variables into 2D arrays
bb = reshape(bb,[nx ny 4]);
crx = reshape(crx,[nx ny 4]);
cenx = mean(crx,3);
cry = reshape(cry,[nx ny 4]);
ceny = mean(cry,3);
ffbz = reshape(ffbz,[nx ny 4]);
fpsi = reshape(fpsi,[nx ny 4]);
gs = reshape(gs,[nx ny 3]);
hx = reshape(hx,[nx ny]);
hy = reshape(hy,[nx ny]);
qz = reshape(qz,[nx ny 2]);
qc = reshape(qc,[nx ny]);
vol = reshape(vol,[nx ny]);
fch = reshape(fch,[nx ny 2]);
fch0 = reshape(fch0,[nx ny 2]);
fchp = reshape(fchp,[nx ny 2]);
fhe = reshape(fhe,[nx ny 2]);
fhe0 = reshape(fhe0,[nx ny 2]);
fhep = reshape(fhep,[nx ny 2]);
fhi = reshape(fhi,[nx ny 2]);
fhi0 = reshape(fhi0,[nx ny 2]);
fhip = reshape(fhip,[nx ny 2]);
fhit = reshape(fhit,[nx ny 2]);
fna = reshape(fna,[nx ny 2 ns]);
fna0 = reshape(fna0,[nx ny 2 ns]);
fnap = reshape(fnap,[nx ny 2 ns]);
fne = reshape(fne,[nx ny 2]);
fni = reshape(fni,[nx ny 2]);
na = reshape(na,[nx ny ns]);
na0 = reshape(na0,[nx ny ns]);
nap = reshape(nap,[nx ny ns]);
ne = reshape(ne,[nx ny]);
ne0 = reshape(ne0,[nx ny]);
ne2 = reshape(ne2,[nx ny]);
nep = reshape(nep,[nx ny]);
ni = reshape(ni,[nx ny 2]);
ni0 = reshape(ni0,[nx ny 2]);
pb = reshape(pb,[nx ny]);
po = reshape(po,[nx ny]);
po0 = reshape(po0,[nx ny]);
pop = reshape(pop,[nx ny]);
te = reshape(te,[nx ny]);
te0 = reshape(te0,[nx ny]);
tep = reshape(tep,[nx ny]);
ti = reshape(ti,[nx ny]);
ti0 = reshape(ti0,[nx ny]);
tip = reshape(tip,[nx ny]);
ua = reshape(ua,[nx ny ns]);
ua0 = reshape(ua0,[nx ny ns]);
uap = reshape(uap,[nx ny ns]);
uadia = reshape(uadia,[nx ny 2 ns]);
fchdia = reshape(fchdia,[nx ny 2]);
fmo = reshape(fmo,[nx ny 2 ns]);
fna_32 = reshape(fna_32,[nx ny 2 ns]);
fna_52 = reshape(fna_52,[nx ny 2 ns]);
fni_32 = reshape(fni_32,[nx ny 2]);
fni_52 = reshape(fni_52,[nx ny 2]);
fne_32 = reshape(fne_32,[nx ny 2]);
fne_52 = reshape(fne_52,[nx ny 2]);
wadia = reshape(wadia,[nx ny 2 ns]);
vaecrb = reshape(vaecrb,[nx ny 2 ns]);
facdrift = reshape(facdrift,[nx ny]);
fac_ExB = reshape(fac_ExB,[nx ny]);
fchvispar = reshape(fchvispar,[nx ny 2]);
fchvisper = reshape(fchvisper,[nx ny 2]);
fchin = reshape(fchin,[nx ny 2]);
fna_nodrift = reshape(fna_nodrift,[nx ny 2 ns]);
fac_vis = reshape(fac_vis,[nx ny]);
resco = reshape(resco,[nx ny ns]);
reshe = reshape(reshe,[nx ny]);
reshi = reshape(reshi,[nx ny]);
resmo = reshape(resmo,[nx ny ns]);
resmt = reshape(resmt,[nx ny]);
respo = reshape(respo,[nx ny]);
sch = reshape(sch,[nx ny 4]);
she = reshape(she,[nx ny 4]);
shi = reshape(shi,[nx ny 4]);
smo = reshape(smo,[nx ny 4 ns]);
smq = reshape(smq,[nx ny 4 ns]);
sna = reshape(sna,[nx ny 2 ns]);
rsana = reshape(rsana,[nx ny ns]);
rsahi = reshape(rsahi,[nx ny ns]);
rsamo = reshape(rsamo,[nx ny ns]);
rrana = reshape(rrana,[nx ny ns]);
rrahi = reshape(rrahi,[nx ny ns]);
rramo = reshape(rramo,[nx ny ns]);
rqahe = reshape(rqahe,[nx ny ns]);
rqrad = reshape(rqrad,[nx ny ns]);
rqbrm = reshape(rqbrm,[nx ny ns]);
rcxna = reshape(rcxna,[nx ny ns]);
rcxhi = reshape(rcxhi,[nx ny ns]);
rcxmo = reshape(rcxmo,[nx ny ns]);
b2stbr_sna = reshape(b2stbr_sna,[nx ny ns]);
b2stbr_smo = reshape(b2stbr_smo,[nx ny ns]);
b2stbr_she = reshape(b2stbr_she,[nx ny]);
b2stbr_shi = reshape(b2stbr_shi,[nx ny]);
b2stbr_sch = reshape(b2stbr_sch,[nx ny]);
b2stbc_sna = reshape(b2stbc_sna,[nx ny ns]);
b2stbc_smo = reshape(b2stbc_smo,[nx ny ns]);
b2stbc_she = reshape(b2stbc_she,[nx ny]);
b2stbc_shi = reshape(b2stbc_shi,[nx ny]);
b2stbc_sch = reshape(b2stbc_sch,[nx ny]);
b2sihs_divue = reshape(b2sihs_divue,[nx ny]);
b2sihs_divua = reshape(b2sihs_divua,[nx ny]);
b2sihs_exbe = reshape(b2sihs_exbe,[nx ny]);
b2sihs_exba = reshape(b2sihs_exba,[nx ny]);
b2sihs_visa = reshape(b2sihs_visa,[nx ny]);
b2sihs_joule = reshape(b2sihs_joule,[nx ny]);
b2sihs_fraa = reshape(b2sihs_fraa,[nx ny]);
b2sihs_str = reshape(b2sihs_str,[nx ny]);
ext_sna = reshape(ext_sna,[nx ny ns]);
ext_smo = reshape(ext_smo,[nx ny ns]);
ext_she = reshape(ext_she,[nx ny]);
ext_shi = reshape(ext_shi,[nx ny]);
ext_sch = reshape(ext_sch,[nx ny]);
calf = reshape(calf,[nx ny 2]);
cdna = reshape(cdna,[nx ny 2 ns]);
cdpa = reshape(cdpa,[nx ny 2 ns]);
ceqp = reshape(ceqp,[nx ny]);
chce = reshape(chce,[nx ny 2]);
chci = reshape(chci,[nx ny 2]);
chve = reshape(chve,[nx ny 2]);
chvemc = reshape(chvemx,[nx ny]);
chvi = reshape(chvi,[nx ny 2]);
chvimc = reshape(chvimx,[nx ny]);
csig = reshape(csig,[nx ny 2]);
cvla = reshape(cvla,[nx ny 2 ns]);
cvsa = reshape(cvsa,[nx ny 2 ns]);
cthe = reshape(cthe,[nx ny ns]);
cthi = reshape(cthi,[nx ny ns]);
csigin = reshape(csigin,[nx ny 2 ns ns]);
cvsa_cl = reshape(cvsa_cl,[nx ny 2 ns]);
fllime = reshape(fllime,[nx ny]);
fllimi = reshape(fllimi,[nx ny]);
sig0 = reshape(sig0,[nx ny]);
hce0 = reshape(hce0,[nx ny]);
alf0 = reshape(alf0,[nx ny]);
hci0 = reshape(hci0,[nx ny]);
hcib = reshape(hcib,[nx ny ns]);
dpa0 = reshape(dpa0,[nx ny ns]);
dna0 = reshape(dna0,[nx ny ns]);
vsa0 = reshape(vsa0,[nx ny ns]);
vla0 = reshape(vla0,[nx ny 2 ns]);
csig_an = reshape(csig_an,[nx ny 2]);

fclose(fid);
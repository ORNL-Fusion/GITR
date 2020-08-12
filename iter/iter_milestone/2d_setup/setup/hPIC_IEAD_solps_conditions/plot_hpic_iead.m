% IEAD
clear all 
close all
clc
plot = 0;
% Physical constants
qe   =  1.60217662e-19;         % Fundamental charge
me   =  9.10938356e-31;         % Electron mass
mp   =  1.67262190e-27;         % Proton mass
KB   =  1.38064852e-23;         % Boltzmann constant
lux  =  299792458;              % Speed of light
eps0 =  1.0/lux/lux/4/pi*1e7;   % Dielectric const
d2r  =  pi/180;

% Target Folder 
TargetFolder = 'hPIC_IEAD_DATA';

% Read SOLPS data
SOLPS = csvread([TargetFolder '/ITER_SOLPS.csv']);

% List of SOLPS Points resolved with hPIC 
LocationList = [ 1:1:38 ];

minEEall = 5;
dE = minEEall/90;
maxEEall = 1250;
Egrid = 0:dE:maxEEall;
Agrid = 0.5:1:89.5;
[xq,yq] = meshgrid(Egrid,Agrid);
IEADs = zeros(length(Egrid),length(Agrid),38);
% minE = 1000;
% maxE = 0;
for i = 1:numel( LocationList )
    
    ii = LocationList(i);
    
    Ai  = SOLPS(ii,1);    % Ion Atomic Mass 
    Zi  = SOLPS(ii,2);    % Ion Charge Number 
    
    Te  = SOLPS(ii,3);      % Electron Temperature [eV] 
    Ti  = SOLPS(ii,4);      % Ion Temperature [eV] 
    B0  = SOLPS(ii,5);      % Magnetic Field [T] 
    psi = SOLPS(ii,6);  % Magnetic Angle [deg] 
    n0  = SOLPS(ii,7);      % Ion Density [m^-3] 
        
    ID = ['SolpsPoint' num2str(ii)];
if plot
    % Figure IEAD - Ion Energy-Angle Distribution
    close all

    figure(1)
    FontSizeAxes = 12;          %  Size of Axes font
    set(gcf,'defaultaxesfontsize',FontSizeAxes)
    set(gcf,'defaultaxesfontname','Arial')
    set(gcf,'defaulttextcolor','black')
end
    EE   = dlmread( [TargetFolder '/' ID '_EnergyCenters.dat'] );
    AA   = dlmread( [TargetFolder '/' ID '_AngleCenters.dat'] );
    IEAD = dlmread( [TargetFolder '/' ID '_IEAD.dat'] );
         [x,y] = meshgrid(EE,AA);
        IEADs(:,:,i) = griddata(x,y,IEAD',xq,yq,'linear')';
        IEADs(find(Egrid > EE(end)),:,i) = 0;
%     if(EE(end)<minE)
%         minE = EE(end);
%     end
%     
%     if(EE(end)>maxE)
%         maxE = EE(end);
%     end
    
if plot
    h  = pcolor(AA,EE,IEAD); 
    hold on
    set(h, 'EdgeColor', 'none');
    shading interp
    colorbar
    % Classical Phi_floating and Presheath, just as a Guide for the Eye
    Mi = Ai * mp;
    phi_floating = Te * log( (Mi/2/pi/me/(1+Ti/Te)).^0.5 );
    plot ( [0 90], phi_floating*[1 1], 'r', 'Linewidth',1.0)
    plot ( [0 90], (Te/2+phi_floating)*[1 1], 'r', 'Linewidth',1.0)
    xlim([0 90])
    xlabel('\theta [deg]')
    ylabel('E [eV]')
    title(['Ions Mi=' num2str(Mi/mp,'%2.2g') ', Te=',num2str(Te,'%2.4g') 'eV, Ti=',num2str(Ti,'%2.4g') 'eV, B=',num2str(B0,'%2.4g') 'T, \psi=',num2str(psi,'%2.4g') 'deg, n=',num2str(n0,'%2.2e') 'm^{-3}'])
    print('-f1', '-dpng', [ ID '_IEAD'] )
end
        
end

z=table2array(solpsTarg(:,3));
r=table2array(solpsTarg(:,2));
rSep = 5.5543;
zSep = -4.391;

[val ind] = min(abs(z-zSep))
rSepIndx = ind;
len = zeros(length(r)-1,1);
for i=1:length(len)
    len(i) = sqrt((r(i+1) - r(i))^2 + (z(i+1) - z(i))^2);
end

rMrs = 0*r;
offset = sum(len(1:rSepIndx))+(zSep - z(rSepIndx));
rMrs = [0 cumsum(len)'] - offset;

nL = 38;
nS = 2;
nE = length(Egrid);
nA = 90;
As = [4 4];
ncid = netcdf.create(['./ieadsHpicMq3.nc'],'NC_WRITE')
 
dimL = netcdf.defDim(ncid,'nLocations',nL);
% dimS = netcdf.defDim(ncid,'nSpecies',nS);
dimE = netcdf.defDim(ncid,'nE',nE);
dimA = netcdf.defDim(ncid,'nA',nA);

% AVar = netcdf.defVar(ncid,'A','float',[dimS]);
% ZVar = netcdf.defVar(ncid,'Z','float',[dimS]);
rmrsVar = netcdf.defVar(ncid,'rMrs','float',[dimL]);

gridEVar = netcdf.defVar(ncid,'gridE','float',[dimE]);
gridAVar = netcdf.defVar(ncid,'gridA','float',[dimA]);
% fluxVar = netcdf.defVar(ncid,'flux','float',[dimS dimL]);
% densVar = netcdf.defVar(ncid,'dens','float',[dimS dimL]);

distVar = netcdf.defVar(ncid,'IEAD','float',[dimE dimA dimL]);


netcdf.endDef(ncid);
% netcdf.putVar(ncid, AVar, As);
% netcdf.putVar(ncid, ZVar, Zis);
netcdf.putVar(ncid, gridEVar, Egrid);
netcdf.putVar(ncid, gridAVar, Agrid);
netcdf.putVar(ncid, rmrsVar, rMrs);
% netcdf.putVar(ncid, densVar, nn);
% netcdf.putVar(ncid, fluxVar, ff);
netcdf.putVar(ncid, distVar, IEADs);


netcdf.close(ncid);
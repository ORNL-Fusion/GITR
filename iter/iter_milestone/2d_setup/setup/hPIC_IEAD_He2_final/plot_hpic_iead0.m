% IEAD
clear all 
close all
clc

% Physical constants
qe   =  1.60217662e-19;         % Fundamental charge
me   =  9.10938356e-31;         % Electron mass
mp   =  1.67262190e-27;         % Proton mass
KB   =  1.38064852e-23;         % Boltzmann constant
lux  =  299792458;              % Speed of light
eps0 =  1.0/lux/lux/4/pi*1e7;   % Dielectric const
d2r  =  pi/180;

% Specify the following three inputs for the figures: 
% 
%   Target Folder
%   Ai, Ion Atomic Mass
%   Zi, Charge State
%
TargetFolder = 'hPIC_IEAD_He2'; % 'hPIC_IEAD_He1'
Ai  = 4;    % Ion Atomic Mass 
Zis  = [1 2];    % Ion Charge Number 
plotting=0;
% Read SOLPS data (skip first row)
% SOLPS = csvread('solpsTarg.csv',1,0);
% SOLPS = csvread([TargetFolder '/ITER_SOLPS.csv']);
readSolpsTarg
% List of SOLPS Points resolved with hPIC
LocationList = [ 1:38 ];
minEEall = 5.6890;
dE = minEEall/90;
maxEEall = 1000;
Egrid = 0:dE:maxEEall;
Agrid = 0.5:1:89.5;
[xq,yq] = meshgrid(Egrid,Agrid);
IEADs = zeros(length(Egrid),length(Agrid),38);

for i = 1:numel( LocationList )
    i
    for j=1:length(Zis)
        Zi = Zis(j);
        ii = LocationList(i);
        
        Te  = SOLPS(ii,4);      % Electron Temperature [eV]
        Ti  = SOLPS(ii,5);      % Ion Temperature [eV]
        B0  = SOLPS(ii,12);     % Magnetic Field [T]
        psi = SOLPS(ii,13);     % Magnetic Angle [deg]
        if Zi==1
            n0  = SOLPS(ii,10); % Ion Density He1+ [m^-3]
        elseif Zi==2
            n0  = SOLPS(ii,9); % Ion Density He2+ [m^-3]
        end
        
        ID = ['SolpsPoint' num2str(ii) ];
        if plotting
        % Figure IEAD - Ion Energy-Angle Distribution
        close all
        
        figure(1)
        FontSizeAxes = 12;          %  Size of Axes font
        set(gcf,'defaultaxesfontsize',FontSizeAxes)
        set(gcf,'defaultaxesfontname','Arial')
        set(gcf,'defaulttextcolor','black')
        end
        EE   = dlmread( [TargetFolder '/' ID '_EnergyCenters.dat'] );
%         maxEE = max(EE(end));
%         if (maxEE > maxEEall)
%             maxEEall = maxEE;
%         end
        AA   = dlmread( [TargetFolder '/' ID '_AngleCenters.dat'] );
        IEAD = dlmread( [TargetFolder '/' ID '_IEAD.dat'] );
        [x,y] = meshgrid(EE,AA);
        IEADs(:,:,i) = griddata(x,y,IEAD',xq,yq,'linear')';
        IEADs(find(Egrid > EE(end)),:,i) = 0;
        
        
        if plotting
            h  = pcolor(AA,EE,IEAD);
            hold on
            set(h, 'EdgeColor', 'none');
            shading interp
            colorbar
            % Classical Phi_floating and Presheath, just as a Guide for the Eye
            Mi = Ai * mp;
            phi_floating = Te * ( 0.5*log(Mi/2/pi/me/(1+Ti/Te)) - log(Zi));
            phi_total = phi_floating + 0.5*Te;
            
            plot([0, 90],[Zi*phi_floating, Zi*phi_floating], 'r')
            plot([0, 90],[Zi*phi_total, Zi*phi_total], 'r')
            xlim([0 90])
            xlabel('\theta [deg]')
            ylabel('E [eV]')
            title(['Ai=' num2str(Ai,'%2.2g') ', Zi=' num2str(Zi,'%2.2g') ' Te=',num2str(Te,'%2.4g') 'eV, Ti=',num2str(Ti,'%2.4g') 'eV, B=',num2str(B0,'%2.4g') 'T, \psi=',num2str(psi,'%2.4g') 'deg, n=',num2str(n0,'%2.2e') 'm^{-3}'])
            print('-f1', '-dpng', [ ID '_He' num2str(Zi) '_IEAD'] )
        end
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

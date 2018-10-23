clear all
close all
nEdistBins = 100;
Eb = 8;
energy = linspace(0,Eb,nEdistBins);
eDist = zeros(nEdistBins);
eDist(end)=1;
% alpha = 2;
% thompson = energy./(energy+Eb).^3;
% thompson = alpha*(alpha-1)*Eb^(alpha-1)*energy./(energy+Eb).^(alpha+1);
% 
% figure(6)
% plot(energy,thompson)
% 
% cumulativeProbE = cumsum(thompson);
% eCDF = cumulativeProbE./cumulativeProbE(end);
% eCDF=1;
% figure(7)
% plot(energy,eCDF)

nAdistBins = 200;
% eAngle = linspace(0,pi/2,nAdistBins);
% eAngleDist = cos(eAngle);
% eAngle=[0,0];
% figure(9)
% plot(eAngle,eAngleDist)
% 
% eAngleCDF = cumsum(eAngleDist);
% eAngleCDF = eAngleCDF./eAngleCDF(end);

% figure(10)
% plot(eAngle,eAngleCDF)

% eAngleCDF(2:end) = eAngleCDF(1:end-1);
% eAngleCDF(1) = 0;
nE = 2;
nA = 3;
E=[0.001, 10000];
A=[0,45,89.99];

spyld = 0.1*ones(nE,nA);
rfyld = 0.5*ones(nE,nA);

angDists=zeros(nAdistBins,nA,nE);
eDists=zeros(nEdistBins,nA,nE);
eDists(end,:,:)=1;
phiGrid=linspace(0,90,nAdistBins);
phiDists = angDists;
phiDists(100,:,:) = 1;
thetaGrid=linspace(0,180,nAdistBins);
thetaDists = angDists;
thetaDists(100,:,:) = 1;
ncid = netcdf.create(['./simpleSurfaceModel8ev.nc'],'NC_WRITE')
 
dimE = netcdf.defDim(ncid,'nE',nE);
dimA = netcdf.defDim(ncid,'nA',nA);
dimEdistBins = netcdf.defDim(ncid,'nEdistBins',nEdistBins);
dimEdistBinsRef = netcdf.defDim(ncid,'nEdistBinsRef',nEdistBins);
dimAdistBins = netcdf.defDim(ncid,'nAdistBins',nAdistBins);

EVar = netcdf.defVar(ncid,'E','float',[dimE]);
AVar = netcdf.defVar(ncid,'A','float',[dimA]);
spyldVar= netcdf.defVar(ncid,'spyld','float',[dimE dimA]);
rfyldVar= netcdf.defVar(ncid,'rfyld','float',[dimE dimA]);
cosXdistVar=netcdf.defVar(ncid,'cosXDist','float',[dimE dimA dimAdistBins]);
cosYdistVar=netcdf.defVar(ncid,'cosYDist','float',[dimE dimA dimAdistBins]);
cosXdistRefVar=netcdf.defVar(ncid,'cosXDistRef','float',[dimE dimA dimAdistBins]);
cosYdistRefVar=netcdf.defVar(ncid,'cosYDistRef','float',[dimE dimA dimAdistBins]);

edistVar=netcdf.defVar(ncid,'energyDist','float',[dimE dimA dimEdistBins]);
edistRefVar=netcdf.defVar(ncid,'energyDistRef','float',[dimE dimA dimEdistBinsRef]);

eDistGridVar = netcdf.defVar(ncid,'eDistEgrid','float',[dimEdistBins]);
eDistGridRefVar = netcdf.defVar(ncid,'eDistEgridRef','float',[dimEdistBinsRef]);
phiGridVar = netcdf.defVar(ncid,'phiGrid','float',[dimAdistBins]);
thetaGridVar = netcdf.defVar(ncid,'thetaGrid','float',[dimAdistBins]);

netcdf.endDef(ncid);
 
netcdf.putVar(ncid, EVar, E);
netcdf.putVar(ncid, AVar, A);
netcdf.putVar(ncid,spyldVar,spyld);
netcdf.putVar(ncid,rfyldVar,rfyld);
netcdf.putVar(ncid,cosXdistVar,phiDists);
netcdf.putVar(ncid,cosYdistVar,thetaDists);

netcdf.putVar(ncid,cosXdistRefVar,phiDists);
netcdf.putVar(ncid,cosYdistRefVar,thetaDists);

netcdf.putVar(ncid, edistVar, eDists);
netcdf.putVar(ncid, edistRefVar, eDists);
netcdf.putVar(ncid,eDistGridVar,energy);
netcdf.putVar(ncid,eDistGridRefVar,energy);

netcdf.putVar(ncid, phiGridVar, phiGrid);
netcdf.putVar(ncid, thetaGridVar, thetaGrid);

netcdf.close(ncid);
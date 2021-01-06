close all
clear all

file = 'ftridynBackgroundDSi.nc';
ncid = netcdf.open(file,'NC_NOWRITE');
[dimname, nE] = netcdf.inqDim(ncid,0);
[dimname, nA] = netcdf.inqDim(ncid,1);
if strcmp(file,'ftridynBackground.nc')
[dimname, nS] = netcdf.inqDim(ncid,2);
else
    nS = 1;
end
energy = ncread(file,'E');
angle = ncread(file,'A');
spyld = ncread(file,'spyld');
rfyld = ncread(file,'rfyld');
cosxDist = ncread(file,'cosXDist');
cosxDistRef = ncread(file,'cosXDistRef');
cosyDist = ncread(file,'cosYDist');
% coszDist = ncread(file,'cosZDist');
eDist = ncread(file,'energyDist');
% eDistRef = ncread(file,'energyDistRef');
eDistEgrid = ncread(file,'eDistEgrid');
phiGrid = ncread(file,'phiGrid');
thetaGrid = ncread(file,'thetaGrid');
% eDistEgridRef = ncread(file,'eDistEgridRef');
% thisEdistRef = reshape(eDistRef(:,1,:),500,[]);
% figure(100)
% plot(eDistEgridRef,thisEdistRef)
% 
% thisEdistRef = eDistRef(:,:,10)
% figure(101)
% plot(eDistEgridRef,thisEdistRef)
if( nS > 1)
    for i=1:nS
    this_spyld = spyld(:,:,i)
    figure(i)
h = pcolor(energy,angle,log10(this_spyld))
h.EdgeColor = 'none';
colorbar
 set(gca, 'XScale', 'log')
    end
else
        this_spyld = spyld
    figure(1)
h = pcolor(energy,angle,this_spyld)
h.EdgeColor = 'none';
colorbar
 set(gca, 'XScale', 'log')
 set(gca, 'ColorScale', 'log')
  xlabel('Energy [eV]')
 ylabel('Angle [degrees]')

title({'Sputtering Yield','D on Si'})
 figure(2)
 semilogy(energy,this_spyld(1,:),'LineWidth',2)
 hold on
 for i=2:length(angle)
     semilogy(energy,this_spyld(i,:),'LineWidth',2)
 end
 xlabel('Energy [eV]')
 ylabel('Yield')
legend(cellstr(num2str(energy)))
title({'Sputtering Yield','D on Si'})
 legend
 axis([0 max(energy) 2e-4 3e-1])
end
figure(10)
plot(eDistEgrid,eDist(:,10,43,1))
figure(11)
plot(phiGrid,cosxDist(:,10,43,1))
figure(12)
plot(thetaGrid,cosyDist(:,35,50,1))

Y0=interpn(energy,angle,spyld',250,0);
R0=interpn(energy,angle,rfyld',250,0);
yr = Y0+R0;
Y0/yr
R0/yr
figure(113)
h = pcolor(energy,angle,log10(spyld))
h.EdgeColor = 'none';
colorbar
% set(gca, 'YDir', 'normal')
 set(gca, 'XScale', 'log')
title({'Sputtering Yield Kr on W','As a Function of Energy and Angle'})
xlabel('E [eV]') % x-axis label
ylabel('Angle [degrees]') % y-axis label
set(gca,'fontsize',16)
ax = h.Parent;   % Important
set(ax, 'XTick', [ 0 50 100 200 400 600 1000])
set(gca,'TickDir','out');
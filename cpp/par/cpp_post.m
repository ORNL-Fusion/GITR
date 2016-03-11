clear variables

Deposition
Depo = Dep;
figure(1)
imagesc(log10(fliplr(Depo)))
colorbar
caxis([-1 max(log10(Depo(:)))])
colormap jet

Charge
Charges = Dep;
MQ = Charges./Depo;
figure(2)
imagesc((fliplr(MQ)))
colorbar
caxis([0 max((MQ(:)))])
colormap jet

Energy
Energies = Dep;
ME = Energies./Depo;
figure(3)
imagesc((fliplr(ME)))
colorbar
caxis([-1 max((ME(:)))])
colormap jet

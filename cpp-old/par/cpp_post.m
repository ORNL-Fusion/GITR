clear variables

Deposition
Depo = Dep;
figure(31)
imagesc(linspace(-0.03,0.03,150),linspace(-0.03,0.03,150),log10(fliplr(Depo)))
colorbar
caxis([-1 max(log10(Depo(:)))])
colormap jet

Charge
Charges = Dep;
gMQ = Charges./Depo;
figure(32)
imagesc((fliplr(gMQ)))
colorbar
caxis([0 max((gMQ(:)))])
colormap jet

Energy
Energies = Dep;
ME = Energies./Depo;
ME(isnan(ME)) = 0;
figure(33)
imagesc((fliplr(ME)))
colorbar
caxis([0 max((ME(:)))])
colormap jet

Tallys = Depo;
gitr_MQ = gMQ;
gitr_ME = ME;
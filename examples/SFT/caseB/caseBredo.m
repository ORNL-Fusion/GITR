close all
clear all

nT = 700;
TD0 = linspace(1,1300,nT);
nS = 1000;
s = linspace(0,20,nS);
n0 = 1e20;
sinj = 0.15;
sv = 1.2;
fcond = 0.1;
Mi = 0.1;
k0e=2000;
gamma=7;
mD = 2;
mz=12;
Z=4;
phi_in = 1.73e23;

mu = mz/(mz+mD);
Beta_i = 3*(mu+5*sqrt(2)*Z^2*(1.1*mu^(5/2) - 0.35*mu^(3/2)) - 1)/(2.6-2*mu + 5.4*mu^2)

cs0=sqrt(2*TD0*1.602e-19/mD/1.66e-27);
P=gamma*n0*cs0.*TD0.*1.602e-19;
T_s = zeros(nT,nS);
dTids = zeros(nT,nS);
for i=1:nT
T_s(i,:) = TD0(i)*(1 + 7/2*fcond*P(i)*s./k0e./TD0(i)^(7/2)).^(2/7);
end
css=sqrt(2*T_s*1.602e-19/mD/1.66e-27);
for i=1:nT
a = 7/2*fcond*P(i)/k0e./TD0(i).^(7/2);
dTids(i,:) = TD0(i).*(2*a./7./(a*s+1).^(5/7));
end

tau_s0 = mz*TD0.*(TD0/2).^0.5./(6.8e4*(1+mD/mz)*n0/1e18*(Z^2)*15);
tau_s = mz*T_s.*(T_s/2).^0.5./(6.8e4*(1+mD/mz)*n0/1e18*(Z^2)*15);
dTds_sv = zeros(1,nT);
dTds_sinj = zeros(1,nT);
tau_sv = zeros(1,nT);
tau_sinj = zeros(1,nT);
for i=1:nT
dTds_sv(i) = interp1(s,dTids(i,:),sv);
dTds_sinj(i) = interp1(s,dTids(i,:),sinj);
tau_sv(i) = interp1(s,tau_s(i,:),sv);
tau_sinj(i) = interp1(s,tau_s(i,:),sinj);
end
Ts_sv = zeros(1,nT);
Ts_sinj = zeros(1,nT);
for i=1:nT
Ts_sv(i) = interp1(s,T_s(i,:),sv);
Ts_sinj(i) = interp1(s,T_s(i,:),sinj);
end
vTh = sqrt(TD0*1.602e-19./mz./1.66e-27);
vdiff = vTh.*vTh.*tau_s0./sinj;
vB = -Mi*cs0;

vTi = 1.602e-19/1.66e-27/mz*Beta_i*tau_s0.*dTds_sinj;

vpl = vdiff - vB + vTi;
np_pl = phi_in./vpl;
np_th = phi_in./vTh;
figure(1)
plot(TD0,np_pl./1e20,'lineWidth',2)
hold on
plot(TD0,np_th./1e20,'lineWidth',2)
% title('Case B Peak Impurity Density')
xlabel('T_{i0} [eV]')
ylabel('n_p [10^{20} m^{-3}]')
axis([0 150 0 1])
set(gca,'fontsize',16)
% M = csvread('../SFT-2/caseB-2-Table 1.csv',1,0)
% TGITR = M(:,1);
% nGITR = M(:,5);
% % nGITR_Scatt = M(:,5);
% zeroVal = find(nGITR == 0);
% % nGITR_noScatt(zeroVal) = [];
% nGITR(zeroVal) = [];
% TGITR(zeroVal) = [];
caseBnp
npGITR(1:2) = npGITR(1:2)/10;
hold on
s1=scatter(TGITR,npGITR/npGITR(1)*.4)
s1.MarkerEdgeColor = 'k';
legend('\phi_{in}/v_{prompt loss}','\phi_{in}/v_{thermal}','GITR')
tickT = linspace(0,150,7);
tau_sTick = mz*tickT.*(tickT/2).^0.5./(6.8e4*(1+mD/mz)*n0/1e18*(Z^2)*15);
vThTick = sqrt(tickT*1.602e-19./mz./1.66e-27);
xticks(tickT);
ax1 = gca;
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none','YColor','none','fontsize',16);
ax2.XLabel.String = '\lambda_{mfp}/s_{inj}';
% line(vTh,0*vTh,'Parent',ax2,'Color','k')
ax2.XLim = [0 tickT(end)]
% ax2.YLim = [-2 -1]
set(ax2,'XTick',tickT,'XTickLabel',round(vThTick.*tau_sTick./sinj,2))

pbaspect([1 1.25 1])
tau_sv = zeros(1,nT);
tau_sinj = zeros(1,nT);
for i=1:nT
tau_sv(i) = interp1(s,tau_s(i,:),sv);
tau_sinj(i) = interp1(s,tau_s(i,:),sinj);
end

% Ts_sv = Ti0.*(1+7/2*fcond*P*sv/k0e./Ti0.^(7/2)).^(2/7);
% Ts_sinj = Ti0.*(1+7/2*fcond*P*sinj/k0e./Ti0.^(7/2)).^(2/7);
Ts_sv = zeros(1,nT);
cs_sv = zeros(1,nT);
Ts_sinj = zeros(1,nT);
for i=1:nT
Ts_sv(i) = interp1(s,T_s(i,:),sv);
cs_sv(i) = interp1(s,css(i,:),sv);
Ts_sinj(i) = interp1(s,T_s(i,:),sinj);
end
FFf_sv = mz*1.66e-27*(-Mi*cs0)./Ts_sv./1.602e-19./tau_sv.*(sv-sinj);
FFf_csv = mz*1.66e-27*(-Mi*cs_sv)./Ts_sv./1.602e-19./tau_sv.*(sv-sinj);
FFf_sinj = mz*1.66e-27*(-Mi*cs0)./Ts_sinj./1.602e-19./tau_sinj.*(sv-sinj);
FFf_s0 = mz*1.66e-27*(-Mi*cs0)./TD0./1.602e-19./tau_s0.*(sv-sinj);
FiGf = (Beta_i-1)*log(Ts_sv./Ts_sinj);
nz_np = exp(FFf_s0+FiGf);
figure(2)
% semilogy(1./(TD0.^2)*1e4,nz_np.*np_pl.*vTi./1.73e23)
hold on
% semilogy(1./(TD0.^2)*1e4,nz_np.*np_pl.*vTh./1.73e23)
axis([0 3.5 1e-5 1e-2])

xlabel('(T_0 [eV])^{-2} x 10^{-4}')
ylabel('\phi_{leak}/\phi_{in}')
title('Divertor Leakage as a Function of T^{-2}')
set(gca,'fontsize',16)
pbaspect([1 1.25 1])
legend('Simple Fluid Theory','GITR')
FFfI = zeros(1,nT);
FiGfI = zeros(1,nT);
ds = s(2) - s(1);
intInds = find(s >= sinj & s<= sv);
for i=1:nT
FFf0 = mz*1.66e-27*(-Mi*cs0(i))./T_s(i,intInds)./1.602e-19./tau_s(i,intInds);
FFfI(i) = sum(FFf0)*ds;
FiGf0 = (Beta_i-1)*dTids(i,intInds)./T_s(i,intInds);
FiGfI(i) = sum(FiGf0)*ds;
end
nz_np = exp(FFfI+FiGfI);
% semilogy(1./(TD0.^2)*1e4,nz_np.*np_pl.*vTi./1.73e23)

% GITR_leak = M(:,7);
% TGITR = M(:,1);
% zeroVal = find(GITR_leak == 0);
% GITR_leak(zeroVal) = [];
% TGITR(zeroVal) = [];
% hold on
x = linspace(0,5e-4);
semilogy(1e4*x,0.03*exp(-2.65e4*x),'lineWidth',2)
set(gca, 'YScale', 'log')
caseBleakage
s1=scatter(1./TGITR./TGITR*1e4,leakage)
s1.MarkerEdgeColor = 'k';
x = linspace(0,5e-4);
% semilogy(1e4*x,0.03*exp(-2.65e4*x))
% set(gca, 'YScale', 'log')
legend('Simple Fluid Theory','GITR')
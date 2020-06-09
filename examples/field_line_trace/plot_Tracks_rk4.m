clear all
close all
load('D3DrzZ2D2Tiles10Res.mat')
rGeom = r;
zGeom = z;
Zgeom = Z;

circle_radius = 0.4;
angle = linspace(0,2*pi,1e3);
cr = circle_radius*sin(angle)+1.6;
cz = circle_radius*cos(angle);

file = 'output/history.nc';
x = ncread(file,'x');
y = ncread(file,'y');
z = ncread(file,'z');
vx = ncread(file,'vx');
vy = ncread(file,'vy');
vz = ncread(file,'vz');
charge = ncread(file,'charge');
weight = ncread(file,'weight');
sizeArray = size(x);
nP = sizeArray(2);


r = sqrt(x.^2+y.^2);
figure(16)
hold on
dt = 1e-8;
nT = 1e7;
Sub = 1;
tloop = 7850*1e-6;
nloops = floor(dt*nT/tloop);
timerange = floor(nT/nloops/Sub);
r_vec = zeros(1,nloops);
if nP ==1
for j=1:nloops
for i=1:nP
    trange = (j-1)*timerange+1:j*timerange;
plot(r(trange),z(trange))

end
zdiff = 0.01;
rsub = r(trange);
averager = mean(rsub(find((z(trange) < zdiff & z(trange) > -zdiff) & (r(trange) >1.6) )));
r_vec(j) = averager;
end
else
    
for j=1:nloops
for i=1:nP
    trange = (j-1)*timerange+1:j*timerange;
plot(r(trange,i),z(trange,i))
end
end
end

plot(rGeom,zGeom)
hold on
scatter(rGeom,zGeom)
plot(cr,cz,'k')

thisr_vec = r_vec;

names = {'boris_dt1en6_nt1e5_nocyl.mat','boris_dt1en7_nt1e6_nocyl.mat','boris_dt1en8_nt1e7_nocyl.mat'};
figure(1)
hold on
%clipboard('copy',r_vec)
r_vec_dt1en8_rk4=[1.99933513998014 1.99934235483463 1.99933549225843 1.99933280035192 1.99935479861583 1.99929134185881 1.99928930100024 1.9992206650752 1.99922681225799 1.99920440811411 1.99920003193382 1.99918549559521];
r_vec_dt1en7_rk4=[1.99929567035712 1.99922609523692 1.99914896779418 1.99910366195667 1.99904595491167 1.99894796357594 1.99888181570288 1.99877534395704 1.99868107536646 1.99859770080908 1.99852569010651 1.99845360241215];
r_vec_dt1en6_rk4=[1.9990080691153 1.99834015634325 1.99766612052917 1.99699633929037 1.9963258658686 1.99565448684077 1.99473218162461 1.99381760663764 1.99296294696747 1.99229059488543 1.99161315348841 1.99094428939204];
r_vec_dt1en8_rk4_cyl=[2.03177533438376 2.03177946632221 2.03180259748218 2.03186440212532 2.03182841303089 2.03187384827021 2.03185861353191 2.03192526810856 2.03191849499467 2.03196317768749 2.03193145468375 2.03193040285011];
r_vec_dt1en7_rk4_cyl=[2.03193646814961 2.03208140527387 2.03230559951492, 2.03245791540609 2.03264615711788 2.03283877628675 2.03298892368748 2.03316162292011 2.03333916244615 2.03352093601421 2.03371536874791 2.03386854080532];
r_vec_dt1en6_rk4_cyl=[2.0327049443161 2.03446040790441 2.0362331066184 2.03802901197962 2.039847186323 2.04168386382017 2.04354837296075 2.0454322235352 2.04733800475899 2.0492657682334 2.05122354603294 2.05320057227651];
% plot(thisr_vec-thisr_vec(1))
hold on
plot(r_vec_dt1en8_rk4-r_vec_dt1en8_rk4(1),'-o')
plot(r_vec_dt1en7_rk4-r_vec_dt1en7_rk4(1),'-o')
plot(r_vec_dt1en6_rk4-r_vec_dt1en6_rk4(1),'-o')
plot(r_vec_dt1en8_rk4_cyl-r_vec_dt1en8_rk4_cyl(1),'-*')
plot(r_vec_dt1en7_rk4_cyl-r_vec_dt1en7_rk4_cyl(1),'-*')
plot(r_vec_dt1en6_rk4_cyl-r_vec_dt1en6_rk4_cyl(1),'-*')
% for i=1:length(names)
%     load(names{i});
%     plot(r_vec-r_vec(1))
% end
legend('Straight dt=1e-8','Straight dt=1e-7','Straight dt=1e-6','Curved dt=1e-8','Curved dt=1e-7','Curved dt=1e-6')
% 
title({'Relative Outward Drift of Particles in Helical Field','In Straight and Cylindrical Geometries For RK4 Integrator'})
xlabel('Poloidal Loop Number')
ylabel('$\bar{r}_i - \bar{r}_1$ [m]','Interpreter','latex')
set(gca,'FontSize',14)
clear all
close all

nu=3e3;
nu_E=100;
load('nu.mat')
dt = 1e-6;
nP=1e3;

nT = 20000;
vy = zeros(1,nP);
vx = zeros(1,nP);
E=10;
vz = sqrt(2*E*1.602e-19/(184*1.66e-27))*ones(1,nP);
v=vz;
perp_direction1 = zeros(nP,3);
perp_direction2 = zeros(nP,3);

vColl = zeros(nP,3);
    v_par = [vx',vy',vz'];
tic
B_unit = [0,0,1];
for i=1:nT

    v0 = sqrt(v_par(:,1).^2+v_par(:,2).^2+v_par(:,3).^2);
    E = 0.5*184*1.66e-27/1.602e-19*v0.^2;
    parallel_direction = v_par./v0;
    s1 = parallel_direction(:,3);
    s2 = sqrt(1-s1.*s1);
    perp_direction1(:,1) = 1.0./s2.*(s1.*parallel_direction(:,1) - B_unit(1));
    perp_direction1(:,2) = 1.0./s2.*(s1.*parallel_direction(:,2)- B_unit(2));
    perp_direction1(:,3) = 1.0./s2.*(s1.*parallel_direction(:,3)- B_unit(3));
    perp_direction2(:,1) = 1.0./s2.*(B_unit(3)*parallel_direction(:,2));% - parallel_direction[2]*B_unit[1]);
    perp_direction2(:,2) = 1.0./s2.*( - B_unit(3)*parallel_direction(:,1));
    s20ind = find(s2==0);
    perp_direction1(s20ind,1) =  s1(s20ind);
    perp_direction1(s20ind,2) =  s2(s20ind);
    perp_direction1(s20ind,3) = 0;
    perp_direction2(s20ind,1) = parallel_direction(s20ind,3);
    perp_direction2(s20ind,2) = parallel_direction(s20ind,1);
    perp_direction2(s20ind,3) = parallel_direction(s20ind,2);
    s1(s20ind) = parallel_direction(s20ind,1).*perp_direction1(s20ind,1)+parallel_direction(s20ind,2).*perp_direction1(s20ind,2)+parallel_direction(s20ind,3).*perp_direction1(s20ind,3);
    s2(s20ind) = sqrt(1-s1(s20ind).*s1(s20ind));
    perp_direction1(s20ind,1) = -1.0./s2(s20ind).*(parallel_direction(s20ind,2).*perp_direction2(s20ind,3)-parallel_direction(s20ind,3).*perp_direction2(s20ind,2));
    perp_direction1(s20ind,2) = -1.0./s2(s20ind).*(parallel_direction(s20ind,3).*perp_direction2(s20ind,1)-parallel_direction(s20ind,1).*perp_direction2(s20ind,3));
    perp_direction1(s20ind,3) = -1.0./s2(s20ind).*(parallel_direction(s20ind,1).*perp_direction2(s20ind,2)-parallel_direction(s20ind,2).*perp_direction2(s20ind,1));
    
    nu_d0 = interp1(Eparticle,nu_d,E);
    nu_E0 = interp1(Eparticle,nu_E,E);
    
    n2 = abs(normrnd(0,1.0,1,nP));
    xsi = rand(1,nP);
    n3 = abs(normrnd(0,1.0,1,nP));
    cosXsi = cos(2*pi*xsi);
    sinXsi = sin(2*pi*xsi);
    vColl(:,1) = (parallel_direction(:,1) + sqrt(dt*0.5*nu_d0).*n2'.*(perp_direction1(:,1).*cosXsi' + perp_direction2(:,1).*sinXsi'));
    vColl(:,2) = (parallel_direction(:,2) + sqrt(dt*0.5*nu_d0).*n2'.*(perp_direction1(:,2).*cosXsi' + perp_direction2(:,2).*sinXsi'));
    vColl(:,3) = (parallel_direction(:,3) + sqrt(dt*0.5*nu_d0).*n2'.*(perp_direction1(:,3).*cosXsi' + perp_direction2(:,3).*sinXsi'));
    vCollMag = sqrt(vColl(:,1).^2+vColl(:,2).^2+vColl(:,3).^2);
%     vNew = (1-dt*0.5*nu_E0'.*n3)'.*v0./vCollMag.*vColl;
      vNew = (1-dt*0.5*nu_E0).*v0.*vColl;
    v_par = vNew;
    
end
toc
vx = v_par(:,1);
vy = v_par(:,2);
vz = v_par(:,3);
figure(1)
histogram(vx)
figure(2)
histogram(vy)
figure(3)
histogram(vz)

E = 0.5*184*1.66e-27/1.602e-19*(vx.^2+vy.^2+vz.^2);
figure(4)
histogram(E)
v=sqrt(vx.^2+vy.^2+vz.^2);
figure(5)
histogram(v)
averageV = mean(v);
energy = 0.5*184*1.66e-27*averageV^2/1.602e-19;

%Define some constants
e0 = 8.85e-12;
me = 9.11e-31;
mi = 1.66e-27;
q = 1.602e-19;


%Background Plasma parameters
mb = mi;
Lc = 10; %Connection length
width = 0.05; %SOL width
%Define Fields
B0 = 1;
B = B0*[0 -.141 .99]; %Magnitude multiplied by a normalized vector
Tu = 25;
Tt = 1;
shtc_e = 7;
shtc_i = 3.5;
n0 = 1e19;
nt = n0/2;
cst = sqrt(2*Tt*q/mb);
k0e = 1.8e3; %W/mev^-7/2 heat conductivities for hydrogenic species
k0i = 60;
D_perp = 0.2;
Z = 1;
cs0 = sqrt((2*Tu)*q/mb);

q0e = shtc_e*nt*cst*Tt*q;
q0i = shtc_e*nt*cst*Tt*q;




%Test particle parameters
ion = 0; %ionization state 1 = true, 0=false
m = 10*mi;
wc = q/m*B0;
Z_test = 1;
    mu = mb/(mb+m);
    alpha = Z_test^2*0.71;
    beta = -3*(1- mu - 5*sqrt(2)*Z_test^2)*(1.1*mu^(5/2) - 0.35*mu^(3/2))/(2.6 - 2*mu+ 5.4*mu^2);
file = '~/GitHub/gimp/matlab/ionization_rates/eH.txt';

fileID = fopen(file, 'r');

C = textscan(fileID, '%f32 %f32 %f32 %f32');

C = [C{1}, C{2}];

fclose('all');
loglog(C(:,1),C(:,2))

%Run specific parameters
Nt = 1e6;
dt = 2/wc/1e2;

%Initialize Energy-Angle for particle
%Ka = [1, pi/4,pi/4]; 
% Energy (eV), azimuthal angle (phi - 0 to pi), polar angle(theta - 0 to pi)


%Set initial position
%r = [0,0,0]; %(x,y,z), y=0 indicates the particle is on the surface







        
%History 
r_hist = zeros(Nt+1,3);
r_hist(1,:) = r;
        v0 = sqrt(Ka(1)*2*q/m);
        v = zeros(1,3);
        v_minus = zeros(1,3);
        v(1) = v0*sin(Ka(2))*cos(Ka(3));
        v(2) = v0*sin(Ka(2))*sin(Ka(3));
        v(3) = v0*cos(Ka(2));
        
        
          %Main Loop
for j=1:Nt
    if ion ==0
        %Neutral particle mover
        r = r + v*dt;

        % Ionization module
        %ioniz
        
        if (rand > 0.99)
           ion = 1;
        end


        
        
    else %particle mover for ions
        v0 = norm(v);

Ero_T
Ero_E
Ero_Relax_max
E_ext
        %Find parallel and perp velocity directions
        direct
 
dv = dt/m*(Z_test*q*([E_par E_perp 0]+ Eext + cross(v,B)) ...
    + alpha*[grad_Ts_e grad_Tr_e 0] + beta*[grad_Ts_i grad_Tr_i 0]) ...
    + v0/tau_s*dt*dir1  ...
            + randn*sqrt(v0^2/tau_e*dt)*dir1 ...
            + randn*sqrt(1/2*v0^2/tau_d*dt)*dir2 ...
            + randn*sqrt(1/2*v0^2/tau_d*dt)*dir3;

        v = v+ dv;
        r = r + v*dt - sqrt(D_perp*dt)*[0 1 0];

    end

     r_hist(j+1,:) = r;
         if r(2) < 0 
        break
         end

end
    
    %Plot position history
figure(1)
plot3(r_hist(1:j,1),r_hist(1:j,3),r_hist(1:j,2))
grid on
set(gca, 'YDir', 'reverse')
axis equal
xlabel('x')
ylabel('z')
zlabel('y')
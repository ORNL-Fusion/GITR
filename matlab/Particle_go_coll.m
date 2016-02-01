clear variables
%Initialize Energy-Angle for particle
Ka = [1, pi/4,pi/4]; 
% Energy (eV), azimuthal angle (phi - 0 to pi), polar angle(theta - 0 to pi)


%Set initial position
r = [0,0,0]; %(x,y,z), y=0 indicates the particle is on the surface
Lc = 10; %Connection length

ion = 0; %ionization state 1 = true, 0=false

%Define Fields
B0 = 1;
B = B0*[0 -1/sqrt(2) 1/sqrt(2)]; %Magnitude multiplied by a normalized vector

%Background plasma parameters
Te = 25;
D_perp = 1e-4;

%Set constants and time steps
Nt = 1e6;
q = 1.602e-19;
m = 1.66e-27;
wc = q/m*B0;
dt = 2/wc/1e3

        %Constants used in Boris method Lorentz Integrator
        q_prime = q/m*dt/2;
        coeff = 2*q_prime/(1+(q_prime*B0)^2);
        
%History 
r_hist = zeros(Nt+1,3);
r_hist(1,:) = r;
        v0 = sqrt(Ka(1)*2*q/m);
        v = zeros(1,3);
        v_minus = zeros(1,3);
        v(1) = v0*sin(Ka(2))*cos(Ka(3));
        v(2) = v0*sin(Ka(2))*sin(Ka(3));
        v(3) = v0*cos(Ka(2));
       
%Calculate relaxation times for the given particle,
% velocity and background plasma
relax_time_max
        
        %Main Loop
for j=1:Nt
    if ion ==0
        %Neutral particle mover
        r = r + v*dt;

        % Made up ionization
        s = rand(1);
        if rand > 0.9999
            ion = 1;
        end
        
        
    else %particle mover for ions
        %update fields
        E = [0 -1e3 1/Lc*1e3];
v0 = norm(v);
        %Find parallel and perp velocity directions
        direct
        %Boris Method Lorentz Integrator
        v_minus = v + q_prime*E + v0/tau_s*dt/2*dir1  ...
            + randn*sqrt(v0^2/tau_e*dt/2)*dir1 ...
            + randn*sqrt(1/2*v0^2/tau_d*dt/2)*dir2 ...
            + randn*sqrt(1/2*v0^2/tau_d*dt/2)*dir3;

        v = v_minus + q_prime*cross(v_minus,B);
        
        v = v_minus + coeff*cross(v,B);
        
        v = v + q_prime*E + v0/tau_s*dt/2*dir1  ...
            + randn*sqrt(v0^2/tau_e*dt/2)*dir1 ...
            + randn*sqrt(1/2*v0^2/tau_d*dt/2)*dir2 ...
            + randn*sqrt(1/2*v0^2/tau_d*dt/2)*dir3;
     
        r = r + v*dt - sqrt(D_perp*dt)*[0 1 0];


    end
    
    r_hist(j+1,:) = r;
    
    %Break Loop if particle strikes material surface
    if r(2) < 0 
        break
    end
end

%Plot position history
figure(1)
plot3(r_hist(1:j,1),r_hist(1:j,3),r_hist(1:j,2))
grid on
set(gca, 'YDir', 'reverse')

xlabel('x')
ylabel('z')
zlabel('y')
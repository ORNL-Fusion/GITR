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
        
        ion = 0;
          %Main Loop
for j=1:Nt
    if ion ==0
        %Neutral particle mover
        r = r + v*dt;

        % Ionization module
        %ioniz
        
        if (rand > 0.9)
           ion = 1;
        end


        
        
    else %particle mover for ions
        v0 = norm(v);
r1 = isnan(r(1));
r2 = isnan(r(2));
r3 = isnan(r(3));
if ((r1 >0) | (r2>0) | (r3>0) )
    r = [-1e-5,-1e-5,-1e-5];
    break
end


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
%{     
     figure(1)

plot3(r_hist(1:j,1),r_hist(1:j,3),r_hist(1:j,2))
grid on
set(gca, 'YDir', 'reverse')
axis equal
xlabel('x')
ylabel('z')
zlabel('y')
pause(0.1)
%}
    if ((r(1) >= Npol*Lp) | (r(3)>= Ntor*Lt) )
    r
    r = [-1e-5,-1e-5,-2e-5];
    break
         if r(2) < 0 
        break
         end

         end
if ((r(1) <0) | (r(2)<0) | (r(3)<0) )
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

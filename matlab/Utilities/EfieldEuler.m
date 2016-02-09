            ME = 9.10938356e-31;
            MI = 1.6737236e-27;
            Q = 1.60217662e-19;
            EPS0 = 8.854187e-12;
            
Z = 2;
amu = 184;
xp = [ 0 0 5e-3];
vp = [0 0 -3.2e3];
dt = 1e-8;
tspan = 1e-5;
ld = 1e-5;
B = [0 0 -2];
BMagPart = norm(B);

surface_dz_dx = 1.73205;
surface_zIntercept = 0;
i = 0;
maxStep = 0;
dvMax = 0;
dv_threshold = 0.001;
perpDistanceToSurface = ( -surface_dz_dx*xp(1) + xp(3) + surface_zIntercept)/sqrt(surface_dz_dx^2+1);
surfaceDirection = [surface_dz_dx 0 -1];
surfaceDirection_unit = surfaceDirection/norm(surfaceDirection);
time = 0;
tolerance = 1e-16;
tic
while abs(tspan - time) > tolerance && perpDistanceToSurface > 0
   
        E = 45/(2*ld)*exp(-perpDistanceToSurface/(2*ld)) * surfaceDirection_unit;

        % Constants used in Boris method Lorentz Integrator
        q_prime = Z*Q/(amu*MI)*dt/2;
        coeff = 2*q_prime/(1+(q_prime*BMagPart).^2);
        
        % Boris Method Lorentz Integrator
        v = vp;
        r = xp;
        v_start = v;
        v_minus = v + q_prime*E;
        
        v = v_minus + q_prime*[v_minus(2)*B(3) - v_minus(3)*B(2), v_minus(3)*B(1) - v_minus(1)*B(3),v_minus(1)*B(2) - v_minus(2)*B(1)];
        
        v = v_minus + coeff*[v(2)*B(3) - v(3)*B(2), v(3)*B(1) - v(1)*B(3),v(1)*B(2) - v(2)*B(1)];
        
        v = v + q_prime*E;
        
        step = v*dt;
        r = r + step;
        
        if abs(norm(v) - norm(vp))/norm(vp) > dv_threshold
            dt = dt/10;
        else
            perpDistanceToSurface = ( -surface_dz_dx*r(1) + r(3) + surface_zIntercept)/sqrt(surface_dz_dx^2+1);
            if perpDistanceToSurface < 0
                disp('hit surface')
                time
                perpDistanceToSurface0 = ( -surface_dz_dx*xp(1) + xp(3) + surface_zIntercept)/sqrt(surface_dz_dx^2+1);
                t = (-surface_zIntercept + surface_dz_dx*xp(1) - xp(3))/(-surface_dz_dx*(xp(1) -r(1)) + (xp(3) -r(3)));
                xp = xp + (xp -r)*t;
                vp = vp + (vp -v)*perpDistanceToSurface0/(perpDistanceToSurface0-perpDistanceToSurface);
            else
            xp = r;
            vp = v;
            end
            time = time + dt;
        end
        %perpDistanceToSurface = ( -surface_dz_dx*xp(1) + xp(3) + surface_zIntercept)/sqrt(surface_dz_dx^2+1);
        i = i+1;
end
toc
0.5*184*1.66e-27*norm(v)^2/1.602e-19
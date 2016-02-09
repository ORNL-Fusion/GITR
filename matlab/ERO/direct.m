exist_v = exist('v');

if exist_v == 0
%Initialize Energy-Angle for particle
Ka = [1, pi/4,pi/4]; 

        v0 = sqrt(Ka(1)*2*q/m);
        v = zeros(1,3);
        v(1) = v0*sin(Ka(2))*cos(Ka(3));
        v(2) = v0*sin(Ka(2))*sin(Ka(3));
        v(3) = v0*cos(Ka(2));
        stop
end
        dir1 = v/norm(v);
        rand_dir = dir1;
        while rand_dir == dir1
            rand_dir = rand(1,3);
            rand_dir = rand_dir/norm(rand_dir);
        end
            dir2 = cross(dir1,rand_dir);
            dir2 = dir2/norm(dir2);
            
            dir3 = cross(dir1,dir2);
            dir3 = dir3/norm(dir3);

        
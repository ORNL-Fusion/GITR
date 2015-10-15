function dydt = myode(t,y,p,E,B,xyz)
    
    constants
    
    dydt = zeros(6,1);
    
    thisE = [0 0 0];
    thisB = [0 0 0];
         
%     thisB(1) = interpn(x,y,z,B.x,p.x,p.y,p.z);
%     thisB(2) = interpn(x,y,z,B.y,p.x,p.y,p.z);
%     thisB(3) = interpn(x,y,z,B.z,p.x,p.y,p.z);
   
    % these are 1-D interpolations since our profiles
    % have d/dy = d/dz =0

    thisE(1) = interpn(xyz.x,E.x(:,1,1),y(1));
    thisE(2) = interpn(xyz.x,E.y(:,1,1),y(1));
    thisE(3) = interpn(xyz.x,E.z(:,1,1),y(1));

    thisB(1) = interpn(xyz.x,B.x(:,1,1),y(1));
    thisB(2) = interpn(xyz.x,B.y(:,1,1),y(1));
    thisB(3) = interpn(xyz.x,B.z(:,1,1),y(1));
    
    if y(1) > 0
        dydt(:) = 0; % this is a hack to catch particles hitting the wall.
    else
        dydt(1:3) = y(4:6);    
        dydt(4:6) = (p.Z*Q/(p.amu*MI)) * (thisE + cross(y(4:6),thisB) );
    end
    
end


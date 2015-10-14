function [ dy ] = myode(t,v,p,B,xyz)
    
    constants
    
    dy = zeros(6,1);
    
    E = [0 0 0]';
    thisB = [1 0 0]';
    
    x=xyz.x;
    y=xyz.y;
    z=xyz.z;
    
%     thisB(1) = interpn(x,y,z,B.x,p.x,p.y,p.z);
%     thisB(2) = interpn(x,y,z,B.y,p.x,p.y,p.z);
%     thisB(3) = interpn(x,y,z,B.z,p.x,p.y,p.z);
%     
%      thisB(1) = interpn(x,B.x(:,1,1),v(1));
%      thisB(2) = interpn(x,B.y(:,1,1),v(1));
%      thisB(3) = interpn(x,B.z(:,1,1),v(1));
    
    dy(4:6) = (p.Z*q/(p.amu*mi)) * (E + cross(v(4:6),thisB) );
    dy(1:3) = v(1:3);
    
end


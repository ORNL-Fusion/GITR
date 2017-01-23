function OutputField = gitrInterpVector1D(particle,xyz,field,s,varargin)
dimensions = ndims(field.x);
x = xyz.x;
switch dimensions
    case 4
        OutputField(1) = interp1q(x,field.x(:,1,1,s),particle.x);
        OutputField(2) = interp1q(x,field.y(:,1,1,s),particle.x);
        OutputField(3) = interp1q(x,field.z(:,1,1,s),particle.x);
        
    case 3        
        OutputField(1) = interp1q(x,field.x(:,1,1),particle.x);
        OutputField(2) = interp1q(x,field.y(:,1,1),particle.x);
        OutputField(3) = interp1q(x,field.z(:,1,1),particle.x);
        
end


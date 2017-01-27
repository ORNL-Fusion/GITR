function OutputField = gitrInterpScalar1D(particle,xyz,field,s,varargin)
dimensions = ndims(field);
x = xyz.x;
switch dimensions
    case 4
        OutputField = interp1q(x,field(:,1,1,s),particle.x);
    case 3
    OutputField = interp1q(x,field(:,1,1),particle.x);

end


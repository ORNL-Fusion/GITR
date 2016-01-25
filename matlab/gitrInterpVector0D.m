function OutputField = gitrInterpVector0D(particle,xyz,field3D,s,varargin)
% nVarargs = length(varargin);
% n = nargin - nVarargs;

dimensions = ndims(field3D.x);
switch dimensions
    case 4
    OutputField(1) = field3D.x(1,1,1,s);
    OutputField(2) = field3D.y(1,1,1,s);
    OutputField(3) = field3D.z(1,1,1,s);
    
    case 3
    OutputField(1) = field3D.x(1,1,1);
    OutputField(2) = field3D.y(1,1,1);
    OutputField(3) = field3D.z(1,1,1);
end

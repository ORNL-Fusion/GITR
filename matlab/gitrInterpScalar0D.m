function Output = gitrInterpScalar0D(particle,xyz,field,s,varargin)
dimensions = ndims(field);
switch dimensions
    case 4
        Output = field(1,1,1,s);
    case 3
        Output = field(1,1,1);
end


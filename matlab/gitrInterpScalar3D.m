function OutputField = gitrInterpScalar3D(particle,xyz,field3D,s, varargin)
dimensions = ndims(field3D);
switch dimensions
    case 4
    OutputField = interpn(xyz.x,xyz.y,xyz.z,field3D(:,:,:,s),particle.x,particle.y,particle.z);
    case 3
    OutputField = interpn(xyz.x,xyz.y,xyz.z,field3D,particle.x,particle.y,particle.z);

end

function OutputField = gitrInterpVector3D(particle,xyz,field3D,s, varargin)
dimensions = ndims(field3D.x);
switch dimensions
    case 4

    OutputField(1) = interpn(xyz.x,xyz.y,xyz.z,field3D.x(:,:,:,s),particle.x,particle.y,particle.z);
    OutputField(2) = interpn(xyz.x,xyz.y,xyz.z,field3D.y(:,:,:,s),particle.x,particle.y,particle.z);
    OutputField(3) = interpn(xyz.x,xyz.y,xyz.z,field3D.z(:,:,:,s),particle.x,particle.y,particle.z);
    
    case 3
    OutputField(1) = interpn(xyz.x,xyz.y,xyz.z,field3D.x,particle.x,particle.y,particle.z);
    OutputField(2) = interpn(xyz.x,xyz.y,xyz.z,field3D.y,particle.x,particle.y,particle.z);
    OutputField(3) = interpn(xyz.x,xyz.y,xyz.z,field3D.z,particle.x,particle.y,particle.z);
end

function OutputField = gimpInterpVector3D(particle,xyz,field3D)

    OutputField(1) = interpn(xyz.x,xyz.y,xyz.z,field3D.x,particle.x,particle.y,particle.z);
    OutputField(2) = interpn(xyz.x,xyz.y,xyz.z,field3D.y,particle.x,particle.y,particle.z);
    OutputField(3) = interpn(xyz.x,xyz.y,xyz.z,field3D.z,particle.x,particle.y,particle.z);

end

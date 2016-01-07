function OutputField = gimpInterpScalar3D(particle,xyz,field3D)

    OutputField = interpn(xyz.x,xyz.y,xyz.z,field3D,particle.x,particle.y,particle.z);

end

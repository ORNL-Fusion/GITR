function OutputField = gitrInterp1D(particle,xyz,field)

    OutputField(1) = interp1(xyz.x,field.x(:,1,1),particle.x);
    OutputField(2) = interp1(xyz.x,field.y(:,1,1),particle.x);
    OutputField(3) = interp1(xyz.x,field.z(:,1,1),particle.x);

end


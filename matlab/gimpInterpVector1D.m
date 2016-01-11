function OutputField = gimpInterpVector1D(particle,xyz,field)

    x = xyz.x;

    OutputField(1) = interp1q(x,field.x(:,1,1),particle.x);
    OutputField(2) = interp1q(x,field.y(:,1,1),particle.x);
    OutputField(3) = interp1q(x,field.z(:,1,1),particle.x);

end


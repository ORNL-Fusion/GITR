function OutputField = gimpInterpVector1D(particle,xyz,field)

    x = xyz.x;
    y = field.x(:,1,1);

    OutputField(1) = interp1q(x,y,particle.x);
    OutputField(2) = interp1q(x,y,particle.x);
    OutputField(3) = interp1q(x,y,particle.x);

end


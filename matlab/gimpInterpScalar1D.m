function OutputField = gimpInterpScalar1D(particle,xyz,field)

    OutputField = interp1q(xyz.x,field(:,1,1),particle.x);

end


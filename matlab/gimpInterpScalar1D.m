function OutputField = gimpInterpScalar1D(particle,xyz,field)

    OutputField = interp1(xyz.x,field(:,1,1),particle.x);

end


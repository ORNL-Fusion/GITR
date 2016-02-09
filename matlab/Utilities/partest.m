dens = zeros(10,1);

parfor p=1:10
    for t=1:5
    dens(p+t*0) = p+t;
    end
end

dens
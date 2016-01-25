
if EfieldInterpolator_number ==0
    EfieldInterpolator = @gitrInterpVector0D;
else if EfieldInterpolator_number == 1
        EfieldInterpolator = @gitrInterpVector1D;
    else if EfieldInterpolator_number == 2
            EfieldInterpolator = @gitrEfieldAnalytic;
        else if EfieldInterpolator_number == 3
                EfieldInterpolator = @gitrInterpVector3D;
            end
        end
    end
end



if EfieldInterpolator_number == 0
    Efield3D.x(:) = Efield_in(1);
    Efield3D.y(:) = Efield_in(2);
    Efield3D.z(:) = Efield_in(3);
else if EfieldInterpolator_number == 1
        Efield3D.x = reshape(dlmread('Efield_x_out.txt','\t'),nXv, nYv, nZv);
        Efield3D.y = reshape(dlmread('Efield_y_out.txt','\t'),nXv, nYv, nZv);
        Efield3D.z = reshape(dlmread('Efield_z_out.txt','\t'),nXv, nYv, nZv);
    else if EfieldInterpolator_number == 2
            initPart = particle;
            for i=1:nXv
                for j=1:nYv
                    for k=1:nZv
                        initPart.x = xyz.x(i);
                        initPart.y = xyz.y(j);
                        initPart.z = xyz.z(k);
                        initPart.hitWall = 0;
                        initPart.leftVolume = 0;
                        initPart.PerpDistanceToSurface(surface_dz_dx,surface_zIntercept);
                        
                        if initPart.perpDistanceToSurface > 0
                            
                            EfieldInterpolator = interpolators{1};
                            potential = -3*temp_eV(i,j,k,1);
                            B = [Bfield3D.x(i,j,k) Bfield3D.y(i,j,k) Bfield3D.z(i,j,k)];
                            
                            E = EfieldInterpolator(initPart, debyeLength, potential,surface_dz_dx,xyz,B, ...
                                background_Z,background_amu,maxTemp_eV);
                            
                            if E(1) > 1e10
                                stop
                            end
                            
                            Efield3D.x(i,j,k) = E(1);
                            Efield3D.y(i,j,k) = E(2);
                            Efield3D.z(i,j,k) = E(3);
                        end
                        
                    end
                end
            end
            
            
            Efield3D.x(isinf(Efield3D.x)) = 0;
            Efield3D.y(isinf(Efield3D.y)) = 0;
            Efield3D.z(isinf(Efield3D.z)) = 0;
        else if EfieldInterpolator_number == 3
                Efield3D.x = reshape(dlmread('Efield_x_out.txt','\t'),nXv, nYv, nZv);
                Efield3D.y = reshape(dlmread('Efield_y_out.txt','\t'),nXv, nYv, nZv);
                Efield3D.z = reshape(dlmread('Efield_z_out.txt','\t'),nXv, nYv, nZv);
            end
        end
    end
end
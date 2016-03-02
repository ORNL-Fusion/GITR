
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
            if printProfiles
            initPart(nXv) = particle;
            fieldx = zeros(nXv,nYv,nZv);
            fieldy = zeros(nXv,nYv,nZv);
            fieldz = zeros(nXv,nYv,nZv);
            for i=1:nXv
                tmpx = zeros(nXv,nYv,nZv);
                tmpy = zeros(nXv,nYv,nZv);
                tmpz = zeros(nXv,nYv,nZv);
                for j=1:nYv
                    for k=1:nZv
                        initPart(i).x = xyz.x(i);
                        initPart(i).y = xyz.y(j);
                        initPart(i).z = xyz.z(k);
                        initPart(i).hitWall = 0;
                        initPart(i).leftVolume = 0;
                        initPart(i).PerpDistanceToSurface(surface_dz_dx,surface_zIntercept);
                        
                        if initPart(i).perpDistanceToSurface > 0
                            
                           
                            potential = -3*temp_eV(i,j,k,1);
                            B = [Bfield3D.x(i,j,k) Bfield3D.y(i,j,k) Bfield3D.z(i,j,k)];
                            
                            E = EfieldInterpolator(initPart(i),xyz,Efield3D, debyeLength, potential,surface_dz_dx,B, ...
                                background_Z,background_amu,maxTemp_eV);
                            
                            if E(1) > 1e10
                                stop
                            end
                            tmpx(i,j,k) = tmpx(i,j,k) +E(1);
                            tmpy(i,j,k) = tmpy(i,j,k) +E(2);
                            tmpz(i,j,k) = tmpz(i,j,k) +E(3);
                        end
                        
                    end
                end
                fieldx = fieldx + tmpx;
                fieldy = fieldy + tmpy;
                fieldz = fieldz + tmpz;
            end
            
            
            Efield3D.x(isinf(Efield3D.x)) = 0;
            Efield3D.y(isinf(Efield3D.y)) = 0;
            Efield3D.z(isinf(Efield3D.z)) = 0;
            
            Efield3D.x = fieldx;
            Efield3D.y = fieldy;
            Efield3D.z = fieldz;
            clear fieldx fieldy fieldz
            end
        else if EfieldInterpolator_number == 3
                Efield3D.x = reshape(dlmread('Efield_x_out.txt','\t'),nXv, nYv, nZv);
                Efield3D.y = reshape(dlmread('Efield_y_out.txt','\t'),nXv, nYv, nZv);
                Efield3D.z = reshape(dlmread('Efield_z_out.txt','\t'),nXv, nYv, nZv);
            end
        end
    end
end
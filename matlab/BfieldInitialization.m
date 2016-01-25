

if BfieldInterpolator_number ==0
    BfieldInterpolator = @gitrInterpVector0D;
else if BfieldInterpolator_number == 1
        BfieldInterpolator = @gitrInterpVector1D;
    else if BfieldInterpolator_number == 2
            BfieldInterpolator = @gitrBfieldAnalytic;
        else if BfieldInterpolator_number == 3
                BfieldInterpolator = @gitrInterpVector3D;
            end
        end
    end
end

if BfieldInterpolator_number == 0
    Bfield3D.x(:) = Bx_in; % Create background B strucutre
    Bfield3D.y(:) = By_in;
    Bfield3D.z(:) = Bz_in;  
    Bfield3D.mag(:) = sqrt( Bfield3D.x.^2 + Bfield3D.y.^2 + Bfield3D.z.^2 );
    
else if BfieldInterpolator_number == 1
        Bfield3D.x = reshape(dlmread('Bfield_x_out.txt','\t'),nXv, nYv, nZv);
        Bfield3D.y = reshape(dlmread('Bfield_y_out.txt','\t'),nXv, nYv, nZv);
        Bfield3D.z = reshape(dlmread('Bfield_z_out.txt','\t'),nXv, nYv, nZv);
        Bfield3D.mag = reshape(dlmread('Bfield_mag_out.txt','\t'),nXv, nYv, nZv);
        
    else if BfieldInterpolator_number == 2
            error('no analytical model for Bfield implemented')
            
        else if BfieldInterpolator_number == 3
                Bfield3D.x = reshape(dlmread('Bfield_x_out.txt','\t'),nXv, nYv, nZv);
                Bfield3D.y = reshape(dlmread('Bfield_y_out.txt','\t'),nXv, nYv, nZv);
                Bfield3D.z = reshape(dlmread('Bfield_z_out.txt','\t'),nXv, nYv, nZv);
                Bfield3D.mag = reshape(dlmread('Bfield_mag_out.txt','\t'),nXv, nYv, nZv);
                
            end
        end
    end
end
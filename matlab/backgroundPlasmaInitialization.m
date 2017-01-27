if temperatureInterpolator_number ==0
    temperatureInterpolator = @gitrInterpScalar0D;
else if temperatureInterpolator_number == 1
        temperatureInterpolator = @gitrInterpScalar1D;
    else if temperatureInterpolator_number == 2
            temperatureInterpolator = @gitrTemperatureAnalytic;
        else if temperatureInterpolator_number == 3
                temperatureInterpolator = @gitrInterpScalar3D;
            end
        end
    end
end

if densityInterpolator_number ==0
    densityInterpolator = @gitrInterpScalar0D;
else if densityInterpolator_number == 1
        densityInterpolator = @gitrInterpScalar1D;
    else if densityInterpolator_number == 2
            densityInterpolator = @gitrDensityAnalytic;
        else if densityInterpolator_number == 3
                densityInterpolator = @gitrInterpScalar3D;
            end
        end
    end
end

if FlowVelocityInterpolator_number ==0
    FlowVelocityInterpolator = @gitrInterpVector0D;
else if FlowVelocityInterpolator_number == 1
        FlowVelocityInterpolator = @gitrInterpVector1D;
    else if FlowVelocityInterpolator_number == 2
            FlowVelocityInterpolator = @gitrFlowVelocityAnalytic;
        else if FlowVelocityInterpolator_number == 3
                FlowVelocityInterpolator = @gitrInterpVector3D;
            end
        end
    end
end

if temperatureInterpolator_number == 0
    for s=1:nS
        temp_eV(:,:,:,s) = maxTemp_eV(s);
    end
else if temperatureInterpolator_number == 1
        temp_eV = reshape(dlmread('temp_eV_out.txt','\t'),nXv, nYv, nZv,nS);
    else if temperatureInterpolator_number == 2
            error('no analytical model for Temperature implemented')
        else if temperatureInterpolator_number == 3
                temp_eV = reshape(dlmread('temp_eV_out.txt','\t'),nXv, nYv, nZv,nS);
            end
        end
    end
end

if densityInterpolator_number == 0
    for s=1:nS
        density_m3(:,:,:,s) = maxDensity(s);
    end
else if densityInterpolator_number == 1
        density_m3 = reshape(dlmread('density_m3_out.txt','\t'),nXv, nYv, nZv,nS);
    else if densityInterpolator_number == 2
            error('no analytical model for density implemented')
        else if densityInterpolator_number == 3
                density_m3 = reshape(dlmread('density_m3_out.txt','\t'),nXv, nYv, nZv,nS);
            end
        end
    end
end

if FlowVelocityInterpolator_number == 0
    for s=1:nS
        thermalVelocity = 9823.8*sqrt((temp_eV(:,:,:,1)+temp_eV(:,:,:,2))/background_amu(2))...
            *background_flow(s);
        flowVelocity_ms.x(:,:,:,s) = thermalVelocity(1).*Bfield3D.x./Bfield3D.mag;
        flowVelocity_ms.y(:,:,:,s) = thermalVelocity(2).*Bfield3D.y./Bfield3D.mag;;
        flowVelocity_ms.z(:,:,:,s) = thermalVelocity(3).*Bfield3D.z./Bfield3D.mag;;
    end
else if FlowVelocityInterpolator_number == 1
        flowVelocity_ms.x = reshape(dlmread('flowVelocity_x_out.txt','\t'),nXv, nYv, nZv,nS);
        flowVelocity_ms.y = reshape(dlmread('flowVelocity_y_out.txt','\t'),nXv, nYv, nZv,nS);
        flowVelocity_ms.z = reshape(dlmread('flowVelocity_z_out.txt','\t'),nXv, nYv, nZv,nS);
    else if FlowVelocityInterpolator_number == 2
            if printProfiles
            initPart = particle;
            for i=1:nXv
                for j=1:nYv
                    for k=1:nZv
                        initPart.x = xyz.x(i);
                        initPart.y = xyz.y(j);
                        initPart.z = xyz.z(k);
                        
                        B_local = [Bfield3D.x(i,j,k) Bfield3D.y(i,j,k) Bfield3D.z(i,j,k)];
                        
                        
                        
            for s=1:nS

            flowVelocity = FlowVelocityInterpolator(initPart,xyz,flowVelocity_ms,s,B_local, ...
                connectionLength,temp_eV,background_amu,surface_dz_dx,surface_zIntercept);
            end 
                        flowVelocity_ms.x(i,j,k,:) = flowVelocity(1);
                        flowVelocity_ms.y(i,j,k,:) = flowVelocity(2);
                        flowVelocity_ms.z(i,j,k,:) = flowVelocity(3);
                    end
                end
            end
            end
        else if FlowVelocityInterpolator_number == 3
                flowVelocity_ms.x = reshape(dlmread('flowVelocity_x_out.txt','\t'),nXv, nYv, nZv,nS);
                flowVelocity_ms.y = reshape(dlmread('flowVelocity_y_out.txt','\t'),nXv, nYv, nZv,nS);
                flowVelocity_ms.z = reshape(dlmread('flowVelocity_z_out.txt','\t'),nXv, nYv, nZv,nS);
            end
        end
    end
end
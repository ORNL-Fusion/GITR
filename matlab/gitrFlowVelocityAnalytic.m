      function flowVelocity = gitrFlowVelocityAnalytic(this,xyz,flowVelocity_ms,species,Bfield,connectionLength, ...
                temp_eV,background_amu,surface_dz_dx,surface_zIntercept)
            
            B_unit = Bfield/norm(Bfield);
            thermalVelocity = 9823.8*sqrt((temp_eV(1,1,1,1) + temp_eV(1,1,1,2))/background_amu(2));
            
            surfaceNormal = [-surface_dz_dx 0 1];
            
            B_dot_sN = B_unit(1)*surfaceNormal(1)+B_unit(2)*surfaceNormal(2)+B_unit(3)*surfaceNormal(3);
            
            if B_dot_sN < 0
                plu_min = 1;
            else
                plu_min = -1;
            end
            
            B_unit = B_unit*plu_min;
            
            s = (surface_zIntercept - this.z + surface_dz_dx*this.x)/(B_unit(3) - surface_dz_dx*B_unit(1));
            s = connectionLength/2 - s;
            flowVelocity = thermalVelocity(1)*(connectionLength/(2*s)- sqrt((connectionLength/(2*s))^2 - 1))*B_unit;
            
            if s >= connectionLength/2
                flowVelocity = [0 0 0];
            end

        end
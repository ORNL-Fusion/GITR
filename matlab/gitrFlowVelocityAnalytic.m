        function flowVelocity = gitrFlowVelocityAnalytic(this,xyz,Bfield,connectionLength, ...
                maxTemp_eV,background_amu,surface_dz_dx,surface_zIntercept)
            
            B_unit = Bfield/norm(Bfield);
            thermalVelocity = 9823.8*sqrt((maxTemp_eV + maxTemp_eV)/background_amu(2));
            
            s = (surface_zIntercept - this.z + surface_dz_dx*this.x)/(B_unit(3) - surface_dz_dx*B_unit(1));
            s = connectionLength/2 - s;
            flowVelocity = thermalVelocity*(connectionLength/(2*s)- sqrt((connectionLength/(2*s))^2 - 1))*B_unit;
            
            if s >= connectionLength/2
                flowVelocity = [0 0 0];
            end

        end
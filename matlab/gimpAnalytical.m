        function Efield = getAnalyticalEfield(particle, decayLength, potential,surface_dz_dx,Bfield, ...
                VectorInterpolatorHandle,background_Z,background_amu,maxTemp_eV)
            B_local = VectorInterpolatorHandle(this,xyz,Bfield);
            B_unit = B_local/norm(B_local);
            
            surfaceDirection = [surface_dz_dx 0 1];
            surfaceDirection_unit = surfaceDirection/norm(surfaceDirection);
            
            normal_B_angle = acosd(dot(surfaceDirection_unit,B_unit))
            
            fd = -4E-11*normal_B_angle^6 + 1E-08*normal_B_angle^5 - 8E-07*normal_B_angle^4 + ...
                3E-05*normal_B_angle^3 - 0.0007*normal_B_angle^2 + 0.0051*normal_B_angle + 0.9899
            
            E_sheath = potential*fd/(2*decayLength)*exp(-particle.perpDistanceToSurface/(2*decayLength));
            
            larmor_radius = 1.44e-4*sqrt(background_amu(2)*maxTemp_eV)/(background_Z(2)*norm(B_local));
            E_magneticPresheath = potential*(1-fd)*exp(-particle.perpDistanceToSurface/(larmor_radius));
            
            Efield = (E_sheath+E_magneticPresheath)*surfaceDirection_unit;
            
        end
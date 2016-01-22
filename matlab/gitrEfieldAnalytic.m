function Efield = gitrEfieldAnalytic(this, decayLength, potential,surface_dz_dx,xyz,B_local, ...
                background_Z,background_amu,maxTemp_eV)
            
            B_unit = B_local/norm(B_local);
            
            surfaceDirection = [surface_dz_dx 0 1];
            surfaceDirection_unit = surfaceDirection/norm(surfaceDirection);
            
            normal_B_angle = acosd(dot(surfaceDirection_unit,B_unit));
            
            fd = -4.2682E-11*normal_B_angle^6 + 9.5856E-09*normal_B_angle^5 - ...
                8.2917E-07*normal_B_angle^4 + 3.3591E-05*normal_B_angle^3 - ...
                7.0040E-04*normal_B_angle^2 + 5.1220E-03*normal_B_angle + 9.8992E-01;

            
            E_sheath = potential*fd/(2*decayLength)*exp(-this.perpDistanceToSurface/(2*decayLength));
            
            larmor_radius = 1.44e-4*sqrt(background_amu(2)*maxTemp_eV)/(background_Z(2)*norm(B_local));
            E_magneticPresheath = potential*(1-fd)*exp(-this.perpDistanceToSurface/(larmor_radius));
            
            Efield = (E_sheath+E_magneticPresheath)*surfaceDirection_unit;
            
end
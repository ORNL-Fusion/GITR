function Efield = gitrEfieldAnalytic(this,xyz,Efield, decayLength, potential,surface_dz_dx,B_local, ...
                background_Z,background_amu,maxTemp_eV)
            norm_Blocal = realsqrt(B_local(1)*B_local(1) + B_local(2)*B_local(2) +B_local(3)*B_local(3));
            B_unit = B_local/norm_Blocal;
            
            surfaceDirection = [surface_dz_dx 0 -1];
            surfaceDirection_unit = surfaceDirection/realsqrt(surfaceDirection(1)*surfaceDirection(1)+surfaceDirection(2)*surfaceDirection(2)+surfaceDirection(3)*surfaceDirection(3));
            
            normal_B_angle = acosd(abs(surfaceDirection_unit(1)*B_unit(1)+surfaceDirection_unit(2)*B_unit(2)+surfaceDirection_unit(3)*B_unit(3)));
            
            fd = -4.2682E-11*normal_B_angle^6 + 9.5856E-09*normal_B_angle^5 - ...
                8.2917E-07*normal_B_angle^4 + 3.3591E-05*normal_B_angle^3 - ...
                7.0040E-04*normal_B_angle^2 + 5.1220E-03*normal_B_angle + 9.8992E-01;
            fd = 0.9398;
            decayLength =1.0513e-5;
            E_sheath = -potential(1)*(fd)/(2*decayLength)*exp(-this.perpDistanceToSurface/(2*decayLength));
            
            larmor_radius = 1.44e-4*realsqrt(background_amu(2)*maxTemp_eV(1)/2)/(background_Z(2)*norm_Blocal);
            larmor_radius = 3.64479e-4;
            E_magneticPresheath = -potential(1)*(1-fd)/(larmor_radius)*exp(-this.perpDistanceToSurface/(larmor_radius));


            ElectricPresheath = 0.0;
            Efield = (E_sheath+E_magneticPresheath)*surfaceDirection_unit + ElectricPresheath*[0 0 -1];


end
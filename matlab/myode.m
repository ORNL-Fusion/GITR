function dydt = myode(t,y,this,E,B,xyz,Einterpolator,Binterpolator, decayLength, potential,surface_dz_dx, ...
               background_Z,background_amu,maxTemp_eV)
    
    constants
    
    dydt = zeros(6,1);
    
    thisE = [0 0 0];
    thisB = [0 0 0];

    thisB = Binterpolator(this,xyz,B);
    thisE = Einterpolator(this,xyz,E, decayLength, potential,surface_dz_dx,thisB, ...
               background_Z,background_amu,maxTemp_eV);

    

        dydt(1:3) = y(4:6);    
        dydt(4:6) = (this.Z*Q/(this.amu*MI)) * (thisE + cross(y(4:6),thisB) );

    
end


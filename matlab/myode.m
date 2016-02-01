function dydt = myode(t,y,this,E,B,xyz,interpolatorHandle)
    
    constants
    
    dydt = zeros(6,1);
    
    thisE = [0 0 0];
    thisB = [0 0 0];


    thisE = interpolatorHandle(this,xyz,E);
    thisB = interpolatorHandle(this,xyz,B);
    

        dydt(1:3) = y(4:6);    
        dydt(4:6) = (this.Z*Q/(this.amu*MI)) * (thisE + cross(y(4:6),thisB) );

    
end


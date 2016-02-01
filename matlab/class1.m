classdef class1 < handle
    
    properties
        x
        y
        z
    end
    
    methods
        function Eclass = getE(this,xyz,Efield)
            nargin
            Eclass = interp1(xyz.x,Efield,this.x);
        end

    end
end
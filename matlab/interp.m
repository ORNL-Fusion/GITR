function E = interp(obj,xyz,Efield)
nargin
E = interpn(xyz.x,Efield,obj.x)
end
function perpDistanceToSurface = perpDist(xp,surface_dz_dx,surface_zIntercept)
perpDistanceToSurface = ( -surface_dz_dx*xp(1) + xp(3) + surface_zIntercept)/sqrt(surface_dz_dx^2+1);
end
function [value,isterminal,direction] = myEvent(t,y,surface_dz_dx,surface_zIntercept)
% Locate the time when distance from surface passes through zero in a 
% decreasing direction and stop integration.
value = ( -surface_dz_dx*y(1) + y(3) + surface_zIntercept)/sqrt(surface_dz_dx^2+1);     % Detect height = 0
isterminal = 1;   % Stop the integration
direction = -1;   % Negative direction only
end
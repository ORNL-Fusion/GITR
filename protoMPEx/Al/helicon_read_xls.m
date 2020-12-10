function [x y z voltage ne] = helicon_read_xls(surface_file)
% T = readtable('Beers_Helicon_3D_100kW_LowDensity_HeliconData_LayerMiddle9.xlsx');
% T = readtable('Beers_Helicon_3D_100kW_HighDensity_HeliconData_LayerMiddle6.xlsx');
%T = readtable(surface_file);
T=xlsread(surface_file);
voltage = T(8:end,27);

x = T(8:end,1);
y = T(8:end,2);
z = T(8:end,3);

ne = T(8:end,7);
scatter3(x,y,z,60,ne, 'filled')
colorbar
end

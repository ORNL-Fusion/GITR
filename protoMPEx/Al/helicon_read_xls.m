function [x y z voltage ne] = helicon_read_xls(surface_file)
% T = readtable('Beers_Helicon_3D_100kW_LowDensity_HeliconData_LayerMiddle9.xlsx');
% T = readtable('Beers_Helicon_3D_100kW_HighDensity_HeliconData_LayerMiddle6.xlsx');
T = readtable(surface_file);
voltage = T.IMABS;

x = T.x_m_;
y = T.y_m_;
z = T.z_m_;

ne = T.ne_1_m_3_;
scatter3(x,y,z,60,ne, 'filled')
colorbar
end

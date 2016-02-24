% clear variables
run('ERO_output/Wtile_test_Surface_step0.m')
x = dMdlX(1,:);
y = dMdlY(:,1);


% Deposition
NC = SurfCell(1,2,:,:);
NC = reshape(NC,iNX,iNY);

figure(1)
surf(x,y,zeros(iNX,iNY),log10(NC))

xlabel('x axis [m]')
ylabel('y axis [m]')
zlabel('z axis [m]')
title('Deposition [# Particles] Log10 Scale')
colormap jet
colorbar
view([0 90]);
%axis([-0.02 0.005 -0.01 0.006])


% Mean Energy
ME = SurfCell(1,5,:,:);
ME = reshape(ME,iNX,iNY);

figure(2)
surf(x,y,zeros(iNX,iNY),ME)

xlabel('x axis [m]')
ylabel('y axis [m]')
zlabel('z axis [m]')
title('Mean Energy [eV]')
colormap jet
colorbar
view([0 90]);
%axis([-0.02 0.005 -0.01 0.006])

% Mean Energy
MQ = SurfCell(1,6,:,:);
MQ = reshape(MQ,iNX,iNY);

figure(3)
surf(x,y,zeros(iNX,iNY),MQ)

xlabel('x axis [m]')
ylabel('y axis [m]')
zlabel('z axis [m]')
title('Mean Charge [#]')
colormap jet
colorbar
view([0 90]);
%axis([-0.02 0.005 -0.01 0.006])
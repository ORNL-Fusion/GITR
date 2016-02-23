%load('output/gitrEfield.mat')
Ez = Efield3D.z;
Ez = Ez(:,1,:);
Ez = reshape(Ez, [nXv,nZv]);

Ex = Efield3D.x;
Ex = Ex(:,1,:);
Ex = reshape(Ex, [nXv,nZv]);
figure(20)
surf(log10(abs(Ez)))
colorbar

% xslice = 0;
% yslice = 0.0;
% zslice = 0;
% figure(4)
% slice(xV_1D,yV_1D,zV_1D,log10(abs(Efield3D.z)),xslice,yslice,zslice)
%     xlabel('x axis')
%     ylabel('y axis')
%     zlabel('z axis')
%     title('GITR Efield')
% colorbar
% x1 =  (zV_1D - surface_zIntercept) / surface_dz_dx;
% y1 = yV_1D*0.0 + yslice;
% 
%      surface_dz_dx = 0.5774;
% surface_zIntercept = 0;
% surface_z = run_param.surf_x2D(1,:)*surface_dz_dx + surface_zIntercept;
%  hold on
%  plot3(run_param.surf_x2D(1,:),0*run_param.surf_z1D,run_param.surf_z1D)
%  hold off
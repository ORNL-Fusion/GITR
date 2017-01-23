
edges{1} = surf_y1D;
edges{2} = surf_z1D;

if plot3D
    [h,bins] = hist3(end_pos(:,2:3), 'Edges',edges);
    figure(3)
    set(gca,'ZDir','reverse');
    surf(surf_z1D,surf_y1D,surf_x2D,h)
    xlabel('z axis')
    ylabel('y axis')
    zlabel('x axis')
    title('Deposition [# particles]')
    axis equal
    colorbar
    ax.XDir = 'normal';
    hold on
    
    set(gca,'ZDir','reverse');
    view([180 90]);
    %scatter3(end_pos(:,3),end_pos(:,2),end_pos(:,1))
    hold on
end

%%%%%%%% Mean energy
EnergyBins = zeros(nY,nZ);
ChargeBins = zeros(nY,nZ);
for p=1:nP
    
    yIndex =  round((end_pos(p,2)- yMin)/(yMax - yMin)*nY);
    zIndex = round((end_pos(p,3)- zMin)/(zMax - zMin)*nZ);
    EnergyBins(yIndex,zIndex) = EnergyBins(yIndex,zIndex) + 0.5*impurity_amu*(end_pos(p,4)^2 + end_pos(p,5)^2 + end_pos(p,6)^2)*MI/Q;
    ChargeBins(yIndex,zIndex) = ChargeBins(yIndex,zIndex) + end_pos(p,7);
    
end
for j=1:nY
    for k=1:nZ
        
        if EnergyBins(j,k)
            EnergyBins(j,k) = EnergyBins(j,k)/h(j,k);
        end
        
        if ChargeBins(j,k)
            ChargeBins(j,k) = ChargeBins(j,k)/h(j,k);
        end
    end
    
end




figure(5)
surf(surf_z1D,surf_y1D,surf_x2D,EnergyBins)
xlabel('z axis')
ylabel('y axis')
zlabel('x axis')
title('Mean Energy [eV]')
axis equal
colorbar
ax.XDir = 'normal';
hold on

set(gca,'ZDir','reverse');
view([180 90]);



figure(6)
surf(surf_z1D,surf_y1D,surf_x2D,ChargeBins)
xlabel('z axis')
ylabel('y axis')
zlabel('x axis')
title('Mean Charge')
axis equal
colorbar
ax.XDir = 'normal';
hold on

set(gca,'ZDir','reverse');
view([180 90]);



figure(1)
hold on
if plotTracks
    for p=1:nP
        
        plot_p = p;
        plot3(zHistory(:,plot_p),yHistory(:,plot_p),xHistory(:,plot_p))
    end
end
plot_fac = 5e-3;
Bx_in = Bfield3D.x(1,1,1);
By_in = Bfield3D.y(1,1,1);
Bz_in = Bfield3D.z(1,1,1);

B_in = [Bx_in By_in Bz_in];
x_start = 0.00;
y_start = 0.00;
z_start = 0.00;
Bvec =  quiver3(z_start,y_start,x_start-0.01,Bz_in*plot_fac/norm(B_in),By_in*plot_fac/norm(B_in),Bx_in*plot_fac/norm(B_in),'Color','r');
Energy =  quiver3(z_start,y_start,x_start,energy_eV_z_start*plot_fac/norm(B_in),energy_eV_y_start*plot_fac/norm(B_in),energy_eV_x_start*plot_fac/norm(B_in),'Color','g');
legend([Bvec Energy],{'B-Field','Initial Energy'});
hold off

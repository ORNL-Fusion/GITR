density = zeros(volumeGridSize);

for p=1:nP
    
    for tt = 1:nT
        if impurityDensity(tt,p)
        [i, j, k] = ind2sub(volumeGridSize, impurityDensity(tt,p));
        density(i,j,k) = density(i,j,k) + dt;
        end
    end
end
xslice = 0;
yslice = -0.005;
zslice = 0;
figure(4)
slice(xV_1D,yV_1D,zV_1D,density,xslice,yslice,zslice)
colorbar
x1 =  (zV_1D - surface_zIntercept) / surface_dz_dx;
y1 = yV_1D*0.0 + yslice;
 hold on
 plot3(x1,y1,zV_1D)
 hold off
 
     xlabel('x axis')
    ylabel('y axis')
    zlabel('z axis')
    title('Surface')
    axis equal
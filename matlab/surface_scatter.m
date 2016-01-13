
part_pos = [];
for i=1:nP
    if particles(i).amu == 1e20
   part_pos = [part_pos; [particles(i).x,particles(i).y,particles(i).z]];
    end

end
edges{1} = surf_y1D;
edges{2} = surf_z1D;

if size(part_pos)

[h,bins] = hist3(part_pos(:,2:3), 'Edges',edges);
figure(3)
set(gca,'ZDir','reverse');
surf(surf_z1D,surf_y1D,surf_x2D,h)
            xlabel('z axis')
            ylabel('y axis')
            zlabel('x axis')
            title('Surface')
            axis equal
            colorbar
            ax.XDir = 'normal';
 hold on

 set(gca,'ZDir','reverse');
 view([75 10]);
 scatter3(part_pos(:,3),part_pos(:,2),part_pos(:,1))
 hold on 
 
 end
 figure(1) 
 hold on
         for p=1:nP

            plot_p = p;
            plot3(zHistory(:,plot_p),yHistory(:,plot_p),xHistory(:,plot_p))
         end
         plot_fac = 5e-3;
         B_in = [Bx_in By_in Bz_in];
        Bvec =  quiver3(z_start,y_start,x_start,Bz_in*plot_fac/norm(B_in),By_in*plot_fac/norm(B_in),Bx_in*plot_fac/norm(B_in),'Color','r');
        Energy =  quiver3(z_start,y_start,x_start,energy_eV_z_start*plot_fac/norm(B_in),energy_eV_y_start*plot_fac/norm(B_in),energy_eV_x_start*plot_fac/norm(B_in),'Color','g');
         legend([Bvec Energy],{'B-Field','Initial Energy'});
         hold off
         
         figure(4)
         a = zeros(nYv,nZv);
for i=1:nYv
    a(i,:) = impurityDensity(12,i,:);
end
surf(zV_1D,yV_1D,a)

colorbar
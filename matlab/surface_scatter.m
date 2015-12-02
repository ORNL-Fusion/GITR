
part_pos = zeros(nP,3);
for i=1:nP
   part_pos(i,:) = [particles(i).x,particles(i).y,particles(i).z]

end
edges{1} = y;
edges{2} = z;
[h,bins] = hist3(part_pos(:,2:3), 'Edges',edges)
figure(3)
set(gca,'ZDir','reverse');
surf(z,y,surfx,h)
            xlabel('z axis')
            ylabel('y axis')
            zlabel('x axis')
            title('Surface')
            axis equal
            colorbar
            ax.XDir = 'normal'
 hold on
 set(gca,'ZDir','reverse');
 view([75 10]);
 scatter3(part_pos(:,3),part_pos(:,2),part_pos(:,1))
 hold on 
 
         for p=1:nP

            plot_p = p;
            plot3(zHistory(1:n_steps,plot_p),yHistory(1:n_steps,plot_p),xHistory(1:n_steps,plot_p))
         end
         plot_fac = 5e-3;
         quiver3(0,0,-0.01,B_local(3)*plot_fac/norm(B_local),B_local(2)*plot_fac/norm(B_local),B_local(1)*plot_fac/norm(B_local),'Color','r')
         quiver3(0,0,-0.01,Ez_dat*plot_fac/norm(B_local),Ey_dat*plot_fac/norm(B_local),Ex_dat*plot_fac/norm(B_local),'Color','g')
         hold off
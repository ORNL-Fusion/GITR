% 
%                     magx = [e3(1) B_unit(1) e1(1) e2(1)];
%                     magy = [e3(2) B_unit(2) e1(2) e2(2)];
%                     magz = [e3(3) B_unit(3) e1(3) e2(3)];
%                     pause(0.1)
%                     figure(7)
%                     quiver3(xs(1),ys(1),zs(1),magx(1),magy(1),magz(1),'Color','r')
%                     hold on
%                     quiver3(xs(2),ys(2),zs(2),magx(2),magy(2),magz(2),'Color','g')
%                     hold on
%                     quiver3(xs(3),ys(3),zs(3),magx(3),magy(3),magz(3),'Color','b')
%                     hold on
%                     quiver3(xs(4),ys(4),zs(4),magx(4),magy(4),magz(4),'Color','m')
%                     xlabel('x axis')
%                     ylabel('y axis')
%                     zlabel('z axis')
%                     title('Vectors')
%                     legend('parallel', 'Bfield', 'perp1', 'perp2')
%                     hold off
%                     pause(0.1)
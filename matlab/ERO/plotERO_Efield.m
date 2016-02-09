figure (10)
xslice = 0; 
yslice = 0; 
zslice = 20;
slice(dX,dY,dZ,log10(abs(dE_z))+3,xslice,yslice,zslice)
colorbar
    xlabel('x axis')
    ylabel('y axis')
    zlabel('z axis')
    title('ERO Efield')
 caxis([-5 8])   
     
     